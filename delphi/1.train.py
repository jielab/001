#%%
import os
import time
import math
import pickle
from contextlib import nullcontext
import numpy as np
import torch
from model import Delphi, DelphiConfig
from utils import get_p2i, get_batch
import importlib.util

dataset   = "/mnt/d/data/ukb/phe/delphi"
device    = "cuda" # e.g. 'cpu', 'cuda', 'cuda:0'
init_from = "scratch" # 'scratch' or 'resume' or 'gpt2*'
seed      = 42
eval_only = False
wandb_log = False
out_dir   = "out"

# load config.py and copy its public names into globals()
config_path = os.path.abspath("config.py")
spec = importlib.util.spec_from_file_location("user_config", config_path)
cfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cfg)
globals().update({k: v for k, v in vars(cfg).items() if not k.startswith("_")})
config = {k: v for k, v in globals().items() if not k.startswith('_') and isinstance(v, (int, float, bool, str))}

if wandb_log:
    import wandb
    wandb.init(project = "delphi", name = "run" + str(time.time()), config = config)


#%%
os.makedirs(out_dir, exist_ok=True)
torch.manual_seed(seed)
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
device_type = 'cuda' if 'cuda' in device else 'cpu'
ptdtype = {'float32': torch.float32, 'float64': torch.float64, 'bfloat16': torch.bfloat16, 'float16': torch.float16}[dtype]
ctx = nullcontext() if device_type == 'cpu' else torch.amp.autocast(device_type=device_type, dtype=ptdtype)
torch.set_default_dtype(ptdtype)

if os.path.isabs(dataset): data_dir = dataset
else: data_dir = os.path.join('data', dataset) 

train_data = np.memmap(os.path.join(data_dir, 'train.bin'), dtype=np.uint32, mode='r').reshape(-1, 3)
val_data   = np.memmap(os.path.join(data_dir, 'val.bin'),   dtype=np.uint32, mode='r').reshape(-1, 3)
train_p2i = get_p2i(train_data)
val_p2i   = get_p2i(val_data)
if data_fraction < 1.0: train_p2i = train_p2i[:int(data_fraction * len(train_p2i))]

# init these up here, can override if init_from='resume' (i.e. from a checkpoint)
iter_num = 0
best_val_loss = 1e9
print(f"found vocab_size = {vocab_size}")
model_args = dict(
    n_layer=n_layer, n_head=n_head, n_embd=n_embd, block_size=block_size,
    bias=bias, vocab_size=vocab_size, dropout=dropout, token_dropout=token_dropout, t_min=t_min,
    mask_ties=mask_ties, ignore_tokens=ignore_tokens
)

if init_from == 'scratch':
    print("Initializing a new model from scratch")
    gptconf = DelphiConfig(**model_args)
    model = Delphi(gptconf)
elif init_from == 'resume':
    print(f"Resuming training from {out_dir}")
    ckpt_path = os.path.join(out_dir, 'ckpt.pt')
    checkpoint = torch.load(ckpt_path, map_location=device)
    checkpoint_model_args = checkpoint['model_args']
    for k in ['n_layer', 'n_head', 'n_embd', 'block_size', 'bias', 'vocab_size']: model_args[k] = checkpoint_model_args[k]
    gptconf = DelphiConfig(**model_args)
    model = Delphi(gptconf)
    state_dict = checkpoint['model']
    unwanted_prefix = '_orig_mod.'
    for k, v in list(state_dict.items()):
        if k.startswith(unwanted_prefix): state_dict[k[len(unwanted_prefix):]] = state_dict.pop(k)
    model.load_state_dict(state_dict)
    iter_num = checkpoint['iter_num']
    best_val_loss = checkpoint['best_val_loss']
model.to(device)

# initialize a GradScaler. If enabled=False scaler is a no-op
scaler = torch.cuda.amp.GradScaler(enabled=(dtype == 'float16'))
# optimizer
optimizer = model.configure_optimizers(weight_decay, learning_rate, (beta1, beta2), device_type)
if init_from == 'resume': optimizer.load_state_dict(checkpoint['optimizer'])
if compile:
    print("compiling the model... (takes a ~minute)")
    unoptimized_model = model
    model = torch.compile(model)  # requires PyTorch 2.0


#%%
@torch.no_grad()
def estimate_loss():
    out = {}
    model.eval()
    for split in ['train', 'val']:
        losses = torch.zeros(eval_iters, 2)
        data = train_data if split == 'train' else val_data
        p2i = train_p2i if split == 'train' else val_p2i
        for k in range(eval_iters):
            ix = torch.randint(len(p2i), (batch_size,))
            X, A, Y, B = get_batch(ix, data, p2i, block_size=block_size, device=device, select='left', no_event_token_rate=no_event_token_rate,  cut_batch=True)
            with ctx: logits, loss, _ = model(X, A, Y, B, validation_loss_mode=True)
            losses[k] = torch.stack([loss['loss_ce'], loss['loss_dt']])
        out[split] = losses.mean(0)
    model.train()
    return out

# learning rate decay scheduler (cosine with warmup)
def get_lr(it):
    if it < warmup_iters: return learning_rate * it / warmup_iters
    if it > lr_decay_iters: return min_lr
    decay_ratio = (it - warmup_iters) / (lr_decay_iters - warmup_iters)
    assert 0 <= decay_ratio <= 1
    coeff = 0.5 * (1.0 + math.cos(math.pi * decay_ratio))
    return min_lr + coeff * (learning_rate - min_lr)


# training loop
ix = torch.randint(len(train_p2i), (batch_size,))
X, A, Y, B = get_batch(ix, train_data, train_p2i, block_size=block_size, device=device, padding='random', lifestyle_augmentations=True, select='left', no_event_token_rate=no_event_token_rate)
t0 = time.time()
local_iter_num = 0

val_loss = None
while True:
    lr = get_lr(iter_num) if decay_lr else learning_rate
    for param_group in optimizer.param_groups: param_group['lr'] = lr

    if iter_num % eval_interval == 0 and iter_num > 0:
        losses = estimate_loss()
        if val_loss is None: val_loss_unpooled = losses['val']
        val_loss_unpooled = 0.1 * losses['val'] + 0.9 * val_loss_unpooled
        val_loss = val_loss_unpooled.sum().item()
        print(f"step {iter_num}: train loss {losses['train'].sum().item():.4f}, val loss {losses['val'].sum().item():.4f} ({val_loss:.4f})")
        if wandb_log: wandb.log({"iter": iter_num, "train/agg_loss": losses['train'].sum().item(), "val/loss": val_loss, "val/loss_ce": val_loss_unpooled[0].item(), "val/loss_dt": val_loss_unpooled[1].item()})

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            if iter_num > 0:
                checkpoint = {'model': model.state_dict(), 'optimizer': optimizer.state_dict(), 'model_args': model_args, 'iter_num': iter_num, 'best_val_loss': val_loss, 'config': config}
                print(f"saving checkpoint to {out_dir}")
                torch.save(checkpoint, os.path.join(out_dir, 'ckpt.pt'))

        if iter_num % 10_000 == 0:
            checkpoint = {'model': model.state_dict(), 'optimizer': optimizer.state_dict(), 'model_args': model_args, 'iter_num': iter_num, 'best_val_loss': best_val_loss, 'config': config}
            print(f"saving checkpoint to {out_dir}")
            torch.save(checkpoint, os.path.join(out_dir, f'ckpt_{iter_num}.pt'))

    if iter_num == 0 and eval_only: break
    for micro_step in range(gradient_accumulation_steps):
        with ctx: logits, loss, att = model(X, A, Y, B)
        ix = torch.randint(len(train_p2i), (batch_size,))
        X, A, Y, B = get_batch(ix, train_data, train_p2i, block_size=block_size, device=device, padding='random', lifestyle_augmentations=True, select='left', no_event_token_rate=no_event_token_rate, cut_batch=True)
        loss = loss['loss_ce'] + loss['loss_dt']
        scaler.scale(loss).backward()

    if grad_clip != 0.0: scaler.unscale_(optimizer); torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
    scaler.step(optimizer)
    scaler.update()
    optimizer.zero_grad(set_to_none=True)

    t1 = time.time()
    dt = t1 - t0
    t0 = t1
    if iter_num % log_interval == 0:
        lossf = loss.item()
        print(f"iter {iter_num}: loss {lossf:.4f}, time {dt*1000:.2f}ms")
        if wandb_log: wandb.log({"iter": iter_num,"train/loss": loss, "lr": lr, "weights": wandb.Histogram(model.transformer.wte.weight.cpu().detach().numpy()), "logits": wandb.Histogram(logits.cpu().detach().numpy())})

    iter_num += 1
    local_iter_num += 1
    if iter_num > max_iters: break
# %%
