
# --------- data / batching ---------
gradient_accumulation_steps = 1
batch_size   = 128
block_size   = 48         # (your previous config had 48)
data_fraction = 1.0
ignore_tokens = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
no_event_token_rate = 5
t_min = 0.1
token_dropout = 0.0
mask_ties = True

# --------- model ---------
n_layer   = 12
n_head    = 12
n_embd    = 120
dropout   = 0.1
bias      = False
vocab_size = 1270         # (your config.py had 1270)

# --------- optimizer ---------
learning_rate = 6e-4
max_iters     = 100000
weight_decay  = 2e-1
beta1         = 0.9
beta2         = 0.99
grad_clip     = 1.0

# --------- lr schedule ---------
decay_lr      = True
warmup_iters  = 1000
lr_decay_iters = 100000
min_lr        = 6e-5

# --------- eval / logging cadence ---------
eval_interval = 250
eval_iters    = 25
log_interval  = 25

# --------- system ---------
dtype   = "float32"       # 'float32' | 'bfloat16' | 'float16' | 'float64'
compile = False
