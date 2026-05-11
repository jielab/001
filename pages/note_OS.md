# Windows, WSL, Python and local Transformer setup


## 1. 常识和常用命令

```powershell
# Restart Outlook
taskkill /f /im outlook.exe
start outlook

# Remove recycle-bin files on a drive; run CMD as administrator
rd /s /q D:\$RECYCLE.BIN

# WSL installation and status
wsl --list --online
wsl --install -d Ubuntu-24.04; wsl --set-default-version 2
wsl -l -v; wsl --status; wsl --version

# 安装Anaconda后
conda init powershell; Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

If a shell window shows `press any key to continue` and the network stack appears broken, try running this in administrator CMD:

```cmd
netsh winsock reset
```

## 2. HPC

```bash
【太乙】: ssh sph-huangj@172.18.6.178
【启明】: ssh -p 18188 sph-huangj@172.18.6.10 
```

后台运行
```
nohup ./assoc.sum.sh & 之后 ps aux | grep ?.sh 之后 kill
```

硬盘额度
```
du -h --max-depth=2; mmlsquota -g sph-huangj --block-size auto
```

bsub
```
queueinfo -gpu -cpu; module avail
```

创园301🖨
```
从富士官网(https://m3support-fb.fujifilm-fb.com.cn/driver_downloads/www/)搜索 ApeosPort C2060 下载安装驱动程序 
 👉“设备类型” 选TCP/IP 👉 打印机IP为 10.20.40.6
 ```
 
创园204🖨 
```
连接 LINK_7204无线网，密码是???2025??04，然后下载安装驱动程序(https://www.canon.com.cn/supports/download/simsdetail/0101228601.html?modelId=1524&channel=4)
```


## 3. New Windows/WSL network setup: Clash, apt, R, Git, and SSH


### 3.1 Clash settings in Windows

In Clash for Windows or Clash Verge:

```text
System Proxy: on
Allow LAN: on
Port: 7897
```

Check the port in Windows PowerShell:

```powershell
netstat -ano | findstr :7897
```

Create or edit the WSL configuration file in Windows PowerShell:

```powershell
notepad $env:USERPROFILE\.wslconfig
```

Recommended content:

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
autoProxy=false
```

Then restart WSL:

```powershell
wsl --shutdown
```
```bash
sudo apt update
```


### 3.3 Add manual proxy_on / proxy_off to ~/.bashrc

Add only the functions below to the end of `~/.bashrc`. Do **not** put `export http_proxy=...` outside the function. Do **not** add a standalone `proxy_on` line.

```bash
cat >> ~/.bashrc <<'BASHRC_EOF'

# ------------------------------------------------------------
# Manual proxy control for WSL + Clash
# Default: proxy is OFF.
# Use proxy_on only when GitHub, Hugging Face, pip, or R downloads fail.
# ------------------------------------------------------------
proxy_on() {
  export http_proxy=http://127.0.0.1:7897
  export https_proxy=http://127.0.0.1:7897
  export all_proxy=http://127.0.0.1:7897
  export HTTP_PROXY=http://127.0.0.1:7897
  export HTTPS_PROXY=http://127.0.0.1:7897
  export ALL_PROXY=http://127.0.0.1:7897
  echo "Proxy on: 127.0.0.1:7897"
}

proxy_off() {
  unset http_proxy https_proxy all_proxy
  unset HTTP_PROXY HTTPS_PROXY ALL_PROXY
  echo "Proxy off"
}
BASHRC_EOF

source ~/.bashrc
```

Use `proxy_on` when direct download fails for:

```bash
git clone https://github.com/xxx/xxx.git
git clone https://huggingface.co/Qwen/Qwen3-8B
pip install xxx
pipx install "huggingface_hub[cli]"
hf download Qwen/Qwen3-8B --local-dir /mnt/d/models/qwen/Qwen3-8B
```

Test Clash proxy manually:

```bash
proxy_on
curl -I --connect-timeout 8 https://huggingface.co
curl -I --connect-timeout 8 https://github.com
proxy_off
```

## 4. PyTorch and common AI packages

Windows Powershell 上安装，可以不用conda
Example conda environment:

```bash
conda env list
conda create -n ai python=3.12
conda activate ai
```

Install PyTorch and common packages. Choose the PyTorch command appropriate for your CUDA/driver version from the official PyTorch installation page.

```bash
pip uninstall -y torch torchvision torchaudio
pip install --upgrade pip
pip install torch torchvision torchaudio
pip install --upgrade numpy tqdm transformers pandas requests openpyxl bitsandbytes accelerate datasets peft evaluate scikit-learn protobuf sentencepiece huggingface_hub tabpfn pytorch-tabular[all]
```

Check GPU availability:

```bash
python -c "import torch; print(torch.cuda.is_available()); print(torch.__version__); print(torch.version.cuda); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU only')"
```

## 5. Hugging Face downloads

For Hugging Face, direct connection may fail. Use `proxy_on` first.

Install the Hugging Face CLI with `pipx`:

```bash
proxy_on
sudo apt update
sudo apt install -y pipx python3.12-venv git git-lfs
pipx ensurepath
source ~/.bashrc
pipx install "huggingface_hub[cli]"
```

Check:

```bash
hf --help
hf version
hf env
```

Login and Download:

```bash
hf auth login
hf download Qwen/Qwen3-8B --local-dir /mnt/d/models/qwen/Qwen3-8B
```

Optional Git LFS clone:

```bash
git lfs install
git clone https://huggingface.co/Qwen/Qwen3-8B
```

再不行，就用 scripts/f/00hf_download.py 


## 6. Background downloads

```bash
sudo apt install -y aria2 screen
screen -dmS download aria2c -x 4 -i url.txt --log-level=info --log=download.log
screen -ls
screen -S download -X quit
```