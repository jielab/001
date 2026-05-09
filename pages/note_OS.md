# Windows, WSL, Python and local Transformer setup

> Public GitHub note: do not commit real passwords, private IPs, server names, access tokens, or institution-specific printer/network details. Use placeholders such as `<SERVER_IP>`, `<USER>`, and `<PROXY_PORT>`.

## 1. Common Windows and WSL commands

```powershell
# Restart Outlook
taskkill /f /im outlook.exe
start outlook

# Remove recycle-bin files on a drive; run CMD as administrator
rd /s /q D:\$RECYCLE.BIN

# WSL installation and status
wsl --list --online
wsl --install -d Ubuntu-24.04
wsl --set-default-version 2
wsl -l -v
wsl --status
wsl --version
```

Useful WSL paths:

```bash
# Windows D: drive in WSL
/mnt/d
```

If a shell window shows `press any key to continue` and the network stack appears broken, try running this in administrator CMD:

```cmd
netsh winsock reset
```

## 2. Background downloads

```bash
sudo apt install -y aria2 screen
screen -dmS download aria2c -x 4 -i url.txt --log-level=info --log=download.log
screen -ls
screen -S download -X quit
```

Split a large command list into smaller scripts:

```bash
awk '{cnt=int(NR/100); print $0 > "download"cnt".sh"}' commands.txt
```

## 3. New Windows/WSL network setup: Clash, apt, R, Git, and SSH

Goal on a new Windows computer:

```text
sudo apt update
R package installation from CRAN/GitHub
GitHub/Hugging Face downloads
ssh <USER>@<SERVER_IP>
```

The clean rule is:

```text
WSL autoProxy: false
apt: no apt proxy file
shell proxy: default off; manually run proxy_on only when needed
```

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

The examples below use port `7897`. If your Clash port is different, replace `7897` with your own port.

### 3.2 WSL2 mirrored networking: keep autoProxy=false

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

Reason:

- `networkingMode=mirrored` helps WSL follow the Windows network and VPN routes.
- `dnsTunneling=true` helps DNS work under VPN or complex networks.
- `autoProxy=false` is important. If it is `true`, Windows proxy settings may be automatically injected into WSL and cause confusing `http_proxy` / `https_proxy` variables in every new terminal.

After reopening Ubuntu, check:

```bash
echo "=== shell proxy ==="
env | grep -Ei '^(http_proxy|https_proxy|all_proxy|HTTP_PROXY|HTTPS_PROXY|ALL_PROXY)=' || echo "No shell proxy"

echo "=== apt proxy ==="
apt-config dump | grep -Ei 'proxy|forceipv4|retries|pipeline' || echo "No apt proxy"
```

The expected result is:

```text
No shell proxy
No apt proxy
```

### 3.3 apt should normally use direct connection

Do **not** create `/etc/apt/apt.conf.d/95proxy` by default.

If this file exists from an old setup, remove it:

```bash
sudo rm -f /etc/apt/apt.conf.d/95proxy
```

Then test:

```bash
sudo apt update
```

If `sudo apt update` works, keep apt direct. Do not set an apt proxy.

Check again:

```bash
apt-config dump | grep -Ei 'proxy|forceipv4|retries|pipeline' || echo "No apt proxy"
```

Expected:

```text
No apt proxy
```

### 3.4 Add manual proxy_on / proxy_off to ~/.bashrc

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

Default terminal should be proxy off:

```bash
env | grep -Ei '^(http_proxy|https_proxy|all_proxy|HTTP_PROXY|HTTPS_PROXY|ALL_PROXY)=' || echo "No shell proxy"
```

Use proxy only when needed:

```bash
proxy_on
# run git / pip / R / hf download commands
proxy_off
```

### 3.5 When to use proxy_on

Usually **do not** use proxy for:

```bash
sudo apt update
sudo apt install xxx
ssh <USER>@<SERVER_IP>
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

### 3.6 SSH and internal servers

Test an internal server directly. Do not send SSH through Clash unless you know you need it.

```bash
proxy_off
ip route get <SERVER_IP>
ping -c 4 <SERVER_IP>
nc -vz <SERVER_IP> 22
ssh <USER>@<SERVER_IP>
```

If `ping` and SSH work but this fails:

```bash
curl -I http://<SERVER_IP>
```

that usually only means the server does not run an HTTP service on port 80. It does not mean SSH is broken.

### 3.7 R package installation

First try direct R package installation. If it fails because of network problems, start R after `proxy_on`.

```bash
proxy_on
R
```

Inside R:

```r
Sys.getenv(c("http_proxy", "https_proxy", "all_proxy"))
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("r-lib/lintr")
```

After installation:

```bash
proxy_off
```

### 3.8 Quick diagnosis commands

```bash
echo "=== WSL config ==="
winhome=$(powershell.exe -NoProfile -Command '[Environment]::GetFolderPath("UserProfile")' | tr -d '\r')
cat "$(wslpath "$winhome")/.wslconfig"

echo "=== shell proxy ==="
env | grep -Ei '^(http_proxy|https_proxy|all_proxy|HTTP_PROXY|HTTPS_PROXY|ALL_PROXY)=' || echo "No shell proxy"

echo "=== apt proxy ==="
apt-config dump | grep -Ei 'proxy|forceipv4|retries|pipeline' || echo "No apt proxy"

echo "=== bashrc proxy lines ==="
grep -nEi 'proxy_on|proxy_off|export .*proxy|http_proxy|https_proxy|all_proxy' ~/.bashrc

echo "=== git proxy ==="
git config --list --show-origin | grep -i proxy || echo "No git proxy"

echo "=== apt test ==="
sudo apt update
```

Expected clean status:

```text
autoProxy=false
No shell proxy
No apt proxy
No git proxy
sudo apt update works
```

## 4. Local Python and Transformer environment

Recommended VS Code extensions:

```text
WSL
Python
Jupyter
```

Check Python paths:

```bash
which python
python --version
python -c "import sys; print(sys.executable)"
```

On Windows CMD/PowerShell:

```powershell
where python
where pip
python -m pip -V
```

## 5. PyTorch and common AI packages

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

## 6. Hugging Face downloads

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

Login:

```bash
hf auth login
```

Test with a small file first:

```bash
mkdir -p /mnt/d/models/qwen

hf download Qwen/Qwen3-8B config.json \
  --local-dir /mnt/d/models/qwen/Qwen3-8B-test
```

Download the full Qwen model:

```bash
hf download Qwen/Qwen3-8B \
  --local-dir /mnt/d/models/qwen/Qwen3-8B
```

If a download stops, run the same `hf download` command again.

Optional Git LFS clone:

```bash
git lfs install
git clone https://huggingface.co/Qwen/Qwen3-8B
```

If `git clone` fails with a timeout, use `hf download` instead.
