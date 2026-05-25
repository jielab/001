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
使用指南：https://hpc.sustech.edu.cn/ref/HPMS_UserGuide.pdf

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
bjobs -wu sph-huangj | awk 'NR>1 {print $6}' | awk -F "*" '{print $NF}' | tr '\n' ' ' | xargs lsload
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


## 3. Clean Windows/WSL network setup: Clash, apt, Git, pip, R, and Hugging Face

目标：**平时不要记 `proxy_on`，也不要记 `sudo -E`。**

最终效果：

```bash
sudo apt update
git clone https://github.com/xxx/xxx.git
pip install xxx
hf download Qwen/Qwen3-8B --local-dir /mnt/d/models/qwen/Qwen3-8B
```

都尽量直接运行。唯一需要记住的是：**先打开 Clash，并确认端口是 7897。**

---

### 3.1 Windows Clash 设置

In Clash for Windows or Clash Verge:

```text
System Proxy: on
Allow LAN: on
Port: 7897
```

Check the port in Windows PowerShell:

```powershell
Test-NetConnection 127.0.0.1 -Port 7897
netstat -ano | findstr :7897
```

Expected:

```text
TcpTestSucceeded : True
0.0.0.0:7897    LISTENING
```

---

### 3.2 WSL2 网络模式设置

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

Open Ubuntu/WSL again.

---

### 3.3 One-time clean setup inside WSL

Run this once in WSL. It does four things:

1. writes a permanent apt proxy, so `sudo apt update` does not need `sudo -E`;
2. replaces TUNA 清华源 with more stable official Ubuntu / CRAN sources, avoiding `403 Forbidden [IP: 127.0.0.1 7897]`;
3. creates an automatic WSL shell proxy script;
4. loads that script automatically from `~/.bashrc`.

```bash
# ------------------------------------------------------------
# 0. Basic check: Windows Clash should be reachable from WSL
# ------------------------------------------------------------
curl -I -x http://127.0.0.1:7897 --connect-timeout 8 https://github.com

# ------------------------------------------------------------
# 1. Permanent apt proxy: sudo apt will use Clash automatically
# ------------------------------------------------------------
sudo tee /etc/apt/apt.conf.d/95proxy >/dev/null <<'EOF'
Acquire::http::Proxy "http://127.0.0.1:7897";
Acquire::https::Proxy "http://127.0.0.1:7897";
Acquire::ForceIPv4 "true";
Acquire::Retries "3";
Acquire::http::Timeout "30";
Acquire::https::Timeout "30";
EOF

# ------------------------------------------------------------
# 2. Backup apt source files
# ------------------------------------------------------------
stamp=$(date +%Y%m%d_%H%M%S)
sudo cp -a /etc/apt/sources.list /etc/apt/sources.list.bak.$stamp 2>/dev/null || true
sudo cp -a /etc/apt/sources.list.d /etc/apt/sources.list.d.bak.$stamp

# ------------------------------------------------------------
# 3. Replace TUNA Ubuntu and CRAN sources
#    This avoids TUNA 403 errors when apt goes through Clash.
# ------------------------------------------------------------
sudo grep -rl "mirrors.tuna.tsinghua.edu.cn/ubuntu" /etc/apt/sources.list /etc/apt/sources.list.d 2>/dev/null | \
  xargs -r sudo sed -i \
  's|http://mirrors.tuna.tsinghua.edu.cn/ubuntu|http://archive.ubuntu.com/ubuntu|g; s|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|http://archive.ubuntu.com/ubuntu|g'

sudo grep -rl "mirrors.tuna.tsinghua.edu.cn/CRAN" /etc/apt/sources.list /etc/apt/sources.list.d 2>/dev/null | \
  xargs -r sudo sed -i \
  's|https://mirrors.tuna.tsinghua.edu.cn/CRAN|https://cloud.r-project.org|g; s|http://mirrors.tuna.tsinghua.edu.cn/CRAN|https://cloud.r-project.org|g'

# ------------------------------------------------------------
# 4. Automatic WSL shell proxy for git, pip, R, Hugging Face, etc.
#    If Clash port is open, proxy is ON automatically.
#    If Clash port is closed, proxy is unset automatically.
# ------------------------------------------------------------
cat > ~/.wsl_proxy.sh <<'EOF'
# ------------------------------------------------------------
# Auto proxy for WSL + Clash
# Default behavior:
#   - if 127.0.0.1:7897 is reachable, proxy is enabled automatically
#   - if not reachable, proxy variables are removed automatically
# This affects normal shell commands: git, curl, pip, R, hf, etc.
# apt is handled separately by /etc/apt/apt.conf.d/95proxy
# ------------------------------------------------------------

export WSL_PROXY_HOST=127.0.0.1
export WSL_PROXY_PORT=7897
export WSL_PROXY_URL=http://${WSL_PROXY_HOST}:${WSL_PROXY_PORT}
export no_proxy=localhost,127.0.0.1,::1
export NO_PROXY=localhost,127.0.0.1,::1

_proxy_set() {
  export http_proxy=${WSL_PROXY_URL}
  export https_proxy=${WSL_PROXY_URL}
  export all_proxy=${WSL_PROXY_URL}
  export HTTP_PROXY=${WSL_PROXY_URL}
  export HTTPS_PROXY=${WSL_PROXY_URL}
  export ALL_PROXY=${WSL_PROXY_URL}
}

_proxy_unset() {
  unset http_proxy https_proxy all_proxy
  unset HTTP_PROXY HTTPS_PROXY ALL_PROXY
}

proxy_on() {
  _proxy_set
  echo "Proxy on: ${WSL_PROXY_URL}"
}

proxy_off() {
  _proxy_unset
  echo "Proxy off"
}

proxy_status() {
  env | grep -Ei '^(http|https|all|no)_proxy=' | sort
}

wsl_proxy_auto() {
  if timeout 1 bash -c ":</dev/tcp/${WSL_PROXY_HOST}/${WSL_PROXY_PORT}" 2>/dev/null; then
    _proxy_set
  else
    _proxy_unset
  fi
}

wsl_proxy_auto
EOF

grep -q "wsl_proxy.sh" ~/.bashrc || cat >> ~/.bashrc <<'EOF'

# Auto proxy for WSL + Clash
[ -f ~/.wsl_proxy.sh ] && . ~/.wsl_proxy.sh
EOF

source ~/.bashrc

# ------------------------------------------------------------
# 5. Test
# ------------------------------------------------------------
curl -I --connect-timeout 8 https://github.com
curl -I --connect-timeout 8 https://r2u.stat.illinois.edu
sudo apt clean
sudo apt update
```

After this setup, normal use should be:

```bash
sudo apt update
sudo apt install -y git curl build-essential
pip install xxx
git clone https://github.com/xxx/xxx.git
```

No `proxy_on` and no `sudo -E` should be needed for `apt`.

---

### 3.4 Why `sudo -E` was needed before

`proxy_on` only exports proxy variables in the current user shell:

```bash
http_proxy=http://127.0.0.1:7897
https_proxy=http://127.0.0.1:7897
all_proxy=http://127.0.0.1:7897
```

However, `sudo` normally removes most user environment variables for safety. Therefore:

```bash
sudo apt update
```

may ignore the proxy variables, while:

```bash
sudo -E apt update
```

keeps them.

The clean solution is **not** to remember `sudo -E`. The clean solution is to give apt its own proxy configuration:

```text
/etc/apt/apt.conf.d/95proxy
```

Then:

```bash
sudo apt update
```

works directly.

---

### 3.5 Quick diagnosis

#### A. Check Windows Clash

Run in Windows PowerShell:

```powershell
Test-NetConnection 127.0.0.1 -Port 7897
netstat -ano | findstr :7897
```

#### B. Check WSL shell proxy

Run in WSL:

```bash
proxy_status
curl -I --connect-timeout 8 https://github.com
```

#### C. Check apt proxy

Run in WSL:

```bash
cat /etc/apt/apt.conf.d/95proxy
sudo apt update
```

#### D. If apt says `403 Forbidden [IP: 127.0.0.1 7897]`

This usually means the proxy works, but the apt mirror rejects the request. Replace that mirror.

Check:

```bash
grep -R "mirrors.tuna.tsinghua.edu.cn" /etc/apt/sources.list /etc/apt/sources.list.d 2>/dev/null
```

Replace TUNA sources:

```bash
sudo grep -rl "mirrors.tuna.tsinghua.edu.cn/ubuntu" /etc/apt/sources.list /etc/apt/sources.list.d 2>/dev/null | \
  xargs -r sudo sed -i \
  's|http://mirrors.tuna.tsinghua.edu.cn/ubuntu|http://archive.ubuntu.com/ubuntu|g; s|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|http://archive.ubuntu.com/ubuntu|g'

sudo grep -rl "mirrors.tuna.tsinghua.edu.cn/CRAN" /etc/apt/sources.list /etc/apt/sources.list.d 2>/dev/null | \
  xargs -r sudo sed -i \
  's|https://mirrors.tuna.tsinghua.edu.cn/CRAN|https://cloud.r-project.org|g; s|http://mirrors.tuna.tsinghua.edu.cn/CRAN|https://cloud.r-project.org|g'

sudo apt clean
sudo apt update
```

#### E. If everything suddenly stops working

Restart Clash, then restart WSL:

```powershell
wsl --shutdown
```

Open WSL again and test:

```bash
proxy_status
curl -I --connect-timeout 8 https://github.com
sudo apt update
```

---

### 3.6 Optional: remove the apt proxy

Only do this if Clash is not needed anymore or if you are on a network where direct apt is stable.

```bash
sudo rm -f /etc/apt/apt.conf.d/95proxy
sudo apt update
```

---

### 3.7 codeX proxy setup

If Windows Codex App shows `Reconnecting...`, `timeout waiting for child process to exit`, or `无法重新安装 Codex 依赖项`, first check whether Windows can reach OpenAI through Clash.

Keep Clash running, and check that the local proxy port is open:

```powershell
Test-NetConnection 127.0.0.1 -Port 7897
netstat -ano | findstr :7897
```

Test OpenAI through the Clash proxy explicitly:

```powershell
curl.exe -I -L -x http://127.0.0.1:7897 https://api.openai.com/v1/models
curl.exe -I -L -x http://127.0.0.1:7897 https://chatgpt.com
```

Expected results:

```text
api.openai.com/v1/models: HTTP/1.1 401 Unauthorized is OK.
chatgpt.com: HTTP 403 Cloudflare challenge may appear in curl and is not the key problem.
```

If direct `curl` fails but `curl -x` works, set proxy variables for the current PowerShell session:

```powershell
$env:HTTP_PROXY="http://127.0.0.1:7897"
$env:HTTPS_PROXY="http://127.0.0.1:7897"
$env:ALL_PROXY="http://127.0.0.1:7897"
$env:NO_PROXY="localhost,127.0.0.1,::1"

curl.exe -I -L https://api.openai.com/v1/models
```

If this returns `401 Unauthorized`, write the proxy variables permanently to the Windows user environment:

```powershell
setx HTTP_PROXY "http://127.0.0.1:7897"
setx HTTPS_PROXY "http://127.0.0.1:7897"
setx ALL_PROXY "http://127.0.0.1:7897"
setx NO_PROXY "localhost,127.0.0.1,::1"
```

Then close all Codex processes, open a new PowerShell, and test again:

```powershell
Get-Process | Where-Object {$_.ProcessName -match "codex|openai|chatgpt"} | Stop-Process -Force
curl.exe -I -L https://api.openai.com/v1/models
```

After `curl` returns `401 Unauthorized`, reopen Windows Codex App and run:

```text
Settings -> Diagnose
Settings -> Reinstall Codex dependencies
```

Optional notes:

```text
1. PowerShell setx does not affect already-open Codex App or PowerShell windows.
2. Keep Clash running before opening Codex App.
3. netsh winhttp set proxy needs administrator PowerShell, but Codex may still rely on HTTP_PROXY/HTTPS_PROXY instead.
4. If Windows networking becomes abnormal, remove the WinHTTP proxy with: netsh winhttp reset proxy
5. For debugging, switch Clash from Rule to Global mode in the Proxies page, then switch back after Codex works.
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

After the clean WSL proxy setup above, Hugging Face commands should work without manually typing `proxy_on`.

Install the Hugging Face CLI with `pipx`:

```bash
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

If Hugging Face still fails, check proxy status first:

```bash
proxy_status
curl -I --connect-timeout 8 https://huggingface.co
```

再不行，就用 scripts/f/00hf_download.py 


## 6. Background downloads

```bash
sudo apt install -y aria2 screen
screen -dmS download aria2c -x 4 -i url.txt --log-level=info --log=download.log
screen -ls
screen -S download -X quit
```