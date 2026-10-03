# Windows, WSL, Python and local Transformer setup


## 1. 常识和常用命令

```powershell
# Restart Outlook
taskkill /f /im outlook.exe
start outlook

sudo mkdir -p /mnt/h; sudo mount -t drvfs H: /mnt/h
aria2c -c -i ckb.url -d raw/ -x 4 -s 4 -j 2
```

```
rm -rf .git
git init -b main
git add -A
git commit -m "Initial commit"
git rev-list --count HEAD # 🏮1
# git config --global user.name "Jie Huang"
# git config --global user.email "jielab@users.noreply.github.com"
# gh auth login
# gh auth setup-git
git remote add origin https://github.com/jielab/analysis.git
git push -u --force origin main

# WSL installation and status
wsl --list --online
wsl --install -d Ubuntu-24.04; wsl --set-default-version 2
wsl -l -v; wsl --status; wsl --version

# 安装Anaconda后
conda init powershell; Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

# 备份文件到J盘
rsync -a --delete --info=progress2 --exclude='/analysis/' --exclude='/System Volume Information/' --exclude='/$RECYCLE.BIN/' /mnt/d/ '/mnt/j/黄捷备份文件/'

# 当出现 `press any key to continue` :
netsh winsock reset

```
---

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
---

## 3. Clean Windows/WSL network setup

VPN 使用 `http://127.0.0.1:7897`，关闭随机端口 📍

### 3.1 Windows 桌面 ChatGPT

`%USERPROFILE%\.codex\.env`：合并以下项，保留其他内容；修改后完全退出并重开桌面应用。

```dotenv
HTTP_PROXY=http://127.0.0.1:7897
HTTPS_PROXY=http://127.0.0.1:7897
NO_PROXY=localhost,127.0.0.1,::1
```

PowerShell 测试 Windows 经过该代理访问认证地址：

```powershell
curl.exe -q -I --proxy http://127.0.0.1:7897 --noproxy localhost `
  --connect-timeout 8 --max-time 20 https://auth.openai.com/oauth/token
```

最终返回 `405 Method Not Allowed` 表示 HEAD 方法不被接受，但已收到认证服务响应；不等于登录成功。连接本地代理失败，先检查 VPN 是否运行及端口是否为 `7897`。

### 3.2 WSL

`%USERPROFILE%\.wslconfig`：合并到同一个 `[wsl2]` 段，保留其他设置。

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
autoProxy=true
```

**仅修改 `.wslconfig` 后**，保存 WSL 工作，再在 PowerShell 执行 `wsl --shutdown`。本次只改代理端口，不需要关闭 WSL。

保留现有 `~/.wsl_proxy.sh` 的代理函数，只确认端口设置：

```bash
export WSL_PROXY_HOST=127.0.0.1
export WSL_PROXY_PORT=7897
```

`~/.bashrc` 只保留一次加载，不重复定义代理：

```bash
[ -f ~/.wsl_proxy.sh ] && . ~/.wsl_proxy.sh
```

当前 Bash 重新加载并检查：

```bash
source ~/.wsl_proxy.sh
proxy_status
curl -I --connect-timeout 8 --max-time 20 https://github.com
```

### 3.3 aria2 / apt / 停用代理

`~/.aria2/aria2.conf`：

```ini
all-proxy=http://127.0.0.1:7897
no-proxy=localhost,127.0.0.1,::1
auto-file-renaming=false
```

`/etc/apt/apt.conf.d/95proxy`：

```text
Acquire::http::Proxy "http://127.0.0.1:7897";
Acquire::https::Proxy "http://127.0.0.1:7897";
```

不再使用代理时：Bash 执行 `proxy_off`；桌面端注释 `.codex\.env` 的代理项并重开；aria2 注释 `all-proxy`；apt 将 `95proxy` 移出 `apt.conf.d`。新 Bash 也不使用代理时，停用 `.bashrc` 中上述加载行。

**`proxy_off` 只影响当前 shell，不会取消桌面端、aria2 或 apt 的固定代理。** 配置修改不影响已运行的程序；下载任务需正常停止后重启，保留续传文件。

---

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
---

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
---

## 6. Background downloads

```bash
sudo apt install -y aria2 screen
screen -dmS download aria2c -x 4 -i url.txt --log-level=info --log=download.log
screen -ls
screen -S download -X quit
```
