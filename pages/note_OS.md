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

先打开 Clash，并设置 `System Proxy: on`、`Allow LAN: on`、端口 `7897`。

### Windows / WSL2

确认 Clash 端口可用：

```powershell
Test-NetConnection 127.0.0.1 -Port 7897
```

在 Windows 的 `%USERPROFILE%\.wslconfig` 中设置：

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
autoProxy=false
```

修改后运行 `wsl --shutdown`，再重新打开 WSL。

### WSL 代理

将以下内容加入 `~/.bashrc`，供 Git、pip、R 和 Hugging Face 等命令使用：

```bash
export http_proxy=http://127.0.0.1:7897
export https_proxy=$http_proxy
export all_proxy=$http_proxy
export no_proxy=localhost,127.0.0.1,::1
```

为 apt 单独设置永久代理：

```bash
sudo tee /etc/apt/apt.conf.d/95proxy >/dev/null <<'EOF'
Acquire::http::Proxy "http://127.0.0.1:7897";
Acquire::https::Proxy "http://127.0.0.1:7897";
EOF
```

测试连接：

```bash
source ~/.bashrc
curl -I --connect-timeout 8 https://github.com
sudo apt update
```

若突然无法联网，先重启 Clash，再在 Windows 运行 `wsl --shutdown`。若 apt 经代理访问清华源出现 `403`，改用 Ubuntu 官方源或 `https://cloud.r-project.org`。

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
