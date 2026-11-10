# Windows, WSL, Python and local Transformer setup


## 1. 常识和常用命令

```powershell
# Restart Outlook
taskkill /f /im outlook.exe

sudo mkdir -p /mnt/h; sudo mount -t drvfs H: /mnt/h
aria2c -c -i ckb.url -d raw/ -x 8 -s 8 -j 8 \
  --all-proxy="" --http-proxy="" --https-proxy="" --no-proxy="ftp.sra.ebi.ac.uk" \
  --auto-file-renaming=false --allow-overwrite=false --always-resume=true --connect-timeout=30 \
  --timeout=120 --max-tries=20 --retry-wait=30 --auto-save-interval=30 --file-allocation=none
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

# 当出现 `press any key to continue` :
netsh winsock reset

```
---

## 2. HPC
使用指南：https://hpc.sustech.edu.cn/ref/HPMS_UserGuide.pdf

```bash
【太乙】: ssh sph-huangj@172.18.6.175[6|8]
【启明】: ssh -p 18188 sph-huangj@172.18.6.10[1|2] 
	shownodes; queueinfo
	CPU一般用 33 [155节点, 40核, 192G, 720h]，中等用38（36节点, 64核, 512G, 720h），大内存需求用1t81c和v3-6t (3节点, 144核, 6T, 336h)
	GPU一般用 4v100-16-e5 (2节点，4'CARD，16G|64G，336h)

bsub -J copy -q 38 -n 1 -W 120:00 -o copy.LOG -e copy.ERR 'rsync -a --partial --info=progress2 /scratch/2026-11-01/sph-huangj/XXX/ /work/sph-huangj/data/XXX/'

du -h --max-depth=1🏮; du -h --max-depth=2; mmlsquota -g sph-huangj --block-size auto
后台运行： nohup ./assoc.sum.sh & 之后 ps aux | grep ?.sh 之后 kill

慧园5栋一层科学与计算中心机房
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

`%USERPROFILE%\.codex\.env`

```dotenv
HTTP_PROXY=http://127.0.0.1:7897
HTTPS_PROXY=http://127.0.0.1:7897
NO_PROXY=localhost,127.0.0.1,::1
```

```powershell
curl.exe -q -I --proxy http://127.0.0.1:7897 --noproxy localhost `
  --connect-timeout 8 --max-time 20 https://auth.openai.com/oauth/token
```

### 3.2 WSL

`%USERPROFILE%\.wslconfig`

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
autoProxy=true
```

---

## 4. PyTorch and common AI packages

```bash
conda env list
conda create -n ai python=3.12
conda activate ai

pip uninstall -y torch torchvision torchaudio
pip install --upgrade pip
pip install torch torchvision torchaudio
pip install --upgrade numpy tqdm transformers pandas requests openpyxl bitsandbytes accelerate datasets peft evaluate scikit-learn protobuf sentencepiece huggingface_hub tabpfn pytorch-tabular[all]

python -c "import torch; print(torch.cuda.is_available()); print(torch.__version__); print(torch.version.cuda); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU only')"
```
---

## 5. Hugging Face downloads

```bash
hf auth login
hf download Qwen/Qwen3-8B --local-dir /mnt/d/models/qwen/Qwen3-8B

git lfs install
git clone https://huggingface.co/Qwen/Qwen3-8B

scripts/f/00hf_download.py 
```
---