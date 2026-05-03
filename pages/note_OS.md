# Windows, WSL, Python and local Transformer setup

> Public GitHub note: do not commit real passwords, private IPs, server names, access tokens, or institution-specific printer/network details. Use placeholders such as `<SERVER_IP>` or `<PROXY_IP>`.

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

## 3. WSL proxy through Clash Verge

In Clash Verge:

```text
System proxy: on
Allow LAN: on
Port: 7897
```

Check the listening port in Windows PowerShell:

```powershell
netstat -ano | findstr :7897
```

Find the working proxy address from WSL. Replace the candidate IPs with your own local network addresses if needed.

```bash
for ip in 127.0.0.1 <WINDOWS_HOST_IP> <LAN_IP>; do
	echo "==== $ip ===="
	curl -I --connect-timeout 5 http://$ip:7897
done
```

Add to `~/.bashrc` when the proxy address is confirmed:

```bash
export PATH="/mnt/d/software/bin:$PATH"
export http_proxy=http://<PROXY_IP>:7897
export https_proxy=http://<PROXY_IP>:7897
```

Remove apt proxy settings if needed:

```bash
sudo rm -f /etc/apt/apt.conf.d/95proxy
unset http_proxy https_proxy all_proxy
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

```bash
hf auth login
hf download google-bert/bert-base-chinese --local-dir ./bert-base-chinese
```

Example model clone:

```bash
git clone https://huggingface.co/Qwen/Qwen3-8B
```

If the download fails because of network issues, use a resumable downloader or a small Python script built around `huggingface_hub`.
