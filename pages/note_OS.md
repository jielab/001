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

## 3. WSL proxy through Clash for Windows

In Clash for Windows:

```text
Port: 7897
Allow LAN: on
System Proxy: on
```

Find the Windows proxy IP from WSL:

```bash
proxy_ip=$(awk '/nameserver/{print $2; exit}' /etc/resolv.conf)
port=7897
echo $proxy_ip
```

Test whether WSL can use the Clash proxy:

```bash
curl -I --connect-timeout 8 -x http://${proxy_ip}:${port} https://mirrors.tuna.tsinghua.edu.cn
curl -I --connect-timeout 8 -x http://${proxy_ip}:${port} https://r2u.stat.illinois.edu
```

If the test works, set proxy for normal shell commands:

```bash
cat >> ~/.bashrc <<'EOF'
export http_proxy=http://<PROXY_IP>:7897
export https_proxy=http://<PROXY_IP>:7897
export all_proxy=http://<PROXY_IP>:7897
EOF
source ~/.bashrc
```

Replace `<PROXY_IP>` with the value printed by `echo $proxy_ip`, for example `192.168.0.1`.

`sudo apt update` may not use the proxy from `.bashrc`. Set the apt proxy directly and force IPv4:

```bash
sudo tee /etc/apt/apt.conf.d/95proxy >/dev/null <<EOF
Acquire::http::Proxy "http://${proxy_ip}:${port}/";
Acquire::https::Proxy "http://${proxy_ip}:${port}/";
Acquire::ForceIPv4 "true";
EOF

apt-config dump | grep -Ei 'proxy|forceipv4'
sudo apt update
```

If `sudo apt update` shows missing GPG keys for CRAN or r2u, for example:

```text
NO_PUBKEY 51716619E084DAB9
NO_PUBKEY A1489FE2AB99A21A
```

add the keys through the same proxy:

```bash
proxy=http://${proxy_ip}:${port}

curl -fsSL -x "$proxy" "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x51716619E084DAB9" -o /tmp/cran.asc
curl -fsSL -x "$proxy" "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0xA1489FE2AB99A21A" -o /tmp/r2u1.asc
curl -fsSL -x "$proxy" "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x67C2D66C4B1D4339" -o /tmp/r2u2.asc

cat /tmp/cran.asc /tmp/r2u1.asc /tmp/r2u2.asc | \
  gpg --dearmor | \
  sudo tee /etc/apt/trusted.gpg.d/cran-r2u.gpg >/dev/null

sudo apt update
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
