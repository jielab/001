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

## 3. New Windows/WSL network setup: Clash proxy, apt, R, and SSH

This section is for a new Windows computer. Follow the order below. The goal is to make all of these work from Ubuntu/WSL:

```text
sudo apt update
R package installation from CRAN/GitHub
ssh <USER>@<SERVER_IP>
```

Do not put real internal IPs or usernames in a public GitHub page. Use placeholders such as `<SERVER_IP>`, `<USER>`, and `<PROXY_PORT>`.

### 3.1 Turn on Clash settings in Windows

In Clash for Windows or Clash Verge:

```text
System Proxy: on
Allow LAN: on
Port: 7897
```

Reason:

- `System Proxy: on` lets Windows programs use Clash.
- `Allow LAN: on` lets WSL reach the Clash proxy when WSL is using a separate network interface.
- The port must match the port shown in Clash. In the examples below, the port is `7897`.

Check the port in Windows PowerShell:

```powershell
netstat -ano | findstr :7897
```

### 3.2 Install WSL2 and enable mirrored networking

In Windows PowerShell:

```powershell
wsl --list --online
wsl --install -d Ubuntu-24.04
wsl --set-default-version 2
wsl --update
wsl --version
wsl --status
```

Create or edit the WSL configuration file:

```powershell
notepad C:\Users\<WINDOWS_USER>\.wslconfig
```

Add:

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
autoProxy=true
```

Then restart WSL:

```powershell
wsl --shutdown
```

Reason:

- WSL2 normally uses NAT networking, which can fail to inherit Windows VPN or internal-network routes.
- `networkingMode=mirrored` makes WSL follow the Windows network more closely. This is important for internal servers such as `ssh <USER>@<SERVER_IP>`.
- `dnsTunneling=true` helps WSL DNS work correctly under VPN.
- `autoProxy=true` helps WSL detect Windows proxy settings.

After reopening Ubuntu, test the internal server:

```bash
unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY

ip route get <SERVER_IP>
ping -c 4 <SERVER_IP>
ssh <USER>@<SERVER_IP>
```

If `ping` works but this fails:

```bash
curl -I http://<SERVER_IP>
```

that usually only means port 80 is not open on the server. It does not mean SSH is broken. Test SSH directly:

```bash
nc -vz <SERVER_IP> 22
ssh <USER>@<SERVER_IP>
```

### 3.3 Find the proxy address that works in WSL

After mirrored networking, `127.0.0.1:7897` usually works. Test it first:

```bash
port=7897

curl -I --connect-timeout 8 -x http://127.0.0.1:${port} https://mirrors.tuna.tsinghua.edu.cn
curl -I --connect-timeout 8 -x http://127.0.0.1:${port} https://cloud.r-project.org
curl -I --connect-timeout 8 -x http://127.0.0.1:${port} https://api.github.com
```

If `127.0.0.1` does not work, test the WSL nameserver address:

```bash
proxy_ip=$(awk '/nameserver/{print $2; exit}' /etc/resolv.conf)
port=7897

echo $proxy_ip
curl -I --connect-timeout 8 -x http://${proxy_ip}:${port} https://mirrors.tuna.tsinghua.edu.cn
curl -I --connect-timeout 8 -x http://${proxy_ip}:${port} https://cloud.r-project.org
curl -I --connect-timeout 8 -x http://${proxy_ip}:${port} https://api.github.com
```

Use the proxy address that works. In mirrored mode, this is often:

```text
http://127.0.0.1:7897
```

### 3.4 Set proxy for normal Linux shell commands

Add this to the end of `~/.bashrc`. Use the proxy address that worked in the previous step.

```bash
cat >> ~/.bashrc <<'EOF'
export http_proxy=http://127.0.0.1:7897
export https_proxy=http://127.0.0.1:7897
export all_proxy=http://127.0.0.1:7897

# Do not send local or internal servers through Clash
export no_proxy=localhost,127.0.0.1,::1,<SERVER_IP>
export NO_PROXY=$no_proxy
EOF

source ~/.bashrc
```

Reason:

- `http_proxy`, `https_proxy`, and `all_proxy` are used by `curl`, `git`, R, Python, and many command-line tools.
- `no_proxy` prevents internal servers from being sent to Clash. This is important for `ssh <USER>@<SERVER_IP>` and other intranet services.

Test:

```bash
env | grep -i proxy
curl -I https://cloud.r-project.org
curl -I https://api.github.com
ssh <USER>@<SERVER_IP>
```

### 3.5 Set apt proxy and force IPv4

Do not rely on `.bashrc` for `sudo apt update`. Set the apt proxy directly:

```bash
sudo tee /etc/apt/apt.conf.d/95proxy >/dev/null <<'EOF'
Acquire::http::Proxy "http://127.0.0.1:7897/";
Acquire::https::Proxy "http://127.0.0.1:7897/";
Acquire::ForceIPv4 "true";
EOF

apt-config dump | grep -Ei 'proxy|forceipv4'
sudo apt update
```

Reason:

- `sudo apt update` may not inherit the proxy variables from `.bashrc`.
- Some mirrors may try IPv6 first and get stuck, for example at an IPv6 address such as `2402:...`.
- `Acquire::ForceIPv4 "true";` avoids that problem.

If `127.0.0.1:7897` does not work on a specific computer, replace it with the working proxy address found above, for example:

```text
http://<PROXY_IP>:7897/
```

### 3.6 Fix missing GPG keys for CRAN and r2u

If `sudo apt update` works but shows warnings like these:

```text
NO_PUBKEY 51716619E084DAB9
NO_PUBKEY A1489FE2AB99A21A
```

add the keys through the proxy:

```bash
proxy=http://127.0.0.1:7897

curl -fsSL -x "$proxy" "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x51716619E084DAB9" -o /tmp/cran.asc
curl -fsSL -x "$proxy" "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0xA1489FE2AB99A21A" -o /tmp/r2u1.asc
curl -fsSL -x "$proxy" "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x67C2D66C4B1D4339" -o /tmp/r2u2.asc

cat /tmp/cran.asc /tmp/r2u1.asc /tmp/r2u2.asc | \
  gpg --dearmor | \
  sudo tee /etc/apt/trusted.gpg.d/cran-r2u.gpg >/dev/null

sudo apt update
```

Reason:

- The network can already be fixed, but apt still refuses CRAN/r2u repositories if their signing keys are missing.
- Using `curl -x "$proxy"` is more stable than `gpg --recv-keys` behind Clash.

### 3.7 R package installation and R update

Start R from the same Ubuntu shell after `source ~/.bashrc`. In R, check that the proxy is not empty:

```r
Sys.getenv(c("http_proxy", "https_proxy", "all_proxy",
             "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY"))
```

The proxy should look like this:

```text
http://127.0.0.1:7897
```

It should not look like this:

```text
http://:7897
```

If R shows the wrong proxy, fix it temporarily inside R:

```r
Sys.setenv(
  http_proxy  = "http://127.0.0.1:7897",
  https_proxy = "http://127.0.0.1:7897",
  all_proxy   = "http://127.0.0.1:7897"
)
Sys.unsetenv(c("HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY"))
```

Test CRAN access from R:

```r
download.file("https://cloud.r-project.org/src/contrib/PACKAGES", tempfile())
```

Install packages from CRAN:

```r
install.packages("lintr", repos = "https://cloud.r-project.org")
install.packages("remotes", repos = "https://cloud.r-project.org")
```

Install packages from GitHub:

```r
remotes::install_github("r-lib/lintr")
```

Update installed R packages:

```r
update.packages(ask = FALSE, checkBuilt = TRUE, repos = "https://cloud.r-project.org")
```

If R often starts without the proxy, create `~/.Renviron`:

```bash
cat > ~/.Renviron <<'EOF'
http_proxy=http://127.0.0.1:7897
https_proxy=http://127.0.0.1:7897
all_proxy=http://127.0.0.1:7897
no_proxy=localhost,127.0.0.1,::1,<SERVER_IP>
NO_PROXY=localhost,127.0.0.1,::1,<SERVER_IP>
EOF
```

Reason:

- R uses the same proxy environment variables as Linux shell tools.
- `install.packages()` fails with `cannot open URL` when R cannot reach CRAN.
- `remotes::install_github()` fails with `Couldn't resolve proxy name` when proxy variables are malformed, for example `http://:7897`.

### 3.8 Quick diagnosis commands

Check proxy variables:

```bash
env | grep -i proxy
```

Check apt proxy:

```bash
apt-config dump | grep -Ei 'proxy|forceipv4'
```

Check internal-server route:

```bash
ip route get <SERVER_IP>
ping -c 4 <SERVER_IP>
nc -vz <SERVER_IP> 22
```

Check whether a server has HTTP on port 80:

```bash
curl -I --connect-timeout 5 http://<SERVER_IP>
```

If `ping` and SSH work but `curl http://<SERVER_IP>` fails, the server probably does not run an HTTP service on port 80.

Remove apt proxy if needed:

```bash
sudo rm -f /etc/apt/apt.conf.d/95proxy
sudo apt update
```

Temporarily remove shell proxy:

```bash
unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY
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
