# Windows, WSL, Python and local Transformer setup


## 1. 常识和常用命令

```powershell
# Restart Outlook
taskkill /f /im outlook.exe
start outlook

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

### 3.1 本机当前配置与职责（2026-10-01 核对）

本机代理客户端 ikuuu 使用 HTTP 代理 `http://127.0.0.1:7897`。保持客户端正常连接；调整 WSL 下载工具不需要修改 VPN、Windows 系统代理或重启桌面应用。

| 层级 | 配置位置 | 本机状态 / 作用 |
| --- | --- | --- |
| Windows 代理客户端 | ikuuu → 我的 → 常规 → 端口 | `7897`，提供本地 HTTP 代理 |
| WSL 网络 | `C:\Users\jiehu\.wslconfig` | `networkingMode=mirrored`、`dnsTunneling=true`、`autoProxy=true`，另有 `memory=72GB` |
| Windows Codex 桌面 | `C:\Users\jiehu\.codex\.env` | 保留现有配置；与 WSL shell 配置分开管理 |
| WSL 交互式命令 | `~/.bashrc` 加载 `~/.wsl_proxy.sh` | 统一使用 `7897`；提供 `proxy_on`、`proxy_off`、`proxy_status` |
| WSL aria2 | `~/.aria2/aria2.conf` | 默认使用 `7897`，关闭自动改名；无需 shell alias |
| WSL apt | `/etc/apt/apt.conf.d/95proxy` | 已配置 `7897`，保留现有重试、超时及 IPv4 设置 |

这里记录的是当前实测配置。旧记录中的 `autoProxy=false` 不是本机现状，不要照旧笔记覆盖 `.wslconfig`。本次只修改 WSL 用户配置，没有执行 `wsl --shutdown`。

### 3.2 WSL shell：只保留一套代理函数

`~/.bashrc` 只保留以下加载入口，删除其中重复定义的代理函数。代理端口统一在 `~/.wsl_proxy.sh` 管理。

```bash
# Auto proxy for WSL + Clash
[ -f ~/.wsl_proxy.sh ] && . ~/.wsl_proxy.sh
```

`~/.wsl_proxy.sh` 的当前实现：

```bash
export WSL_PROXY_HOST=127.0.0.1
export WSL_PROXY_PORT=7897
export WSL_PROXY_URL=http://${WSL_PROXY_HOST}:${WSL_PROXY_PORT}
export no_proxy=localhost,127.0.0.1,::1
export NO_PROXY=$no_proxy

_proxy_set() {
  export http_proxy=$WSL_PROXY_URL https_proxy=$WSL_PROXY_URL all_proxy=$WSL_PROXY_URL
  export HTTP_PROXY=$WSL_PROXY_URL HTTPS_PROXY=$WSL_PROXY_URL ALL_PROXY=$WSL_PROXY_URL
}

_proxy_unset() {
  unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY
}

proxy_on() { _proxy_set; echo "Proxy on: $WSL_PROXY_URL"; }
proxy_off() { _proxy_unset; echo "Proxy off"; }
proxy_status() { env | grep -Ei '^(http|https|all|no)_proxy=' | sort; }

wsl_proxy_auto() {
  if timeout 1 bash -c ":</dev/tcp/${WSL_PROXY_HOST}/${WSL_PROXY_PORT}" 2>/dev/null; then
    _proxy_set
  else
    _proxy_unset
  fi
}

wsl_proxy_auto
```

新交互式终端自动加载。已有终端运行下面的命令即可更新代理，避免重新执行 `.bashrc` 中其他环境初始化：

```bash
source ~/.wsl_proxy.sh
proxy_status
curl -I --connect-timeout 8 --max-time 20 https://github.com
```

自动检查只判断本地端口能否连接，不代表代理节点或远端服务一定正常。运行中的程序不会自动获得新环境变量；更改后需要重新启动对应命令。`proxy_off` 只影响当前 shell 及随后启动的子进程。

### 3.3 aria2：把通用默认值写入专用配置

aria2 不需要每次输入代理和自动改名参数。WSL 用户 `huangj` 的默认配置位于 `/home/huangj/.aria2/aria2.conf`：

```ini
all-proxy=http://127.0.0.1:7897
no-proxy=localhost,127.0.0.1,::1
auto-file-renaming=false
```

原来 shell 配置仍检测旧端口 `7890`，检测失败便清除了代理变量，aria2 因而直连 EBI，反复出现 `SSL/TLS handshake failure: Error in the pull function`。2026-10-01 实测：直连超时，使用 `7897` 代理时 aria2 可用性检查成功，Range 请求收到 HTTP 206 和 1024 字节数据。

修复后，即使清除所有 HTTP/HTTPS/ALL 代理环境变量，aria2 仍能通过自身默认配置完成可用性检查。后台任务同样可以使用此配置；使用其他 Linux 用户、`sudo`、`--no-conf` 或显式指定另一个 `--conf-path` 时，要另行确认读取的配置。

日常下载直接运行：

```bash
cd /mnt/g/metag
aria2c -i emp.url -d emp/ -x 4 -s 4 -j 2 -c \
  --max-tries=20 --retry-wait=10 \
  --log=emp-download.log --log-level=notice
```

`-c` 续传；保留 `.aria2` 文件。关闭自动改名可以避免产生 `.1` 等额外文件，但不等于校验已有文件完整性。并发数、日志和重试参数按任务设置。命令行参数可以覆盖默认配置，具体见 [aria2 官方手册](https://aria2.github.io/manual/en/html/aria2c.html)。

修改配置不影响已经运行的 aria2。需要生效时，在原下载终端按 Ctrl+C，等待退出，再重新运行下载命令，避免两个进程同时写同一目录。

### 3.4 Windows 桌面连接：保留已经成功的配置

此前本机 Codex 桌面登录排查中，为 `%USERPROFILE%\.codex\.env` 明确配置代理并完全退出、重新打开应用后，登录成功。没有做逐项回退实验，因此不能认定某一项是唯一原因，也不能将 Codex 的配置文件当成所有 ChatGPT 应用通用的配置入口。

当时使用的代理项如下（文件如有其他内容，应保留）：

```dotenv
HTTP_PROXY=http://127.0.0.1:7897
HTTPS_PROXY=http://127.0.0.1:7897
NO_PROXY=localhost,127.0.0.1,::1
```

`HTTPS_PROXY` 使用 `http://` 是因为本地服务是 HTTP 代理，HTTPS 流量通过 CONNECT 隧道传输。Windows 桌面应用与 WSL 的 `.bashrc` 分开管理；修复 aria2 不需要编辑这个文件、Windows 环境变量或证书设置。

此前浏览器正常并不能证明桌面登录请求也通过同一条代理路径。排查时收到 `200 Connection established` 后还需要检查实际端点响应；`405 Method Not Allowed` 只说明该端点拒绝相应请求方法，不能单独证明登录成功。

### 3.5 排查和切换网络环境

先确认出错的是哪一层，再调整对应配置：

```bash
# WSL shell 当前代理
proxy_status

# 显式测试本地代理到远端的 HTTPS 路径
curl -I --proxy http://127.0.0.1:7897 \
  --connect-timeout 8 --max-time 20 https://github.com

# aria2 默认配置与下载错误
cat ~/.aria2/aria2.conf
tail -n 40 /mnt/g/metag/emp-download.log
```

本地端口连接失败时查代理客户端；端口可连接但远端请求失败时查节点、路由或远端返回。aria2 有进程但没有新增文件时，先看日志和传输速度；它可能在续传，也可能在报错。

出国或不再使用代理时，按需要分别停用：

| 使用位置 | 停用方法 |
| --- | --- |
| WSL shell | 当前终端执行 `proxy_off`；若希望新终端也不自动启用，注释 `.bashrc` 的 helper 加载行 |
| aria2 | 注释 `~/.aria2/aria2.conf` 的 `all-proxy` 行；同时确认 shell 代理变量已清除 |
| apt | 将 `95proxy` 移到 `/etc/apt/apt.conf.d/` 之外保存，例如 `/etc/apt/95proxy.disabled` |
| Windows Codex | 仅在确实要切换桌面连接时，注释 `.codex\.env` 中代理项，保留其他内容，然后重启应用 |
| Windows 系统代理 / VPN | 仅在需要切换这一层时调整客户端或系统设置 |

aria2 的固定代理独立于 shell，`proxy_off` 不会取消它。单次临时直连可在 `proxy_off` 后使用 `aria2c --all-proxy="" ...`。关闭客户端不会自动删除 aria2 或 apt 的固定代理。

恢复时先确认客户端 `7897` 可用，再恢复需要的配置。换端口时同步检查上述表格，尤其是 `~/.wsl_proxy.sh`、`~/.aria2/aria2.conf` 和 apt 配置；不要残留 `7890`。只有实际修改 WSL 网络配置时才考虑 `wsl --shutdown`，先保存 WSL 内正在运行的工作。

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
