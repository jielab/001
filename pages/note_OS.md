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

### 3.1 本次成功记录（2026-10-01）

出国期间卸载 Clash 并移除代理设置；回来后改用 ikuuu。网页版 ChatGPT 正常，但 Windows 桌面应用登录报错：

```text
token_exchange_failed
Token exchange failed: error sending request for url (https://auth.openai.com/oauth/token)
```

本次排查与结果：

1. 将 ikuuu 本地端口从 `7890` 改为原来的 `7897`，断开后重新连接；仅改端口后仍然报错。
2. 通过显式代理 `http://127.0.0.1:7897` 测试认证端点，收到 `200 Connection established`，随后收到 `405 Method Not Allowed`，确认这条代理路径能完成 HTTPS 请求。
3. 检查 Windows 的 Process/User/Machine 环境变量：未发现 HTTP_PROXY、HTTPS_PROXY、ALL_PROXY、NO_PROXY、SSL_CERT_FILE、CODEX_CA_CERTIFICATE、CODEX_HOME；当时也没有 `%USERPROFILE%\.codex\.env`。
4. 创建 Windows 的 `%USERPROFILE%\.codex\.env`，明确指定代理和 localhost 例外；完全退出并重新打开桌面应用后，登录成功。

**最有力的修复线索是第 4 步：为 Codex 明确配置代理并重启应用。没有做逐项回退实验，不能断言单个设置是唯一原因。此次没有通过修改 WSL 的 `.bashrc` 修复桌面登录。**

### 3.2 代理客户端：统一使用 7897

ikuuu：**我的 → 常规 → 端口 → `7897`**，保存后断开并重新连接。

本次排查使用 **TUN + 全局**。这是成功过程中的配置记录，不代表两者永远必需：

- **全局 / 规则**：决定流量如何分流。规则正确时也可以使用规则模式；排查时先保持已验证的模式。
- **TUN / 系统代理**：决定流量如何被接管，与全局/规则是不同的设置。
- 换回 Clash 或其他客户端时，确认其 HTTP 或混合代理端口仍是 `7897`；不要仅凭客户端名称假定端口相同。

在 **Windows PowerShell** 验证：

```powershell
Test-NetConnection 127.0.0.1 -Port 7897
curl.exe --proxy http://127.0.0.1:7897 --noproxy localhost,127.0.0.1 -I --connect-timeout 10 --max-time 20 https://auth.openai.com/oauth/token
```

`TcpTestSucceeded : True` 只说明本地端口可连接。curl 在代理的 `200 Connection established` 之后返回 `405 Method Not Allowed`，说明认证端点已响应 HEAD 请求；不代表实际 OAuth POST 登录必然成功。

### 3.3 Windows ChatGPT / Codex 桌面版：明确指定代理

**文件位于 Windows 用户目录，不是 WSL 的 Linux 用户目录：**

```text
C:\Users\jiehu\.codex\.env
```

其他 Windows 用户对应 `%USERPROFILE%\.codex\.env`。若自行设置过 `CODEX_HOME`，应检查实际配置目录。

在 **Windows PowerShell** 打开配置文件：

```powershell
$codexDir = Join-Path $env:USERPROFILE '.codex'
New-Item -ItemType Directory -Force -Path $codexDir | Out-Null
notepad.exe (Join-Path $codexDir '.env')
```

写入以下三行并保存；文件已经存在时，只更新同名配置，保留其他内容，避免重复条目。确认文件名是 `.env`，不是 `.env.txt`。

```dotenv
HTTP_PROXY=http://127.0.0.1:7897
HTTPS_PROXY=http://127.0.0.1:7897
NO_PROXY=localhost,127.0.0.1,::1
```

`HTTPS_PROXY` 的值仍以 `http://` 开头：这里指定的是本地 HTTP 代理，HTTPS 请求通过它建立隧道。

保存后：

1. 完全退出 ChatGPT / Codex 桌面应用；任务管理器中确认没有残留的 `ChatGPT.exe`、`Codex.exe` 或 `codex.exe`。先保存正在进行的工作。
2. 保持代理客户端连接，端口为 `7897`。
3. 重新打开桌面应用，重新点击“继续登录”，不要刷新旧的错误页面。

**浏览器能访问 ChatGPT，不等于桌面应用的认证请求走了同一条网络路径。Windows 桌面应用通常不会读取 WSL 的 `.bashrc`。**

参考：[本次采用的 GitHub 排查线索：openai/codex #26764](https://github.com/openai/codex/issues/26764)。其中包含不同原因的用户报告；本节以本机实际成功过程为准。

### 3.4 WSL2 网络与命令行代理（独立配置）

以下保留原有 WSL 方案，供 Git、pip、R、Hugging Face 和 apt 使用；本次桌面登录修复没有重新验证这些 WSL 步骤。

在 Windows 的 `%USERPROFILE%\.wslconfig` 中合并以下设置；若已有 `[wsl2]`，修改对应项目，不要重复创建同名段落或覆盖其他设置：

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
autoProxy=false
```

该方案依赖支持 mirrored networking 的 Windows/WSL 版本。保存 WSL 内正在运行的工作后，在 **Windows PowerShell** 执行：

```powershell
wsl --shutdown
```

重新打开 WSL。在 **WSL 的 `~/.bashrc`** 中更新或加入以下内容，不要保留其他指向旧端口的重复代理设置：

```bash
export http_proxy=http://127.0.0.1:7897
export https_proxy=$http_proxy
export all_proxy=$http_proxy
export no_proxy=localhost,127.0.0.1,::1
```

在 **WSL** 生效并测试：

```bash
source ~/.bashrc
curl -I --connect-timeout 8 --max-time 20 https://github.com
```

apt 使用独立配置，在 **WSL** 执行：

```bash
sudo tee /etc/apt/apt.conf.d/95proxy >/dev/null <<'EOF'
Acquire::http::Proxy "http://127.0.0.1:7897";
Acquire::https::Proxy "http://127.0.0.1:7897";
EOF
sudo apt update
```

### 3.5 出国不需要代理时，以及下次换客户端时

**固定代理仍指向本机端口时，关闭或卸载代理客户端会使依赖它的程序无法联网。切换网络环境时，要同步停用或恢复对应配置。**

| 使用位置 | 需要检查的配置 |
| --- | --- |
| Windows ChatGPT / Codex | `%USERPROFILE%\.codex\.env` 中的 HTTP_PROXY、HTTPS_PROXY；若另有 ALL_PROXY 也要检查 |
| WSL 常用命令 | `~/.bashrc` 中的代理 export，以及其他启动文件中自定义的代理设置 |
| WSL apt | `/etc/apt/apt.conf.d/95proxy` |
| Windows 系统代理 | 设置 → 网络和 Internet → 代理；检查是否还指向已关闭的本地端口 |

不需要代理时：

- Windows Codex：将 `.env` 中的 `HTTP_PROXY`、`HTTPS_PROXY`（以及另行设置的 `ALL_PROXY`）行前加 `#` 注释，保留其他内容；完全退出并重新打开应用。
- WSL：注释 `.bashrc` 中的代理 export，再在当前终端运行下面的 unset；仅修改文件不会清除当前终端已经存在的变量。

```bash
unset http_proxy https_proxy all_proxy no_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY NO_PROXY
```

- apt：将 `95proxy` 备份到 `/etc/apt/apt.conf.d/` 之外，例如：

```bash
sudo mv /etc/apt/apt.conf.d/95proxy /etc/apt/95proxy.disabled
```

重新需要代理时：先启动代理客户端并确认 `7897` 可用，再恢复 `.env` 和 `.bashrc` 中的代理行；运行 `source ~/.bashrc`，将 apt 配置移回原处，并重启桌面应用。

若仍报错，先重复 3.2 的显式代理测试：本地端口连接失败时查客户端；收到端点响应但桌面应用仍失败时查 `.codex\.env`、环境变量及应用日志。不要把所有 `token_exchange_failed` 都当成同一个原因；若详情明确出现 `403 Country, region, or territory not supported`，那与本次“请求发送失败”是不同的返回结果。

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
