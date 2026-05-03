# R setup and package installation notes

## 1. Library paths

On Windows, set `R_LIBS_USER` in environment variables. On Linux/WSL, set it in `~/.Renviron`.

Check current library paths:

```r
.libPaths()
```

Example WSL configuration:

```bash
mkdir -p /mnt/d/R_lib/Linux
echo 'R_LIBS_USER=/mnt/d/R_lib/Linux' >> ~/.Renviron
echo '.libPaths(c("/mnt/d/R_lib/Linux", .libPaths()))' >> ~/.Rprofile
```

## 2. Useful R resources

- [R Graph Gallery](https://r-graph-gallery.com/index.html)
- [COVID-19 R resources](https://statsandr.com/blog/top-r-resources-on-covid-19-coronavirus/)
- [Global Environmental Health code and data](https://globalenvhealth.org/code-data-download/)
- [Jokergoo bioinformatics visualization tools](https://jokergoo.github.io/software/)
- [梁志生 R 包项目](https://gitee.com/sheng0825/projects)

## 3. When `install.packages()` fails

Typical error:

```text
cannot open URL 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/PACKAGES'
package 'pacman' is not available for this version of R
```

This usually means R cannot access the CRAN index. It does not necessarily mean the package is incompatible with your R version.

Try:

```r
options(repos = c(CRAN = "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 600)
```

## 4. Install R and R packages through apt/r2u on Ubuntu 24.04

For WSL/Ubuntu, binary R packages from apt/r2u are often faster and more stable than compiling from source.

```bash
echo "deb [trusted=yes] https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/linux/ubuntu noble-cran40/" | \
	sudo tee /etc/apt/sources.list.d/cran-r.list > /dev/null

sudo apt update
sudo apt-mark unhold r-base-core r-base r-base-dev r-recommended

sudo apt install --no-install-recommends \
	r-base-core=4.6.0-2.2404.0 \
	r-base=4.6.0-2.2404.0 \
	r-base-dev=4.6.0-2.2404.0 \
	r-recommended=4.6.0-2.2404.0

sudo apt-mark hold r-base-core r-base r-base-dev r-recommended
```

Set up the R library path:

```bash
mkdir -p /mnt/d/R_lib/Linux
echo 'R_LIBS_USER=/mnt/d/R_lib/Linux' >> ~/.Renviron
echo '.libPaths(c("/mnt/d/R_lib/Linux", .libPaths()))' >> ~/.Rprofile
```

Add r2u:

```bash
sudo tee /etc/apt/sources.list.d/r2u.list > /dev/null <<'EOF_R2U'
deb [trusted=yes arch=amd64] https://r2u.stat.illinois.edu/ubuntu noble main
EOF_R2U

sudo tee /etc/apt/preferences.d/99r2u > /dev/null <<'EOF_R2U_PIN'
Package: *
Pin: release o=CRAN-Apt Project
Pin: release l=CRAN-Apt Packages
Pin-Priority: 700
EOF_R2U_PIN

sudo apt update
```

Install common R packages:

```bash
sudo apt install --no-install-recommends \
	r-cran-pacman r-cran-data.table r-cran-tidyverse r-cran-survival \
	r-cran-broom r-cran-ggplot2 r-cran-patchwork r-cran-readxl \
	r-cran-writexl r-cran-purrr r-cran-dplyr r-cran-stringr \
	r-cran-lubridate

R --version
```
