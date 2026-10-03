# 安装 OPI

OPI 提供一个 Python 库，用于访问 ORCA 并处理其计算结果。目前通过 [PyPI](https://pypi.org/project/orca-pi/) 和 [GitHub](https://github.com/faccts/opi) 发布。

可以直接从 PyPI 安装 OPI：

```
pip install orca-pi
```

也可以从 GitHub 安装：

```
git clone https://github.com/faccts/opi.git
cd opi
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install .
```

建议在 _venv_ 虚拟环境中安装 OPI，以便正确使用。

## ORCA 与 Open MPI

OPI 的大部分功能要求主机上安装 ORCA。
ORCA 必须从我们的[下载门户](https://www.faccts.de/orca/)（商业用户和学术用户）或 [ORCA 论坛](https://orcaforum.kofo.mpg.de)（学术用户）单独获取。
并行计算还需要 [Open MPI](https://www.open-mpi.org/)。
在大多数现代操作系统上，通常可直接通过系统包管理器安装 Open MPI。
否则，需要从其网站获取合适版本并编译。

**注意：OPI 首次随 ORCA 6.1 引入，不兼容更早版本。
最低支持的 ORCA 版本始终记录在 [ORCA_MINIMAL_VERSION](https://github.com/faccts/opi/blob/main/src/opi/__init__.py) 中。**

## 配置文件

可通过系统路径（`$PATH` 和 `$LD_LIBRARY_PATH`）或配置文件向 OPI 提供 ORCA 和 Open MPI 的位置。
配置文件应位于 `$XDG_CONFIG_HOME/opi/config.toml`（通常 `$XDG_CONFIG_HOME` 指向 `$HOME/.config`）。

`config.toml` 文件应包含 ORCA 和 Open MPI 基目录的绝对路径。
这里的基目录指包含相应 `bin/` 和 `lib/` 目录的父目录。

```
ORCA_PATH = "/path/to/orca/"
# > Path to Open MPI installed with `apt`
MPI_PATH = "/usr"
```
