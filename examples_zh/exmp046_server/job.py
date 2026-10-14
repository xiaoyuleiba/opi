#!/usr/bin/env python3

import shutil
import sys
import warnings
from pathlib import Path

from opi.core import Calculator
from opi.external_methods.server import OpiServer, ServerStatus
from opi.input.blocks import BlockMethod
from opi.input.simple_keywords import ExternalTools, Opt
from opi.input.structures.structure import Structure

# ------------------------------------------------------------------------
# 本示例需要安装服务器脚本和客户端脚本。
# 可从 https://github.com/faccts/orca-external-tools 获取这些脚本
# 请注意，umaserver 等脚本还需要额外步骤，
# 例如登录 HuggingFace 账户。
# 更多信息见以下相应教程：
# https://www.faccts.de/docs/orca/6.1/tutorials/workflows/extopt.html
# -------------------------------------------------------------------------

# 初始化 Server 对象
server = OpiServer(
    # 将下方路径改为服务器脚本的路径，例如来自 orca-external-tools 的脚本
    serverpath="/path/to/working/server/script/from/orca-external-tools"
)

# 启动服务器。将方法改为所需方法，例如 AIMNet2 或 UMA。
status = server.start_server(cmd_arguments="method")
if status is not ServerStatus.RUNNING:
    if status is ServerStatus.ALREADY_RUNNING:
        warnings.warn("Some server was already running. Using old one.")
    elif status == ServerStatus.PORT_IN_USE:
        warnings.warn("Port for server in use. Continue on your own risk.")
    else:
        print(f"Server setup returned with error {status}")
        sys.exit(1)

# 处理工作目录
working_dir = Path("RUN")
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()

# 设置计算器
structure = Structure.from_xyz("h2o.xyz")
orca_calc = Calculator(basename="extopt_ext", working_dir=working_dir)
orca_calc.structure = structure
orca_calc.input.ncores = 2
orca_calc.input.add_simple_keywords(ExternalTools.EXTOPT, Opt.OPT)
# 将下方路径改为客户端脚本的路径，例如来自 orca-external-tools 的脚本
path_client = "/path/to/working/client/script/from/orca-external-tools"
orca_calc.input.add_blocks(BlockMethod(ProgExt=path_client))
orca_calc.write_input()

# 运行计算
print("Running ORCA calculation ...", end="")
orca_calc.run()
print("   Done")

# 终止服务器
server.kill_server()

