#!/usr/bin/env python
# coding: utf-8

# # OPI 使用入门：第一个示例
# 
# 本教程演示如何将 ORCA 与 ORCA Python 接口（OPI）配合使用。以水分子能量计算为例，编写脚本执行 OPI 计算并访问结果；这两项操作均可通过 OPI 完成。
# 
# 本 Notebook 将：
# 1. 导入所需 Python 依赖。
# 2. 定义工作目录。
# 3. 设置输入结构。
# 4. 定义计算设置。
# 5. 执行 ORCA 计算。
# 6. 检查正常终止与收敛。
# 7. 处理结果。

# ## 步骤 1：导入依赖
# 
# 使用 OPI 时，首先导入所需模块。此外，还导入 `py3Dmol` 等可视化模块。可能需要将 `py3Dmol` 安装到 OPI 的 `venv` 中，例如激活 `.venv` 后运行 `uv pip install py3Dmol`。

# In[1]:


from pathlib import Path
import shutil

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import Dft, Task
from opi.input.structures.structure import Structure

# > 用于分子可视化
import py3Dmol


# ## 步骤 2：工作目录
# 
# 使用 OPI 时，建议定义一个子目录来存放实际 ORCA 计算。本例创建 `RUN` 目录。

# In[2]:


# > 在 `RUN` 中执行计算
working_dir = Path("RUN")
# > 自动创建或重新创建 `working_dir`
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()


# ## 步骤 3：设置输入结构
# 
# 以埃为单位定义结构的笛卡尔坐标，并将其可视化：
# 

# In[3]:


# > 将以 angstrom 为单位的笛卡尔坐标定义为 Python 字符串
xyz_data = """\
3

O      0.00000   -0.00000    0.00000
H      0.00000    0.96899    0.00000
H      0.93966   -0.23409    0.03434\n
"""
# > 可视化输入结构
view = py3Dmol.view(width=400, height=400)
view.addModel(xyz_data, 'xyz')
view.setStyle({}, {'stick': {}, 'sphere': {'scale': 0.3}})
view.zoomTo()
view.show()

# > 将输入结构写入文件
with open(working_dir / "struc.xyz","w") as f:
    f.write(xyz_data)
# > 将结构读入对象
structure = Structure.from_xyz(working_dir / "struc.xyz")


# ## 步骤 4：定义计算设置
# 接下来定义 ORCA 优化使用的设置。可以创建一个 calculator 对象，并向其添加简单关键词：

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 1) -> Calculator:
    # > 设置 Calculator 对象，同时设置 ORCA 计算的基本名称
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > 将结构赋给计算器
    calc.structure = structure

    # 定义理论水平：使用 r2SCAN-3c 进行几何优化
    sk_list = [
        Dft.R2SCAN_3C, # > r2SCAN-3c 方法（附带预定义基组）
        Task.SP # > 执行单点计算
    ]

    # > 在计算器中使用简单关键词
    calc.input.add_simple_keywords(*sk_list)

    # > 定义计算使用的 CPU 数量
    calc.input.ncores = ncores # > 本次 ORCA 计算使用的 CPU 数量（默认值：1）

    return calc

calc = setup_calc("single_point", working_dir=working_dir, structure=structure)


# ## 步骤 5：执行计算
# 
# 执行单点计算时，先将输入写成 ORCA 格式，再使用 calculator 类的 `run()` 函数运行作业。可通过其 `get_output()` 函数将 ORCA 输出写入 calculator 对象。

# In[5]:


def run_calc(calc: Calculator) -> Output:
    # > 写出 ORCA 输入文件
    calc.write_input()
    # > 运行 ORCA 计算
    print("Running ORCA calculation ...", end="")
    calc.run()
    print("   Done")

    # > 获取输出对象
    output = calc.get_output()

    return output

output = run_calc(calc)


# ## 步骤 6：检查正常终止与收敛
# 执行自动化计算时，应始终检查是否正常终止和收敛。简单检查示例如下：

# In[6]:


def check_and_parse_output(output: Output):
    # > 检查 ORCA 是否正常终止
    status = output.terminated_normally()
    if not status:
        # > ORCA 未正常终止
        raise RuntimeError("ORCA did not terminate normally. Please check RUN/single_point.out")
    else:
        # > ORCA 已正常终止，可以解析输出
        output.parse()

    # 现在检查 SCF 是否收敛
    if output.results_properties.geometries[0].single_point_data.converged:
        print("SCF CONVERGED")
    else:
        raise RuntimeError("SCF DID NOT CONVERGE")

check_and_parse_output(output)


# ## 步骤 7：访问结果
# 计算成功后，可以访问所得性质。通常可通过 `output` 对象获取：

# In[7]:


# > 以 Hartree 为单位的总能量
Etot = output.results_properties.geometries[0].single_point_data.finalenergy

# > 输出总能量
print("Final energy in Hartee: ",Etot)


# ## 小结
# 
# 本 Notebook 演示了如何使用 ORCA Python 接口（OPI）执行 DFT 单点能计算。

