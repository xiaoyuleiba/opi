#!/usr/bin/env python
# coding: utf-8

# # 原子 Fukui 指数
# 
# 本教程演示如何使用 ORCA Python 接口（OPI）计算原子 Fukui 指数。
# 
# 本 Notebook 将：
# 1. 导入所需依赖。
# 2. 定义工作目录并准备结构。
# 3. 可视化初始分子结构。
# 4. 运行 ORCA 计算。
# 5. 提取电荷并计算 Fukui 函数。
# 6. 对 Fukui 函数进行三维可视化。

# ## 步骤 1：导入依赖
# 
# 首先导入以下用途所需的模块：
# - 访问 ORCA 输入与输出
# - 数值计算与数据处理
# - 绘制结果
# - 目录管理
# 
# > **注意：** 我们还导入 `py3Dmol` 等可视化与绘图模块。可能需要将 `py3Dmol` 安装到 OPI 的 `venv` 中，例如激活 `.venv` 后运行 `uv pip install py3Dmol`。

# In[ ]:


# > 导入 pathlib，用于目录处理
from pathlib import Path
import shutil

# > 导入 numpy，用于数据处理和数值运算
import numpy as np

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, Dft, DispersionCorrection
from opi.input.structures.structure import Structure

# > 导入可视化库
import py3Dmol
from IPython.display import display, HTML


# ## 步骤 2：定义工作目录
# 
# 所有实际计算均在子目录 `RUN` 中执行。

# In[ ]:


# > 在 `RUN` 中执行计算
working_dir = Path("RUN")
# > 自动创建或重新创建 `working_dir`
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()


# ## 步骤 3：设置输入结构
# 
# 以甲醛为示例分子。用 XYZ 格式定义笛卡尔坐标下的三维结构，并将其可视化。

# In[3]:


# > 将分子以 Angstrom 为单位的笛卡尔坐标定义为 Python 字符串
xyz_data = """\
4

C      0.59761    0.92710    0.00015
O      1.78965    0.92659    0.00044
H      0.00000    1.85448    0.00000
H      0.00000    0.00000    0.00000\n
"""

# > 使用 py3Dmol 可视化分子结构
view = py3Dmol.view(width=500, height=500)
view.addModel(xyz_data, "xyz")
view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})
view.setBackgroundColor('white')
view.zoomTo()
view.show()

# > 将 XYZ 结构保存到文件
with open(working_dir / "struc.xyz","w") as f:
    f.write(xyz_data)


# ## 步骤 4：ORCA 计算
# 
# 对分子的三种电子态进行单点能计算：
# - 中性态（默认）
# - 阴离子态（增加一个电子，charge = -1，multiplicity = 2）
# - 阳离子态（移去一个电子，charge = +1，multiplicity = 2）

# 需要检查计算是否成功并解析输出，这由 `check_and_parse_output` 完成：

# In[4]:


def check_and_parse_output(output: Output)-> None:
    # > 检查 ORCA 是否正常终止
    status = output.terminated_normally()
    if not status:
        # > ORCA 未正常终止
        raise RuntimeError(f"ORCA did not terminate normally, see output file: {output.get_outfile()}")
    else:
        # > ORCA 已正常终止，可以解析输出
        output.parse()

    # 现在检查 SCF 是否收敛
    if output.results_properties.geometries[0].single_point_data.converged:
        print("ORCA terminated normally and the SCF converged")
    else:
        raise RuntimeError("SCF DID NOT CONVERGE")


# 先进行中性态计算：

# In[5]:


# > 结构文件的路径
xyz_file = working_dir / "struc.xyz"

# > 创建用于生成 ORCA 输入并执行计算的 Calculator 对象
calc = Calculator(basename="neutral", working_dir=working_dir)

# > 从 XYZ 文件加载分子结构
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 0
calc.structure.multiplicity = 1

# > 添加计算关键词
calc.input.add_simple_keywords(
    Dft.B3LYP,
    BasisSet.DEF2_SVP,
    DispersionCorrection.D3
)

# > 写出 ORCA 输入文件并运行计算
calc.write_input()
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")

# > 获取输出对象
output_neutral = calc.get_output()
check_and_parse_output(output_neutral)


# 接下来进行阴离子计算：

# In[6]:


calc = Calculator(basename="anion", working_dir=working_dir)
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = -1
calc.structure.multiplicity = 2
calc.input.add_simple_keywords(
    Dft.B3LYP,
    BasisSet.DEF2_SVP,
    DispersionCorrection.D3
)
calc.write_input()
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")
output_anion = calc.get_output()
check_and_parse_output(output_anion)


# 最后进行阳离子计算：

# In[7]:


calc = Calculator(basename="cation", working_dir=working_dir)
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 1
calc.structure.multiplicity = 2
calc.input.add_simple_keywords(
    Dft.B3LYP,
    BasisSet.DEF2_SVP,
    DispersionCorrection.D3
)
calc.write_input()
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")
output_cation = calc.get_output()
check_and_parse_output(output_cation)


# ## 步骤 5：提取电荷并计算 Fukui 函数

# Fukui 函数有三种常用形式：
# 
# - f⁺（亲核进攻）：ρ(N+1) - ρ(N)
# - f⁻（亲电进攻）：ρ(N) - ρ(N-1)
# - f⁰（自由基进攻）：½[ρ(N+1) - ρ(N-1)]
# 
# 这里采用 Mulliken 电荷计算原子 Fukui 指数，详见 [J. Am. Chem. Soc. 1986, 108, 19, 5708–5711](https://doi.org/10.1021/ja00279a008)。注意：使用 Mulliken 电荷代替布居或密度时，上述所有公式都需要引入负号。

# In[8]:


def compute_fukui_charge(charges1: list[float], charges2: list[float], mode: str) -> list[float]:
    """根据原子电荷计算 Fukui 函数。"""
    q1 = np.array(charges1)
    q2 = np.array(charges2)

    # > 由于使用电荷而非布居或密度，因此颠倒相减顺序
    if mode == 'plus/minus':
        diff_charge = (q2 - q1)
    elif mode == 'zero':
        diff_charge = (q2 - q1) / 2
    return diff_charge.tolist()

# > 从 Mulliken 布居分析中提取原子电荷
charges_neutral = [c[0] for c in output_neutral.results_properties.geometries[-1].mulliken_population_analysis[-1].atomiccharges]
charges_anion = [c[0] for c in output_anion.results_properties.geometries[-1].mulliken_population_analysis[-1].atomiccharges]
charges_cation = [c[0] for c in output_cation.results_properties.geometries[-1].mulliken_population_analysis[-1].atomiccharges]

# > 根据 Mulliken 电荷计算原子 Fukui 指数，注意负号
# > 因为此处使用的是原子电荷，而非密度
f_plus_charges = compute_fukui_charge(charges1=charges_anion, charges2=charges_neutral, mode='plus/minus')
f_minus_charges = compute_fukui_charge(charges1=charges_neutral, charges2=charges_cation, mode='plus/minus')
f_zero_charges = compute_fukui_charge(charges1=charges_anion, charges2=charges_cation, mode='zero')


# ## 步骤 6：Fukui 函数的三维可视化

# In[9]:


def read_xyz(xyz_path: str) -> tuple[list[tuple[str, float, float, float]], str]:
    """读取 XYZ 文件，返回原子信息和原始字符串。"""
    with open(xyz_path, "r") as f:
        lines = f.readlines()
    num_atoms = int(lines[0])
    atom_lines = lines[2:2 + num_atoms]
    coords = []
    for line in atom_lines:
        parts = line.strip().split()
        atom = parts[0]
        x, y, z = map(float, parts[1:])
        coords.append((atom, x, y, z))
    return coords, "".join(lines)

def create_molecule_view(xyz_file: str, title: str, charges: list[float] = None) -> str:
    """使用 py3Dmol 生成可交互的分子三维视图。"""
    atoms, xyz_str = read_xyz(xyz_file)
    view = py3Dmol.view(width=300, height=300)
    view.addModel(xyz_str, "xyz")
    view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})

    for (atom, x, y, z), charge in zip(atoms, charges):
        view.addLabel(f"{charge:.3f}", {
            'position': {'x': x, 'y': y, 'z': z},
            'backgroundColor': 'black',
            'fontColor': 'white',
            'showBackground': True,
            'fontSize': 12
        })
    title += " - Label"

    view.zoomTo()
    view.zoom(1.8)
    view.setBackgroundColor('white')
    html = view._make_html()

    return f"""
    <div style="display:inline-block; text-align:center; margin:10px;">
        <h4>{title}</h4>
        {html}
    </div>
        """

# 定义 XYZ 文件路径
xyz_path = working_dir / "struc.xyz"

# 生成可视化结果
label_views = [
    create_molecule_view(xyz_path, title="Mulliken charges", charges=charges_neutral),
    create_molecule_view(xyz_path, title="F+", charges=f_plus_charges),
    create_molecule_view(xyz_path, title="F-", charges=f_minus_charges),
    create_molecule_view(xyz_path, title="F0", charges=f_zero_charges),
]

# 以 2x2 网格显示带标签的视图
label_html = f"""
<h3>Atomic Fukui Indices (Labels)</h3>
<div style='display: flex; flex-wrap: wrap; justify-content: left; gap: 10px;'>
    {label_views[0]} {label_views[1]}
    {label_views[2]} {label_views[3]}
</div>
"""
display(HTML(label_html))


# Fukui 指数表明，亲核进攻倾向于发生在碳原子上，亲电进攻倾向于发生在氧原子上，自由基进攻也倾向于发生在氧原子上。

# ## 小结
# 
# 本 Notebook 使用 ORCA Python 接口（OPI），基于 Mulliken 电荷计算并可视化原子 Fukui 指数。也可以使用其他电荷方案，如 Loewdin、Mayer。

