#!/usr/bin/env python
# coding: utf-8

# # ADLD：HFLD
# 
# 本教程演示如何将 ORCA 与 ORCA Python 接口（OPI）配合使用，通过**[伦敦色散能的原子分解（ADLD）](https://doi.org/10.1021/acscentsci.5c00356)**分析，将分子二聚体中 **[HFLD](https://doi.org/10.1021/acs.jctc.9b00425)** 的伦敦色散相互作用能分解为原子贡献。HFLD 即 Hartree-Fock 加伦敦色散，分解其能量可非经验地估计各原子的伦敦色散贡献。关于速度更快但更具经验性的 ADLD(D4) 分解，请参阅相应 Notebook。与 ADLD(D4) 相比，这里只需计算一次完整体系，因为仅包含一个片段的子体系没有贡献。
# 
# 本 Notebook 将：
# 1. 导入所需 Python 依赖。
# 2. 定义工作目录。
# 3. 设置输入结构。
# 4. 对超体系运行 HFLD 计算。
# 7. 处理 ADLD 结果。
# 8. 可视化结果。
# 9. 生成 Cube 文件。

# ## 步骤 1：导入依赖
# 
# 首先导入以下用途所需的模块：
# - 访问 ORCA 输入与输出
# - 数值计算与数据处理
# - 绘制结果
# 
# > **注意：** 我们还导入 `py3Dmol` 等可视化与绘图模块。可能需要将 `py3Dmol` 安装到 OPI 的 `venv` 中，例如激活 `.venv` 后运行 `uv pip install py3Dmol`。

# In[1]:


from pathlib import Path
import shutil

# > 使用 pandas 和 numpy 处理数据
import numpy as np

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, Dlpno, \
    Scf, AuxBasisSet, Approximation
from opi.input.structures.structure import Structure
from opi.input.blocks import BlockFrag

# > 用于分子可视化
import py3Dmol


# ## 步骤 2：工作目录与换算系数
# 
# 定义子目录 `RUN`，用于执行实际的 ORCA 计算。同时定义换算系数，将所得相互作用能转换为 kcal/mol，以便解释。

# In[ ]:


# > 在 `RUN` 中执行计算
working_dir = Path("RUN")
# > 自动创建或重新创建 `working_dir`
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()
# > 原子单位到 kcal/mol 的换算系数
unit_conversion = 627.509 


# ## 步骤 3：设置输入结构
# 
# 以立方烷二聚体为例分解伦敦色散（LD）相互作用。定义笛卡尔坐标下的三维结构，并将其可视化：
# 

# In[3]:


# > 将以 Angstrom 为单位的笛卡尔坐标定义为 Python 字符串
# > 原则上可替换为其他中性非共价相互作用二聚体
xyz_data = """\
32
cubane dimer
C -2.37438447 -0.01391040 0.09884370
H -1.28799640 -0.03512861 0.15212364
C -4.14778752 -0.93512832 -0.83041188
C -4.12750489 1.22660964 -0.39560995
C -3.20643343 0.25290551 -1.19176459
C -4.23643748 -0.22901780 1.25693113
C -3.31541388 -1.20282062 0.46084844
C -3.29511100 0.95905424 0.89568028
C -5.06903241 0.03775498 -0.03358762
H -4.49182155 -1.69883200 -1.52468842
H -4.45488733 2.20468479 -0.74189213
H -2.79034514 0.44620851 -2.17829474
H -4.65224637 -0.42236391 2.24356374
H -2.98771769 -2.18079631 0.80710684
H -2.95068094 1.72266937 1.58984614
H -6.15558078 0.05880427 -0.08645048
C 2.05734830 0.00538088 -0.06151765
H 0.97096030 0.02659937 -0.11479860
C 2.97807526 -0.96758450 -0.85835285
C 2.99837839 1.19429062 -0.42352221
C 3.91940243 0.22048707 -1.21960352
C 2.88939603 -0.26143458 1.22909153
C 3.81046792 -1.23513946 0.43293829
C 3.83075080 0.92659878 0.86773902
C 4.75199612 -0.04628527 0.07091613
H 2.63364560 -1.73119990 -1.55251858
H 2.67068285 2.17226622 -0.76978149
H 4.33521233 0.41383252 -2.20623583
H 2.47330675 -0.45473693 2.21562140
H 4.13784970 -2.21321450 0.77922136
H 4.17478444 1.69030273 1.56201542
H 5.83854442 -0.06733485 0.12378001

\n
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


# ## 步骤 4：对超分子体系运行 HFLD
# 
# 对完整体系运行 HFLD/def2-TZVP 计算。HFLD 的伦敦色散能依赖基组，适中大小的 def2-TZVP 能较好地兼顾计算成本与准确性。ORCA 会自动定义片段，本例为两个片段，每个单体一个。定义并运行三个辅助函数：`setup_calc`、`run_calc` 和 `check_and_parse_output`。
# 
# 首先使用 `setup_calc` 设置计算：

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 4) -> Calculator:
    # > 设置 Calculator 对象
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > 将结构赋给计算器
    calc.structure = structure

    # > 定义 ORCA 计算的简单关键词列表。
    sk_list = [
    Dlpno.HFLD,
    BasisSet.DEF2_TZVP,
    AuxBasisSet.DEF2_TZVP_C,
    Scf.NORMALSCF,
    Approximation.RIJCOSX,
    Dlpno.ADLD,
    Dlpno.LED
    ]

    # > 在计算器中使用简单关键词
    calc.input.add_simple_keywords(*sk_list)

    # > 通过设置 frag 块自动分配片段
    calc.input.add_blocks(BlockFrag(storefrags=True))

    # 定义计算使用的 CPU 数量
    calc.input.ncores = ncores # > 本次 ORCA 计算使用 4 个 CPU

    return calc

calc_supersystem = setup_calc("supersystem", working_dir=working_dir, structure=structure)


# 然后使用 `run_calc` 运行计算并获取输出：

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

output_supersystem = run_calc(calc_supersystem)


# 需要检查计算是否成功并解析输出，这由 `check_and_parse_output` 完成：

# In[6]:


def check_and_parse_output(output: Output):
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

check_and_parse_output(output_supersystem)


# ## 步骤 5：处理 ADLD 能量
# 
# 现在对 ADLD 进行后处理。首先收集伦敦色散（LD）的原子贡献。

# In[7]:


# > 二聚体（产物）中伦敦色散原子贡献的数据路径 [Eh]
product_vdw = output_supersystem.results_properties.geometries[0].mdci_adld.adlddispatomic_loewdin

# > 反应物的贡献为零，因此这些能量已经是相对值
# > 应用换算系数，将结果转换为 kcal/mol
atomic_disp: list[float] = []
for i in range(len(product_vdw)):
    product_x = float(product_vdw[i][0])
    product_x_kcal = (product_x) * unit_conversion
    atomic_disp.append(product_x_kcal)


# 现在得到各原子对相对伦敦色散相互作用能贡献的数组。可定义伦敦色散密度（$\rho_{disp}$），便于可视化和分析原子贡献：

# In[8]:


def rho_disp(r, alpha, epsilons, positions):
    """
    计算位置 r 处的色散密度 rho_disp。

    参数：
    - r: 形状为 (3,) 的类数组，待求值的三维空间点。
    - alpha: 浮点数，参数 alpha（通常设为 0.3-0.5）。
    - epsilons: 形状为 (N,) 的类数组，包含各原子 A 的色散强度。
    - positions: 形状为 (N, 3) 的类数组，包含各原子的位置 R_A。

    返回值：
    - rho: 浮点数，位置 r 处的色散密度。
    """
    # 前因子
    prefactor = (np.pi / alpha)**(-1.5)

    # 对原子求和
    r = np.array(r)
    positions = np.array(positions)
    deltas = r - positions  # 形状为 (N, 3)
    squared_distances = np.sum(deltas**2, axis=1)  # 形状为 (N,)

    contributions = epsilons * np.exp(-alpha * squared_distances)
    total_sum = np.sum(contributions)

    return prefactor * total_sum


# ## 步骤 6：ADLD 贡献的简易可视化
# 
# 原始 ADLD 论文在网格上计算 `rho_disp`，将结果写入 `.cube` 文件用于可视化。为便于在本 Notebook 内简单可视化，先仅在原子位置计算该函数，按原子中心的数值为原子着色。这些图用于定性分析，细节可能与论文中的图不同。下文也会生成用于完整可视化的 cube 文件。有兴趣的读者可参阅 GitHub 上的 [LDDsuite](https://github.com/bistonigroup/LDDsuite)。
# 
# 先在原子位置计算 $\rho_{disp}$：

# In[9]:


# > 二聚体笛卡尔坐标的数据路径
elements = output_supersystem.results_properties.geometries[0].geometry.coordinates.cartesians

# > 元素符号列表
elem_list: list[str] = []
# > 原子位置列表
positions: list[tuple[float,float,float]] = []

# > 填充列表
for elem in elements:
    elem_list.append(elem[0])
    positions.append([elem[1],elem[2],elem[3]])

# > 在原子位置计算伦敦色散密度（rho_disp）
rho_disp_atomic: list[float] = []
for i in range(len(positions)):
    # > 选择 alpha 值为 0.3
    rho_disp_atomic.append(rho_disp(positions[i], 0.3, atomic_disp, positions))


# 接下来定义一个辅助函数，将色散密度数值映射为颜色渐变：

# In[10]:


def value_to_color(value, min_val):
    """
    将 min_val 到 0 之间的数值映射到白色至红色的渐变。
    接近 0 的数值对应白色；接近 min_val 的数值对应红色。
    """
    if value >= 0:
        return '#FFFFFF'

    # 确保 min_val 为负，并避免除以零
    if min_val == 0:
        min_val = -1e-6  # 用于避免除以零的小量 epsilon

    # 归一化：从 0（白色）到 min_val（红色）
    norm = min(abs(value / min_val), 1.0)

    r = 255
    g = int(255 * (1 - norm))
    b = int(255 * (1 - norm))

    return f'#{r:02X}{g:02X}{b:02X}'


# 最后向 py3Dmol 可视化中添加彩色球体，绘制 $\rho_{disp}$。需要自定义范德华半径，以确定添加球体的大小：

# In[11]:


# > 用于可视化的自定义范德华半径
vdw_radii = {
    'H': 1.20, 'He': 1.40, 'Li': 1.82, 'Be': 1.53,
    'B': 1.92, 'C': 1.70,  'N': 1.55,  'O': 1.52,
    'F': 1.47
}

# > 定义颜色渐变的最大值和最小值
min_val = -0.01 # > kcal/mol Bohr³

# > 设置查看器
view = py3Dmol.view(width=400, height=400)
view.addModel(xyz_data)

# > 初始样式：所有原子均以白色棒状模型显示
view.setStyle({}, {'stick': {'color' : 'white'}})

# 叠加球体以实现颜色编码
for i, val in enumerate(rho_disp_atomic):
    color = value_to_color(val, min_val)
    radius = vdw_radii[elem_list[i]] * 0.35 # > 为获得合适球体大小而设置的任意缩放系数
    selection = {'serial': i}
    view.addStyle(selection, {'sphere': {'color': color, 'radius': radius}})

view.zoomTo()
view.show()


# ## 步骤 7：生成 Cube 文件
# 
# 为正确可视化，也可生成 cube 文件，再单独绘图：

# In[12]:


# > 用于获取原子序数的辅助字典
element_dict = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9
}

# > 用于定义包围盒的辅助函数
def get_bounding_box(coords, padding):
    coords = np.array(coords)
    min_coords = coords.min(axis=0) - padding
    max_coords = coords.max(axis=0) + padding
    return min_coords, max_coords

# > 用于写出 Cube 文件的辅助函数
def write_cube_file(filename, coords, atomic_numbers, values, min_coords, voxel_size, n_points):
    with open(filename, 'w') as f:
        f.write("Cube file generated by adld notebook\n")
        f.write("rho_disp\n")
        f.write(f"{len(coords):5d} {min_coords[0]:12.6f} {min_coords[1]:12.6f} {min_coords[2]:12.6f}\n")
        f.write(f"{n_points[0]:5d} {voxel_size[0]:12.6f} 0.0 0.0\n")
        f.write(f"{n_points[1]:5d} 0.0 {voxel_size[1]:12.6f} 0.0\n")
        f.write(f"{n_points[2]:5d} 0.0 0.0 {voxel_size[2]:12.6f}\n")
        for atom_num, (x, y, z) in zip(atomic_numbers, coords):
            f.write(f"{atom_num:5d} 0.0 {x:12.6f} {y:12.6f} {z:12.6f}\n")
        flat_values = values.flatten(order='F')  # Cube 文件使用 Fortran 顺序
        for i in range(0, len(flat_values), 6):
            line = "".join(f"{v:13.5e}" for v in flat_values[i:i+6])
            f.write(line + "\n")


# > 定义包围盒尺寸
padding = 3.0  # 在分子周围留出 3 Å 的边距
min_coords, max_coords = get_bounding_box(positions, padding)

# > 定义网格
n_points = (40, 40, 40)
x_lin = np.linspace(min_coords[0], max_coords[0], n_points[0])
y_lin = np.linspace(min_coords[1], max_coords[1], n_points[1])
z_lin = np.linspace(min_coords[2], max_coords[2], n_points[2])
voxel_size = (max_coords - min_coords) / (np.array(n_points) - 1)

# 在网格上求 rho_disp 的值
values = np.zeros(n_points)
for i, x in enumerate(x_lin):
    for j, y in enumerate(y_lin):
        for k, z in enumerate(z_lin):
            r = (x, y, z)
            # > 将 alpha 设为 0.3
            values[i, j, k] = rho_disp(r, 0.3, atomic_disp, positions)

# > 获取原子序数
atomic_numbers = [element_dict[symbol] for symbol in elem_list]

# 写入 Cube 文件
write_cube_file(working_dir / "rho_disp.cube", positions, atomic_numbers, values, min_coords, voxel_size, n_points)
print("Cube file written: rho_disp.cube")


# ## 小结
# 
# 本 Notebook 演示如何使用 ORCA Python 接口，将 HFLD 的伦敦色散能分解为原子贡献。在 Notebook 内进行了粗略可视化，并生成 cube 文件用于完整分析。

