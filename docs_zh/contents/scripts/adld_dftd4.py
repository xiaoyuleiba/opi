#!/usr/bin/env python
# coding: utf-8

# # ADLD：DFT-D4
# 
# 本教程演示如何将 ORCA 与 ORCA Python 接口（OPI）配合使用，通过**[伦敦色散能的原子分解（ADLD）](https://doi.org/10.1021/acscentsci.5c00356)**分析，将分子二聚体中 **[DFT-D4](https://doi.org/10.1063/1.5090222)** 的伦敦色散相互作用能分解为原子贡献。分解 DFT-D4 能量可以快速估计各原子的伦敦色散贡献。关于 LED 框架下更精细的分解，请参阅相应 Notebook。
# 
# 本 Notebook 将：
# 1. 导入所需 Python 依赖。
# 2. 定义工作目录。
# 3. 设置输入结构。
# 4. 对超体系运行 DFT-D4 计算。
# 5. 自动划分非共价子体系。
# 6. 分别计算各二聚体。
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
    Scf, Approximation, Dft, DispersionCorrection
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
# 以金刚烷二聚体为例分解伦敦色散（LD）相互作用。定义笛卡尔坐标下的三维结构，并将其可视化：
# 

# In[3]:


# > 将以 Angstrom 为单位的笛卡尔坐标定义为 Python 字符串
# > 原则上可替换为其他中性非共价相互作用二聚体
xyz_data = """\
52
adamantane dimer
C -1.2512580945 0.0014447688 -4.5331335448
C -0.0001235489 1.2445251646 -2.7612423124
C 1.2481335316 -0.0013367292 -4.5311413071
C -0.0029283957 -1.2565268787 -2.7691903511
C -1.2503012639 -1.2486921591 -3.6528854863
C -0.0006594466 0.0020950040 -5.4138613643
C -1.2467375929 1.2467239390 -3.6461299866
C 1.2473963684 1.2436834160 -3.6447126890
C 1.2453783869 -1.2511955542 -3.6514603233
C -0.0026505360 -0.0089965776 -1.8844115499
H -2.1450794262 0.0043889516 -5.1648857415
H 0.0009877220 2.1356012499 -2.1257000212
H 2.1418366957 -0.0004543894 -5.1618030462
H -0.0038482437 -2.1525328115 -2.1403679585
H 2.1455059233 -1.2737565650 -3.0272449767
H 1.2692678550 -2.1503564250 -4.2769208835
H -0.0009339824 -0.8777868157 -6.0657312229
H -1.2754300701 -2.1478464803 -4.2782930880
H -2.1504389072 -1.2698005946 -3.0287598493
H -0.8848964964 -0.0102264960 -1.2357477040
H 0.0010895729 0.8848621166 -6.0618708738
H 0.8778642353 -0.0123636823 -1.2336945592
H -2.1479497683 1.2694405004 -3.0236737797
H 1.2698733945 2.1461993692 -4.2654872765
H -1.2663652760 2.1493731136 -4.2667886778
H 2.1487325387 1.2637698441 -3.0222209102
H 3.8360448211 0.0002950370 -6.3613914539
C 4.7298358432 0.0010844825 -6.9919263799
C 4.7303085767 1.2481697816 -7.8756345239
C 4.7330471126 -1.2464687062 -7.8748194344
C 5.9783164308 0.0029481456 -6.1088865644
C 5.9773801566 1.2523639250 -8.7598945071
H 4.7076765449 2.1490254113 -7.2526031839
H 3.8288296955 1.2686817563 -8.4978754295
C 5.9807806158 -1.2475248312 -8.7579867258
H 3.8320676507 -1.2689344965 -8.4978189406
H 4.7112191662 -2.1476227034 -7.2521085138
C 7.2287109083 0.0062218143 -6.9896685254
H 5.9796434423 -0.8795548254 -5.4604383280
H 5.9753753912 0.8830178962 -5.4572372534
C 7.2245631611 1.2521204877 -7.8754596256
C 5.9782475715 0.0016408752 -9.6407077164
H 5.9761534374 2.1461192284 -9.3917461569
C 7.2279643357 -1.2408795855 -7.8736092341
H 5.9837006737 -2.1426289557 -9.3878126810
H 8.1226008273 0.0075858975 -6.3578267144
H 7.2461857452 2.1548254344 -7.2550714782
H 8.1256650763 1.2722286312 -8.4982189724
H 6.8586688519 0.0016709201 -10.2918095343
H 5.0963929708 0.0001414723 -10.2900312493
H 8.1290236854 -1.2603983485 -8.4964381158
H 7.2516713723 -2.1413794735 -7.2504037938
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


# ## 步骤 4：对超分子体系运行 DFTD4
# 
# 先对完整体系运行 B3LYP-D4/MINIX 计算。D4 伦敦色散能不依赖基组，因此选用最小基组。ORCA 会自动定义计算中的片段，本例为两个片段，每个分子一个。定义并运行三个辅助函数：`setup_calc`、`run_calc` 和 `check_and_parse_output`。
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
    Dft.B3LYP,
    DispersionCorrection.D4,
    BasisSet.MINIX,
    Scf.NORMALSCF,
    Approximation.RIJCOSX,
    Dlpno.ADLD
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


# ## 步骤 5：定义子体系
# 接下来提取自动生成的片段并定义子体系：

# In[7]:


# > 片段数据的简写变量
frags = output_supersystem.results_properties.geometries[0].geometry.fragments
# > 获取 LED 数据中的片段数量
nfrags = max(item[0] for item in frags)

# > 将自动片段划分得到的片段数据提取到字典中
# > 假定片段编号连续

# > 生成片段中的原子列表
atoms_in_frag = [[] for _ in range(nfrags)]
for index, value in enumerate(frags):
    frag_id = int(value[0])
    atoms_in_frag[frag_id-1].append(index)
    structure.atoms[index].fragment_id = frag_id

# > 定义子体系中的原子列表
atoms_in_sub = [] 
atoms_in_sub.append(atoms_in_frag[0])   # > 第一子体系
atoms_in_sub.append([atom for frag in atoms_in_frag[1:] for atom in frag])  # > 第二子体系
nsub = len(atoms_in_sub)

# > 定义子体系中的片段列表
frag_in_sub = [[] for _ in range(nsub)]
for ifrag in range(nfrags):
    if ifrag == 0: frag_in_sub[0].append(ifrag) # > 第一个片段为第一子体系
    if ifrag > 0: frag_in_sub[1].append(ifrag)  # > 第二个片段为第二子体系

# > 可视化片段
print("Here you can view the fragment IDs:")
for atom in structure.atoms:
    x, y, z = atom.coordinates.coordinates
    view.addLabel(f"{atom.fragment_id}", {
        'position': {'x': x, 'y': y, 'z': z},
        'backgroundColor': 'black',
        'fontColor': 'white',
        'showBackground': True,
        'fontSize': 12
    })
view.show()

# > 将 ID 设为零
for atom in structure.atoms:
    atom.fragment_id = None


# ## 步骤 6：运行子体系计算
# 
# 提取子体系结构，并对各子体系进行计算。使用此前定义的辅助函数：

# In[8]:


sub_outputs: list[Output] = []
# 遍历两个子体系，获取子结构并执行计算
for i in range(len(atoms_in_sub)):

    # > 准备子体系计算
    subsystem = i + 1
    sub_structure = calc_supersystem.structure.extract_substructure(atoms_in_sub[i])
    calc_subsystem = setup_calc(f"subsystem_{subsystem}", working_dir, sub_structure)

    # > 可视化即将进行计算的子体系
    # > 先移除片段数据，以便可视化
    for atom in calc_subsystem.structure.atoms: atom.fragment_id = None
    xyz_data_sub = calc_subsystem.structure.to_xyz_block()
    # > 可视化子体系
    view_sub = py3Dmol.view(width=400, height=400)
    view_sub.addModel(xyz_data_sub, 'xyz')
    view_sub.setStyle({}, {'stick': {}, 'sphere': {'scale': 0.3}})
    view_sub.zoomTo()
    print(f"Visualizing subsystem {subsystem}")
    view_sub.show()

    # > 运行计算
    output_sub = run_calc(calc_subsystem)

    # > 检查并解析输出
    check_and_parse_output(output_sub)

    sub_outputs.append(output_sub)


# ## 步骤 7：处理 ADLD 能量
# 
# 现在对 ADLD 进行后处理。收集伦敦色散（LD）的总原子贡献，并从组合后二聚体（产物）的贡献中减去独立单体（反应物）的贡献，得到相对原子贡献：

# In[9]:


# > 二聚体（产物）中伦敦色散原子总贡献的数据路径 [Eh]
product_vdw = output_supersystem.results_properties.geometries[0].vdw_correction.vdw_atomic
# > 单体（反应物）中伦敦色散原子总贡献的数据路径 [Eh]
reactant_vdw = sub_outputs[0].results_properties.geometries[0].vdw_correction.vdw_atomic + sub_outputs[1].results_properties.geometries[0].vdw_correction.vdw_atomic

# > 计算相对原子色散贡献
atomic_disp: list[float] = []
for i in range(len(product_vdw)):
    product_x = float(product_vdw[i][0])
    reactant_x = float(reactant_vdw[i][0])
    difference = (product_x - reactant_x) * unit_conversion
    atomic_disp.append(difference)


# 现在得到各原子对相对伦敦色散相互作用能贡献的数组。可定义伦敦色散密度（$\rho_{disp}$），便于可视化和分析原子贡献：

# In[10]:


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


# ## 步骤 8：ADLD 贡献的简易可视化
# 
# 原始 ADLD 论文在网格上计算 `rho_disp`，将结果写入 `.cube` 文件用于可视化。为便于在本 Notebook 内简单可视化，先仅在原子位置计算该函数，按原子中心的数值为原子着色。这些图用于定性分析，细节可能与论文中的图不同。下文也会生成用于完整可视化的 cube 文件。有兴趣的读者可参阅 GitHub 上的 [LDDsuite](https://github.com/bistonigroup/LDDsuite)。
# 
# 先在原子位置计算 $\rho_{disp}$：

# In[11]:


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

# In[12]:


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

# In[13]:


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


# ## 步骤 9：生成 Cube 文件
# 
# 为正确可视化，也可生成 cube 文件，再单独绘图：

# In[14]:


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
# 本 Notebook 演示如何使用 ORCA Python 接口，将 DFT-D4 的伦敦色散能分解为原子贡献。在 Notebook 内进行了粗略可视化，并生成 cube 文件用于完整分析。

