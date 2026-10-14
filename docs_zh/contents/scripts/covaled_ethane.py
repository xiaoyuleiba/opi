#!/usr/bin/env python
# coding: utf-8

# # CovaLED 分析
# 
# 本教程演示如何将 ORCA 与 ORCA Python 接口（OPI）配合使用，通过**[局域能量分解（LED）](https://doi.org/10.1021/acs.jctc.6b00523)**分析，分解非共价主客体复合物的相互作用能。CovaLED 还允许将子体系划分为由共价键连接的片段。此外，还进行**[片段对 LED（fp-LED）](https://doi.org/10.1002/anie.202421922)**后处理，将相互作用能完全分解为两两片段间的贡献。
# 
# 本 Notebook 将：
# 1. 导入所需 Python 依赖。
# 2. 定义工作目录。
# 3. 设置主客体输入结构。
# 4. 对完整体系运行 LED 计算。
# 5. 自动划分主体与客体。
# 6. 运行主体与客体计算。
# 7. 对标准 LED 进行后处理。
# 8. 对 fp-LED 进行后处理。
# 9. 将结果整理为 pandas DataFrame。
# 10. 使用热图可视化 LED 相互作用。
# 11. 解释结果。

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
import pandas as pd
import numpy as np

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, Dlpno, \
    Scf, AuxBasisSet, Approximation, Wft
from opi.input.structures.structure import Structure

# > 导入用于结果绘图和分子可视化的库
import matplotlib.pyplot as plt
import seaborn as sns
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
# 以简单的乙烷二聚体为例。定义笛卡尔坐标下的三维结构，并将其可视化：
# 

# In[3]:


# > 将笛卡尔坐标定义为 Python 字符串
xyz_data = """\
16

   C   -7.95601053588590      3.15195494946836      0.27839423386884
   H   -8.05179200418556      2.32105035333223      0.97643318840757
   H   -8.76671687926839      3.06907938063429     -0.44469508271450
   H   -8.10439460873783      4.07476211185951      0.83796813304307
   C   -6.59679579710356      3.13992001745864     -0.41129316526153
   H   -6.44433232188420      2.21634557454439     -0.96856734466103
   H   -5.78737090643299      3.22474156709677      0.31305681682191
   H   -6.50317719735328      3.96948543870281     -1.11109265496380
   C   -7.98491788051681     -0.47692316195123     -0.27920144845602
   H   -8.79547255792874     -0.40025515615629      0.44475840345098
   H   -8.14024532507479     -1.38893575993245     -0.85442072275227
   H   -8.07387099754762      0.36622931992467     -0.96327738147245
   C   -6.62637230870360     -0.48881716933962      0.41184856873533
   H   -6.54173070676907     -1.33138502402438      1.09714355762301
   H   -6.46548339397797      0.42323949114744      0.98551040416686
   H   -5.81650657862927     -0.56907193276525     -0.31256550583595\n
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


# ## 步骤 4：对超分子体系运行 CovaLED
# 
# 先对完整体系进行 DLPNO-CCSD(T)/def2-SVP 的 LED 计算，采用 `NORMALPNO` 和 `NORMALSCF` 设置。为降低本 Notebook 的计算成本，选用小基组和适中的 PNO、SCF 设置。注意：这些设置通常不足以满足正式计算的要求。ORCA 会自动检测 CovaLED 计算中的片段，每个乙烷分子分为两个 CH<sub>3</sub> 片段，共四个；也可手动指定。定义并运行三个辅助函数：`setup_calc`、`run_calc` 和 `check_and_parse_output`。
# 
# 首先使用 `setup_calc` 设置计算：

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 4) -> Calculator:
    # > 设置 Calculator 对象，将 ORCA 计算的基本名称设为 LED
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > 将结构赋给计算器
    calc.structure = structure

    # > 定义 ORCA 计算的简单关键词列表
    sk_list = [
    Wft.DLPNO_CCSD_T, # > 使用 DLPNO-CCSD(T)
    BasisSet.DEF2_SVP, # > 使用 def2-SVP 基组
    AuxBasisSet.DEF2_SVP_C, # > 使用 def2-svp/c 辅助基组
    AuxBasisSet.DEF2_J, # > 使用 def2-svp/j 辅助基组
    Dlpno.NORMALPNO, # > Normal PNO 设置
    Scf.NORMALSCF, # > Normal SCF 设置
    Approximation.RIJCOSX, # > 使用 RIJCOSX
    Dlpno.LED, # > 进行能量分解
    ]

    # 在计算器中使用简单关键词
    calc.input.add_simple_keywords(*sk_list)

    # > 如需指定片段，可在此赋值，
    # 否则由 ORCA 的自动片段划分程序生成。
    calc.input.add_arbitrary_string("%frag\n STOREFRAGS true\n FragProc  FUNCTIONALGROUPS,EXTEND,CONNECTIVITY\n DOINTERFRAGBONDS true\n end")

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
# 提取自动生成的片段，定义主体和客体子体系。本例中，主体和客体各为一个含两个片段的乙烷分子。

# In[7]:


# > 片段数据的简写变量
frags = output_supersystem.results_properties.geometries[0].geometry.fragments
# > 获取超体系 LED 数据中的片段数量（本例为 4）
nfrags = int(output_supersystem.results_properties.geometries[0].mdci_led.numoffragments)

# > 客体的第一个片段
ofs_sub = 3

# > 生成片段中的原子列表
atoms_in_frag = [[] for _ in range(nfrags)]
for index, value in enumerate(frags):
    frag_id = int(value[0])
    # > ORCA 内部编号从 0 开始，因此偏移 1
    atoms_in_frag[frag_id-1].append(index)
    structure.atoms[index].fragment_id = frag_id

# > 定义子体系中的原子列表
atoms_in_sub = [] 
atoms_in_sub.append([atom for frag in atoms_in_frag[:ofs_sub-1] for atom in frag]) # > 第一子体系包含第一个乙烷分子的两个片段
atoms_in_sub.append([atom for frag in atoms_in_frag[ofs_sub-1:] for atom in frag]) # > 第二子体系包含第二个乙烷分子的两个片段
nsub = len(atoms_in_sub)

# > 定义子体系中的片段列表
frag_in_sub = [[] for _ in range(nsub)]
for ifrag in range(nfrags):
    if ifrag < ofs_sub-1: frag_in_sub[0].append(ifrag)
    if ifrag >= ofs_sub-1: frag_in_sub[1].append(ifrag)

# > 定义客体子体系的索引偏移量
offset = len(frag_in_sub[0])

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


# ## 步骤 6：运行主体与客体计算

# In[8]:


sub_outputs: list[Output] = []
# > 遍历子体系，获取子结构并执行计算
for i in range(len(atoms_in_sub)):

    # > 准备子体系计算
    subsystem = i + 1
    sub_structure = calc_supersystem.structure.extract_substructure(atoms_in_sub[i])
    calc_subsystem = setup_calc(f"subsystem_{subsystem}", working_dir, sub_structure)

    # > 可视化即将进行计算的子体系
    # > 先移除片段数据，以便可视化
    for atom in calc_subsystem.structure.atoms: 
        atom.fragment_id = None

    # > 运行计算
    output_sub = run_calc(calc_subsystem)

    # > 检查并解析输出
    check_and_parse_output(output_sub)

    sub_outputs.append(output_sub)


# ## 步骤 7：标准 LED 的后处理
# 
# 现在对 LED 进行后处理。

# In[9]:


# > 超体系 LED 数据的简写变量
led_data = output_supersystem.results_properties.geometries[0].mdci_led
# > 子体系数量为两个
nsub = 2
# > 片段数量
nfrags = int(output_supersystem.results_properties.geometries[0].mdci_led.numoffragments)

# > 为子体系（主体和客体）的 CCSD(T) 能量创建简写变量
sub_outputs_energy = []
sub_outputs_led = []
for output in sub_outputs:
    # > energy[0] 为 SCF 能量，energy[1] 为 CCSD(T) 能量
    sub_outputs_energy.append(output.results_properties.geometries[0].energy[1])
    sub_outputs_led.append(output.results_properties.geometries[0].mdci_led)

# > 处理 LED 数据
LED_pairwise_Total = {}
LED_pairwise_Ref = {}
LED_pairwise_Corr = {}
LED_pairwise_Ref_Coulomb = {}
LED_pairwise_Ref_Exch = {}
LED_pairwise_Corr_Disp = {}

# > 首先收集超体系 LED 计算的结果
for x in range(nfrags):
    for y in range(x + 1):
        pair = (x, y)
        E_tot = led_data.totint[x][y]
        E_ref = led_data.refint[x][y]
        E_corr = led_data.corrint[x][y]
        LED_pairwise_Ref[pair] = E_ref
        LED_pairwise_Corr[pair] = E_corr
        LED_pairwise_Total[pair] = E_tot
        if x != y:
            E_ref_elstat = led_data.electrostref[x][y]
            E_ref_exch = led_data.exchangeref[x][y]
            E_corr_disp = led_data.dispcontr[x][y] + led_data.dispweak[x][y]
            LED_pairwise_Ref_Coulomb[pair] = E_ref_elstat
            LED_pairwise_Ref_Exch[pair] = E_ref_exch
            LED_pairwise_Corr_Disp[pair] = E_corr_disp

# > 现在减去子体系 LED 计算的结果
for x in range(nfrags):
    for y in range(x + 1):
        pair = (x, y)
        E_ref = 0.
        E_corr = 0.
        foundPair = False
        # > 检查片段是否属于子体系 1
        if x in frag_in_sub[0] and y in frag_in_sub[0]:
            foundPair = True
            if sub_outputs_led[0]: 
                E_ref = sub_outputs_led[0].refint[x][y]
                E_corr = sub_outputs_led[0].corrint[x][y]
                if x != y: 
                    E_ref_elstat = sub_outputs_led[0].electrostref[x][y]
                    E_ref_exch = sub_outputs_led[0].exchangeref[x][y]
                    E_corr_disp = sub_outputs_led[0].dispcontr[x][y] + sub_outputs_led[0].dispweak[x][y]
            else:
                E_ref = sub_outputs_energy[0].refenergy[0][0]
                E_corr = sub_outputs_energy[0].correnergy[0][0]
        # > 检查片段是否属于子体系 2
        elif x in frag_in_sub[1] and y in frag_in_sub[1]:
            foundPair = True
            if sub_outputs_led[1]:
                E_ref = sub_outputs_led[1].refint[x-offset][y-offset]
                E_corr = sub_outputs_led[1].corrint[x-offset][y-offset]
                if x != y and sub_outputs_led[1]:
                    E_ref_elstat = sub_outputs_led[1].electrostref[x-offset][y-offset]
                    E_ref_exch = sub_outputs_led[1].exchangeref[x-offset][y-offset]
                    E_corr_disp = sub_outputs_led[1].dispcontr[x-offset][y-offset] + sub_outputs_led[1].dispweak[x-offset][y-offset]
            else:
                E_ref = sub_outputs_energy[1].refenergy[0][0]
                E_corr = sub_outputs_energy[1].correnergy[0][0]

        # > 减去子体系的计算数据
        if foundPair:
            LED_pairwise_Ref[pair] -= E_ref
            LED_pairwise_Corr[pair] -= E_corr
            LED_pairwise_Total[pair] -= (E_ref + E_corr)
            if x != y:
                LED_pairwise_Ref_Coulomb[pair] -= E_ref_elstat
                LED_pairwise_Ref_Exch[pair] -= E_ref_exch
                LED_pairwise_Corr_Disp[pair] -= E_corr_disp


# ## 步骤 8：fp-CovaLED 的后处理
# 
# 接下来对 CovaLED 结果进行后处理，得到 fp-CovaLED 结果，将相互作用能完全分解为两两贡献。为改善数值稳定性，根据绝对值计算缩放因子。分别对总能量和参考能量进行缩放，再由总能量减去参考能量，得到相关贡献。

# In[10]:


def fp_analysis_abs():
    for x in range(nfrags):
        for y in range(x):
            pairXY = (x,y)
            pairXX = (x,x)
            pairYY = (y,y)
            E_Total_XY_unscaled = LED_pairwise_Total[pairXY]
            dE_Total_X_unscaled = LED_pairwise_Total[pairXX]
            dE_Total_Y_unscaled = LED_pairwise_Total[pairYY]
            E_Ref_XY_unscaled = LED_pairwise_Ref[pairXY]
            dE_Ref_X_unscaled = LED_pairwise_Ref[pairXX]
            dE_Ref_Y_unscaled = LED_pairwise_Ref[pairYY]
            denom_X_Total = 0
            denom_Y_Total = 0
            denom_X_Ref = 0
            denom_Y_Ref = 0
            # > 使用绝对值计算权重的分母
            for k in range(nfrags):
                minXK = min(x, k)
                maxXK = max(x, k)
                pairXK = (maxXK, minXK)
                E_Total_XK_unscaled = LED_pairwise_Total[pairXK]
                E_Ref_XK_unscaled = LED_pairwise_Ref[pairXK]
                if k != x:
                    denom_X_Total += abs(E_Total_XK_unscaled)
                    denom_X_Ref += abs(E_Ref_XK_unscaled)
                minYK = min(y, k)
                maxYK = max(y, k)
                pairYK = (maxYK, minYK)
                E_Total_YK_unscaled = LED_pairwise_Total[pairYK]
                E_Ref_YK_unscaled = LED_pairwise_Ref[pairYK]
                if k != y:
                    denom_Y_Total += abs(E_Total_YK_unscaled)
                    denom_Y_Ref += abs(E_Ref_YK_unscaled)
            scalingX_Total = dE_Total_X_unscaled / denom_X_Total
            scalingY_Total = dE_Total_Y_unscaled / denom_Y_Total
            scaling_Total = scalingX_Total + scalingY_Total
            scalingX_Ref = dE_Ref_X_unscaled / denom_X_Ref
            scalingY_Ref = dE_Ref_Y_unscaled / denom_Y_Ref
            scaling_Ref = scalingX_Ref + scalingY_Ref
            # > 取片段间贡献的绝对值
            E_Total_XY_scaled = scaling_Total * abs(E_Total_XY_unscaled)
            E_Ref_XY_scaled = scaling_Ref * abs(E_Ref_XY_unscaled)
            E_Corr_XY_scaled = E_Total_XY_scaled - E_Ref_XY_scaled
            # > 计算 Elprep 贡献
            LED_pairwise_Ref_Elprep[pairXY] = E_Ref_XY_scaled
            LED_pairwise_Corr_Elprep[pairXY] = E_Corr_XY_scaled
            LED_pairwise_Total_Elprep[pairXY] = E_Total_XY_scaled

    # > 加入 Elprep 贡献
    for x in range(nfrags):
        for y in range(x):
            pairXY = (x, y)
            LED_pairwise_Ref[pairXY] += LED_pairwise_Ref_Elprep[pairXY]
            LED_pairwise_Corr[pairXY] += LED_pairwise_Corr_Elprep[pairXY]
            LED_pairwise_Total[pairXY] += LED_pairwise_Total_Elprep[pairXY]

    # > 将已分配出去的片段内贡献设为零
    for x in range(nfrags):
        pair = (x,x)
        LED_pairwise_Total[pair] = 0.
        LED_pairwise_Ref[pair] = 0.
        LED_pairwise_Corr[pair] = 0.
        LED_pairwise_Ref_Coulomb[pair] = 0.
        LED_pairwise_Ref_Exch[pair] = 0.
        LED_pairwise_Corr_Disp[pair] = 0.

LED_pairwise_Total_Elprep = {}
LED_pairwise_Ref_Elprep = {}
LED_pairwise_Corr_Elprep = {}
# > 运行实际分析
fp_analysis_abs()


# ## 步骤 9：将结果汇集到 pandas 数据框
# 
# 为方便访问，将结果保存为由 pandas DataFrame 组成的字典。

# In[11]:


# > 生成 pandas 数据框的函数
def gen_led_frame(led_dict: dict[tuple[int,int], float]) -> pd.DataFrame:
    # 创建一个全部填充 NaN 的二维数组
    grid = np.full((nfrags, nfrags), np.nan)

    # > 用 LED 结果填充网格
    for (x, y), value in led_dict.items():
        grid[y, x] = value * unit_conversion

    # > 转换为 DataFrame
    df = pd.DataFrame(grid)
    # > 索引偏移 1，使编号从 1 开始
    df.index += 1
    df.columns += 1

    return df

# > 生成数据框
data_frames = {}
data_frames["LED_Total"] = gen_led_frame(LED_pairwise_Total)
data_frames["LED_Corr"] = gen_led_frame(LED_pairwise_Corr)
data_frames["LED_Ref"] = gen_led_frame(LED_pairwise_Ref)
data_frames["LED_Ref_Coulomb"] = gen_led_frame(LED_pairwise_Ref_Coulomb)
data_frames["LED_Ref_Exch"] = gen_led_frame(LED_pairwise_Ref_Exch)
data_frames["LED_Corr_Disp"] = gen_led_frame(LED_pairwise_Corr_Disp)
data_frames["LED_Total_Elprep"] = gen_led_frame(LED_pairwise_Total_Elprep)
data_frames["LED_Ref_Elprep"] = gen_led_frame(LED_pairwise_Ref_Elprep)
data_frames["LED_Corr_Elprep"] = gen_led_frame(LED_pairwise_Corr_Elprep)


# ## 步骤 10：用热图可视化结果
# 
# 使用热图可视化 fp-LED 总相互作用能：

# In[12]:


# 将结果绘制为热图的函数
def plot_led(name: str, df: pd.DataFrame):
    plt.figure(figsize=(6, 5))
    sns.heatmap(df, annot=True,fmt=".1f", cmap="RdBu", center=0)
    plt.title(f"{name} (kcal/mol)")
    plt.xlabel("Fragment X")
    plt.ylabel("Fragment Y")
    plt.show()

# > 再次显示带有片段 ID 的结构
print("Visualization of fragments:")
view.show()

# > 计算片段 1 与其余所有片段之间的总相互作用能
etot = np.nansum(data_frames["LED_Total"].to_numpy())
print(f"Total interaction energy is {etot:.1f} kcal/mol")

# > 绘制热图
plot_led("LED_Total",data_frames["LED_Total"])


# ## 步骤 11：结果解读
# 上面的热图显示片段 1 与其余部分相互作用能的两两贡献。总相互作用能为负。也可绘制不同贡献：

# In[13]:


plot_led("LED_Ref",data_frames["LED_Ref"])
plot_led("LED_Corr",data_frames["LED_Corr"])


# 可以看到，相关效应对于总相互作用能为负至关重要。
# > **注意：** 如需查看其他贡献，可下载本 Notebook，在本地运行，并使用这里提供的功能绘制其他热图。

# ## 小结
# 将 fp-CovaLED 方案应用于乙烷二聚体，演示如何使用 ORCA Python 接口（OPI）执行 CovaLED 计算。它可用于按片段研究蛋白质与配体相互作用等问题。原则上，自动片段划分器可用于任何结构，但应确保子体系之间没有共价键。

