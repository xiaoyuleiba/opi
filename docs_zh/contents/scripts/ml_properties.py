#!/usr/bin/env python
# coding: utf-8

# # 用于机器学习模型的性质
# 
# 本教程演示如何结合 ORCA 与 ORCA Python 接口（OPI）计算电子性质，例如用于训练机器学习原子间势（MLIP）的性质。示例采用 &omega;B97M-V/def2-TZVPD 参考水平，该水平可用于生成可靠的训练数据，例如 [OMol25](https://arxiv.org/abs/2505.08762)。示例分子为 4-甲基吡唑。
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
# 首先导入以下用途所需的模块：
# - 访问 ORCA 输入与输出
# - 数值计算与数据处理
# - 绘制结果
# 

# In[1]:


from pathlib import Path
import shutil

# > 使用 pandas 和 numpy 处理数据
import pandas as pd
import numpy as np

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, \
    Scf, AuxBasisSet, Approximation, Dft, Task, OutputControl, Grid
from opi.input.structures.structure import Structure

# > 用于分子可视化
import py3Dmol


# ## 步骤 2：工作目录
# 
# 定义子目录 `RUN`，用于执行实际的 ORCA 计算。

# In[ ]:


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
12

N  1.3237  0.6845 -0.0002
N  1.4406 -0.6544 -0.0001
C -0.7492 -0.0043  0.0001
C  0.0327  1.1249  0.0002
C -2.2189 -0.0810 -0.0002
C  0.1711 -1.0697  0.0001
H -0.2057  2.1791  0.0003
H -2.6522  0.6920 -0.6433
H -2.5652 -1.0527 -0.3670
H -2.6091  0.0584  1.0131
H -0.0275 -2.1321  0.0003
H  2.1699  1.2401 -0.0003\n
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
# 接下来定义 ORCA 计算使用的设置。可以创建 calculator 对象，并向其添加简单关键词：

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 4) -> Calculator:
    # > 设置 Calculator 对象，将 ORCA 计算的基本名称设为 LED
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > 将结构赋给计算器
    calc.structure = structure

    # 使用简单关键词定义理论水平：wB97M-V/def2-TZVPPD
    sk_list = [
        Dft.WB97M_V, # > wB97M-V 方法
        BasisSet.DEF2_TZVPD, # > def2-TZVPD 基组
        AuxBasisSet.DEF2_J, # > 辅助基组
        Scf.TIGHTSCF, # > 使用严格的 SCF 收敛设置
        Approximation.RIJCOSX, # > 使用 RIJ 和 COSX
        Grid.DEFGRID3, # > 交换相关和 COSX 使用较密的积分网格
        Task.ENGRAD, # > 计算能量和梯度
        OutputControl.PRINTGAP, # > 输出 HOMO-LUMO 能隙
    ]

    # 在计算器中使用简单关键词
    calc.input.add_simple_keywords(*sk_list)

    # 定义计算使用的 CPU 数量
    calc.input.ncores = ncores # > 本次 ORCA 计算使用 4 个 CPU

    return calc

calc = setup_calc("single_point", working_dir=working_dir, structure=structure)


# ## 步骤 5：执行 ORCA 计算
# 
# 先将输入写成 ORCA 格式，再使用 calculator 类的 `run()` 函数运行作业。可通过其 `get_output()` 函数将 ORCA 输出写入 calculator 对象。

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
# 执行自动化计算时，应始终检查正常终止与收敛；计算较大数据集时，还应处理收敛失败的情况。简单检查示例如下：

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


# ## 步骤 7：处理结果
# 计算成功后，可以访问并保存结果。例如：

# In[7]:


def process_output(output: Output) -> dict:
    # > Hartree 到 eV 的换算系数
    h_to_eV = 27.2114
    # > Bohr 到 Angstrom 的换算系数
    bohr_to_A = 0.529177

    # > 能量和梯度信息
    # > 以 eV 为单位的总能量
    Etot = output.results_properties.geometries[0].single_point_data.finalenergy*h_to_eV
    # > 以 eV 为单位的 VV10 能量部分
    # > 注意，VV10 在梯度计算中自洽处理，因此使用 ecnl
    E_VV10 = output.results_properties.geometries[0].dft_energy.ecnl*h_to_eV
    # > 以 eV/Bohr 为单位的梯度（长度为 n_atoms 的列表，每个原子对应一个含三个元素的列表）
    Grad = output.results_properties.geometries[0].nuclear_gradient[0].grad
    Grad = [[x * (h_to_eV/bohr_to_A) for x in nuc_grad] for nuc_grad in Grad]

    # > 一般信息
    # > 分子电荷（浮点数）
    tot_chrg = output.results_properties.calculation_info.charge
    # > 分子自旋（浮点数）
    tot_spin = output.results_properties.calculation_info.mult
    # > 原子数（浮点数）
    n_at = output.results_properties.calculation_info.numofatoms
    # > 电子数（浮点数）
    n_el = output.results_properties.calculation_info.numofelectrons
    # > 基函数数目（浮点数）
    n_bf = output.results_properties.calculation_info.numofbasisfuncts
    # > alpha 电子数（浮点数）
    n_el_alp = output.results_properties.geometries[0].dft_energy.nalphael
    # > beta 电子数（浮点数）
    n_el_beta = output.results_properties.geometries[0].dft_energy.nbetael

    # > 布居分析
    # > Loewdin 电荷（长度为 n_atoms 的列表，其中每项为含 3 个元素的列表）
    chrgs_lowedin = output.results_properties.geometries[0].loewdin_population_analysis[0].atomiccharges
    # > Mulliken 电荷（长度为 n_atoms 的列表，其中每项为含 3 个元素的列表）
    chrgs_mulliken = output.results_properties.geometries[0].mulliken_population_analysis[0].atomiccharges
    # > Mayer 键级（列表长度等于键的数量）
    bond_order = output.results_properties.geometries[0].mayer_population_analysis[0].bondorders
    # > Mayer 总价（长度为 n_atoms 的列表）
    va_mayer = output.results_properties.geometries[0].mayer_population_analysis[0].va

    # > 偶极矩
    # > 总偶极矩（包含 x、y、z 三个分量的列表）
    dip_tot = output.results_properties.geometries[0].dipole_moment[0].dipoletotal
    # > 偶极矩大小（浮点数）
    dip_mag = output.results_properties.geometries[0].dipole_moment[0].dipolemagnitude
    # > 电子偶极矩贡献（包含 x、y、z 三个分量的列表）
    dip_el = output.results_properties.geometries[0].dipole_moment[0].dipoleeleccontrib
    # > 原子核偶极矩贡献（包含 x、y、z 三个分量的列表）
    dip_nu = output.results_properties.geometries[0].dipole_moment[0].dipolenuccontrib

    # > 分子轨道处理
    mos = output.results_gbw.molecule.molecularorbitals.mos # > 保存一个长度为 n_atoms 的分子轨道列表
    mo_energies = []
    mo_occupancy = []
    for mo in mos:
        mo_energies.append(mo.orbitalenergy) # > 提取分子轨道能量
        mo_occupancy.append(mo.occupancy) # > 提取占据数
    mo_energies = np.array(mo_energies) # > 转换为浮点数列表
    mo_occupancy = np.array(mo_occupancy) # > 转换为浮点数列表
    if np.all(np.isin(mo_occupancy, [1.0, 0.0])):
        print("Unrestricted calculation detected.")
    elif np.all(np.isin(mo_occupancy, [2.0, 0.0])):
        print("Restricted calculation detected.")

    HOMO_energy = np.max(mo_energies[(mo_occupancy == 1.0) | (mo_occupancy == 2.0)]) # > 确定 HOMO 能量
    HOMO_energy = HOMO_energy * h_to_eV
    LUMO_energy = np.min(mo_energies[mo_occupancy == 0.0]) # > 确定 LUMO 能量
    LUMO_energy = LUMO_energy * h_to_eV

    summary_data = {
    "Total Energy [eV]": [Etot],
    "Gradient [ev/A]": [str(Grad[0][0])[:10]+', '+str(Grad[1][0])[:10]+', '+str(Grad[2][0])[:10]],
    "VV10 Energy [eV]": [E_VV10],
    "Dipole Magnitude": [dip_mag],
    "Total Charge": [tot_chrg],
    "Spin Multiplicity": [tot_spin],
    "Number of Atoms": [n_at],
    "Number of Electrons": [n_el],
    "Alpha Electrons": [n_el_alp],
    "Beta Electrons": [n_el_beta],
    "Basis Functions": [n_bf],
    "HOMO Energy [eV]": [HOMO_energy],
    "LUMO Energy [eV]": [LUMO_energy],
    }

    return summary_data

summary_data = process_output(output)


# 下面使用 pandas 显示所得数值：

# In[8]:


df_summary = pd.DataFrame(summary_data, index=[""])
df_summary.T


