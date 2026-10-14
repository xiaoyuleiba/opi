#!/usr/bin/env python
# coding: utf-8

# # 计算红外光谱
# 
# 本教程演示如何使用 ORCA Python 接口（OPI）计算红外光谱。
# 
# 本 Notebook 将：
# 1. 导入所需 Python 依赖。
# 2. 定义工作目录。
# 3. 设置计算的输入结构。
# 4. 执行气相频率计算。
# 5. 解析结果并绘图。

# ## 步骤 1：导入依赖
# 
# 首先导入以下用途所需的模块：
# - 访问 ORCA 输入与输出
# - 数值计算与数据处理
# - 绘制结果
# 
# > **注意：** 我们还导入 `py3Dmol` 等可视化与绘图模块。可能需要将 `py3Dmol` 安装到 OPI 的 `venv` 中，例如激活 `.venv` 后运行 `uv pip install py3Dmol`。

# In[1]:


# > 导入 pathlib，用于目录处理
from pathlib import Path
import shutil

# > 导入 numpy，用于数据处理和数值运算
import numpy as np

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.input.simple_keywords import Dft, Task, BasisSet, Approximation, AuxBasisSet, DispersionCorrection
from opi.input.structures.structure import Structure
from opi.input.blocks.block_freq import BlockFreq

# > 导入可视化库
import matplotlib.pyplot as plt
import py3Dmol


# ## 步骤 2：定义工作目录
# 
# 所有实际计算均在子目录 **RUN_09** 中执行。

# In[ ]:


# > 在 `RUN` 中执行计算
working_dir = Path("RUN")
# > 自动创建或重新创建 `working_dir`
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()


# ## 步骤 3：设置输入结构
# 
# 以**甲醇**为示例分子。用 XYZ 格式定义笛卡尔坐标下的三维结构，并将其可视化。
# 

# In[3]:


# > 将分子以 Angstrom 为单位的笛卡尔坐标定义为 Python 字符串
xyz_data = """\
7

O            1.14430000000000        0.24120000000000        0.00000000000000
C           -1.25740000000000        0.18150000000000        0.00000000000000
C            0.11300000000000       -0.42260000000000        0.00000000000000
H           -1.79380000000000       -0.14930000000000        0.89240000000000
H           -1.18650000000000        1.27190000000000        0.00160000000000
H           -1.79280000000000       -0.14680000000000       -0.89380000000000
H            0.14780000000000       -1.52520000000000       -0.00070000000000\n
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


# ## 步骤 4：气相计算
# 
# 
# 执行几何优化和频率计算

# In[4]:


# > 设置隐式模型计算的路径并创建目录
xyz_file = working_dir / "struc.xyz"

# > 创建用于生成 ORCA 输入并执行计算的 Calculator 对象
calc = Calculator(basename="gas", working_dir=working_dir)

# > 从 XYZ 文件加载分子结构
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 0
calc.structure.multiplicity = 1

# > 添加计算关键词
calc.input.add_simple_keywords(
    Dft.B3LYP, # > B3LYP 方法
    BasisSet.AUG_CC_PVDZ, # > daug-cc-pVDZ 基组
    AuxBasisSet.DEF2_J, # > 辅助基组
    Approximation.RIJCOSX, # > 使用 RIJ 和 COSX
    DispersionCorrection.D4, # > 使用 D4 模型
    Task.FREQ # > 计算解析 Hessian 矩阵
)

calc.input.add_blocks(
    BlockFreq(scalfreq=0.970) # > B3LYP/aug-cc-pVDZ 的最佳缩放系数
)

# > 写出 ORCA 输入文件
calc.write_input()
# > 运行 ORCA 计算
print("Running ORCA frequency calculation ...", end="")
calc.run()
print("   Done")


# ## 步骤 5：解析振动光谱并绘图
# 
# 定义一个函数，从 ORCA 中提取振动光谱并绘图。

# In[5]:


def plot_vibrational_spectrum(vibrational_spectrum: Path) -> None:
    wave_numbers = []
    intensities = []

    # > OPI 已能提供频率，
    # > 但目前尚不能提供强度。
    # > 因此在此手动检索并提取。
    with open(vibrational_spectrum) as f:
        in_block = False
        for line in f:
            if not in_block and "IR SPECTRUM" in line:
                in_block = True
                for _ in range(5):
                    next(f)
            elif in_block and not line.strip():  # 块结束
                break
            elif in_block:
                parts = line.split()
                try:
                    wave_number = float(parts[1])
                    intensity = float(parts[3])
                    if wave_number > 0 and intensity > 0:
                        wave_numbers.append(wave_number)
                        intensities.append(intensity)
                except ValueError:
                    continue

    # > 创建零基线并填入强度值，以生成红外光谱
    max_wave = int(max(wave_numbers)) + 100
    full_range = np.arange(0, max_wave, 1)
    spectrum = np.zeros_like(full_range, dtype=float)

    for wn, inten in zip(wave_numbers, intensities):
        index = int(wn)
        if index < len(spectrum):
            spectrum[index] = inten

    # > 绘制振动红外光谱
    plt.figure(figsize=(10, 5))
    plt.plot(full_range, spectrum)
    plt.xlabel('Wavenumber (cm⁻¹)')
    plt.ylabel('IR Intensity (km/mol)')
    plt.title('Vibrational Spectrum')
    plt.grid(True)
    plt.xlim(0, max_wave)
    plt.ylim(0, max(intensities))
    plt.tight_layout()
    plt.show()

# > 绘制隐式溶剂化结果的振动光谱
spectrum_path = working_dir / "gas.out"
plot_vibrational_spectrum(spectrum_path)


# ## 小结
# 
# 本 Notebook 演示如何使用 ORCA Python 接口计算预测红外光谱所需的频率与强度。

