#!/usr/bin/env python
# coding: utf-8

# # Fock 矩阵对角化
# 
# 本教程演示如何使用 ORCA Python 接口（OPI）从 ORCA 计算中提取 Fock 矩阵。提取后将 Fock 矩阵对角化，并利用 Koopmans 定理估算电离能。
# 
# 本 Notebook 将：
# 1. 导入所需依赖。
# 2. 定义工作目录并准备结构。
# 3. 从 PubChem 下载并准备分子结构
# 4. 运行 ORCA 计算。
# 5. 提取矩阵并对 Fock 矩阵进行对角化

# ## 步骤 1：导入依赖
# 
# 首先导入本教程所需的 Python 模块，包括：
# 
# - **文件与目录处理**：管理输入和输出文件。
# - **JSON 解析**：读取 ORCA 生成的 GBW JSON 文件。
# - **数值计算库**：执行矩阵运算、对角化和本征值计算。
# - **PubChem 数据获取**：通过 REST API 直接从 PubChem 下载分子结构。
# - **RDKit**：处理和转换分子，并生成三维几何结构。
# - **ORCA Python 接口（OPI）**：设置、运行并解析 ORCA 量子化学计算。
# - **py3Dmol**：直接在 Notebook 中交互式地显示分子三维结构。
# 
# > **注意：** 我们还导入 `py3Dmol` 等可视化与绘图模块。可能需要将 `py3Dmol` 安装到 OPI 的 `venv` 中，例如激活 `.venv` 后运行 `uv pip install py3Dmol`。

# In[1]:


# > 导入 pathlib，用于目录处理
from pathlib import Path
import shutil

# > 导入用于解析和处理 JSON 数据的库
import json

# > 导入用于数据处理和数值运算的库
import numpy as np
from scipy.linalg import eigh

# > 导入用于从 PubChem 获取分子数据的库
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

# > 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import BasisSet, Scf
from opi.input.simple_keywords.method import Method

# > 导入 py3Dmol，用于分子的三维可视化
import py3Dmol


# ## 步骤 2：定义工作目录
# 
# 所有实际计算均在子目录 `RUN` 中执行。

# In[ ]:


# > 在 `RUN` 中执行计算
working_dir = Path("RUN")
# > 自动创建或重新创建 `working_dir`
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()
# > 原子单位到 eV 的换算系数
unit_conversion = 27.211386


# ## 步骤 3：从 PubChem 下载并准备分子结构
# 
# 本步骤通过 PubChem 的 REST API 直接获取目标分子的结构。
# 

# In[3]:


# > 从 PubChem 下载分子数据
mol_name = "Benzene"
cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{mol_name}/cids/JSON"
cid_response = requests.get(cid_url)
cid_response.raise_for_status()
cid = cid_response.json()['IdentifierList']['CID'][0]

sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
sdf_response = requests.get(sdf_url)
sdf_response.raise_for_status()
sdf_data = sdf_response.text

# > 解析 SDF，并用 RDKit 生成三维构象
mol = Chem.MolFromMolBlock(sdf_data)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()

# > 转换为 XYZ 格式
num_atoms = mol.GetNumAtoms()
lines = [str(num_atoms), f"{mol_name} from PubChem"]
for atom in mol.GetAtoms():
    pos = conf.GetAtomPosition(atom.GetIdx())
    symbol = atom.GetSymbol()
    lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

xyz_data = "\n".join(lines)

# > 使用 py3Dmol 可视化分子结构
view = py3Dmol.view(width=500, height=500)
view.addModel(xyz_data, "xyz")
view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})
view.setBackgroundColor('white')
view.zoomTo()
view.show()

# 将 XYZ 结构保存到文件
with open(working_dir / "struc.xyz","w") as f:
    f.write(xyz_data)


# ## 步骤 4：ORCA 计算
# 
# 现在执行 HF/STO-3G 计算。设置 config_dict，并用它请求输出积分。

# In[4]:


# > 设置结构文件路径
xyz_file = working_dir / "struc.xyz"

# > 定义 GBW JSON 配置
config_dict = {
    "MOCoefficients": True,
    "1elPropertyIntegrals": ["dipole"],
    "1elIntegrals": ["H", "S"],
    "FockMatrix": ["F", "J", "K"],
    "JSONFormats": ["json"]
}

# > 创建用于生成 ORCA 输入并执行计算的 Calculator 对象
calc = Calculator(basename="benzene", working_dir=working_dir)

# > 从 XYZ 文件加载分子结构
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 0
calc.structure.multiplicity = 1

# > 添加计算关键词
calc.input.add_simple_keywords(
            Method.HF,
            BasisSet.STO_3G,
            Scf.TIGHTSCF
)

# > 写出 ORCA 输入文件
calc.write_input()
# > 运行 ORCA 计算
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")

# > 获取输出，并按配置生成 GBW JSON 输出
output = calc.get_output()
status = output.terminated_normally()
if status: 
    output.parse()
    output.create_gbw_json(force=True,config=config_dict)
else:
    raise RuntimeError("ORCA did not terminate normally.")


# ## 步骤 5：提取矩阵并对 Fock 矩阵进行对角化

# In[5]:


# 加载 GBW JSON 文件
with open(working_dir / "benzene.json") as f:
    json_data = json.load(f)

# 从 JSON 文件中提取积分：
h = np.array(json_data['Molecule']['H-Matrix'])
f = np.array(json_data['Molecule']['F-Matrix'][0])
s = np.array(json_data['Molecule']['S-Matrix'])

# 将芯哈密顿量与 Fock 修正相加，构建完整的 AO Fock 矩阵
fao = h + f
# 求解广义本征值问题
eigenvalues, _ = eigh(fao, s)
# 将轨道能量从 Hartree 转换为 eV
eigenvalues_ev = eigenvalues * unit_conversion

# 使用 Koopmans 定理估算电离能
print("Ionization energy (eV):", -max([e for e in eigenvalues_ev if e < 0]))


# ## 小结
# 
# 本教程演示了如何使用 ORCA Python 接口（OPI）获取 ORCA 计算中的积分。我们据此获取了苯分子的 Fock 矩阵，用 `eigh` 将其对角化，并使用 Koopmans 定理估算电离能。

