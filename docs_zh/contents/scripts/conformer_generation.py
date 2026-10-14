#!/usr/bin/env python
# coding: utf-8

# # 构象搜索

# 本 Notebook 使用 ORCA Python 接口（OPI）通过力场搜索构象。由于所得能量排序可能不够准确，再使用 DFT 单点能计算重新排序。先演示 RDKit，再演示 GOAT。
# 
# 工作流程：
# 1. 导入所需 Python 依赖。
# 2. 定义工作目录。
# 3. 通过 SMILES 定义输入结构。
# 4. 使用 RDKit 生成构象。
# 5. 进行 DFT 单点计算，重新排序。
# 6. 可视化所得构象集合。
# 7. 通过 OPI 使用 GOAT 生成构象，并重新排序。

# ## 步骤 1：导入依赖
# 
# 首先导入以下用途所需的模块：
# - 访问 ORCA 输入与输出
# - RDKit
# - 数值计算与数据处理
# - 绘制结果
# 
# > **注意：** 我们还导入 `py3Dmol` 等可视化与绘图模块。可能需要将 `py3Dmol` 安装到 OPI 的 `venv` 中，例如激活 `.venv` 后运行 `uv pip install py3Dmol`。

# In[1]:


import re
import copy
import shutil
from pathlib import Path

# 使用 RDKit 生成构象
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem.rdchem import Mol

# 导入用于执行 ORCA 计算和读取输出的 OPI 模块
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import BasisSet, Dft, ForceField, Scf
from opi.input.simple_keywords.goat import Goat

# 用于结果绘图和分子可视化
import py3Dmol
import pandas as pd
import matplotlib.pyplot as plt


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


# ## 步骤 3：通过 SMILES 定义输入结构
# 
# 将输入结构定义为 SMILES 字符串，加载到 RDKit Mol 中，然后可视化二维结构。接着添加氢原子并识别手性中心。

# In[3]:


input_smiles='C1=CC(=CC=C1C[C@@H](C(=O)O)N)O'
# 将 SMILES 转换为 RDKit Mol，并标记原子以便可视化
mol = Chem.MolFromSmiles(input_smiles)
display(mol)
for atom in mol.GetAtoms():
    atom.SetProp('atomLabel', str(atom.GetIdx()))
display(mol)

# 添加氢原子并识别手性中心
mol = Chem.AddHs(mol)
chiral_centers = Chem.FindMolChiralCenters(mol)
if chiral_centers:
    print('Chiral centers identified:', chiral_centers)


# ## 步骤 4：使用 RDKit 生成构象
# 构象集合无需通过 OPI 生成，也可来自其他来源。这里使用 RDKit 生成的构象：

# In[4]:


def confGenerator_RDKit(smiles: str, e_thresh: float, confs: int, 
                        rms: float, method: str) -> tuple[list[tuple[float, int]], Mol, str]:
    """使用 RDKit 生成并优化构象"""

    # 生成三维构象
    conformer_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=confs,
        pruneRmsThresh=rms,
        randomSeed=1,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True
    )
    print('Number of raw conformers:', len(conformer_ids))

    # 优化每个构象
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=method)
    conformer_energies = []
    for conf_id in conformer_ids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
        ff.Minimize()
        energy_value = float(ff.CalcEnergy())
        conformer_energies.append((energy_value, conf_id))

    # 按能量排序，并保留低于能量阈值的构象
    conformer_energies.sort()
    min_energy = conformer_energies[0][0]

    filtered_mol = copy.deepcopy(mol)
    filtered_mol.RemoveAllConformers()
    selected_conf_ids = []
    filtered_rel_e = []

    data = []
    for idx, (energy_value, conf_id) in enumerate(conformer_energies, start=1):
        delta_e = energy_value - min_energy
        status = "below" if delta_e <= e_thresh else "above"
        data.append({
            "ConfID": conf_id,
            "Energy (kcal/mol)": f"{energy_value:.2f}",
            "ΔE (kcal/mol)": f"{delta_e:.2f}",
            "Status": status
        })
        if delta_e <= e_thresh:
            conf = mol.GetConformer(conf_id)
            filtered_mol.AddConformer(conf)
            selected_conf_ids.append(conf_id)
            filtered_rel_e.append((delta_e, conf_id))

    df = pd.DataFrame(data)
    print(f"\nNumber of conformers below {e_thresh} kcal/mol ({method}): {len(selected_conf_ids)}")
    display(df.style.hide(axis="index"))
    return filtered_rel_e, filtered_mol, method

# > 运行 RDKit 构象搜索
filtered_rel_e_rdkit, filtered_mol_rdkit, method_rdkit = confGenerator_RDKit(smiles=input_smiles, e_thresh=1.5, confs=60, rms=0.02, method='MMFF94s')


# ## 步骤 5：DFT 单点能计算
# 
# 可使用 OPI 执行 DFT 单点能计算，对生成的构象重新排序。下面的函数采用复合 DFT 方法 r²SCAN-3c，并直接使用 matplotlib 可视化排序结果。为降低本 Notebook 的计算成本，执行 DFT 计算的能量窗口设置得比正式计算所需的更小。

# In[5]:


def dft_calculations(smiles: str, e_thresh: float, functional: Dft, basis: BasisSet | None, 
                     rel_e: list[tuple[float, int]], mol: Mol, method: str, working_dir: Path = Path("RUN")) -> tuple[Mol, list[int]]:
    """使用 ORCA 执行 DFT 单点能计算"""
    smiles = re.sub(r'[^A-Za-z0-9]', '_', smiles)
    smiles_dir = working_dir / f"{smiles}_{method}"
    smiles_dir.mkdir(exist_ok=True)

    dft_energies = []
    comparison = []
    print("\nRunning DFT calculations with ORCA...")

    for delta_e, cid in rel_e:

        xyz_file = smiles_dir / f"conf_{cid}.xyz"
        xyz_block = Chem.MolToXYZBlock(mol, confId=cid)

        with open(xyz_file, 'w') as f:
            f.write(xyz_block)

        # > 设置结构
        xyz_file = Path(xyz_file)
        mol_name = xyz_file.stem
        structure = Structure.from_xyz(xyz_file)

        # > 定义计算文件目录
        calc_dir = xyz_file.parent / mol_name
        calc_dir.mkdir(exist_ok=True)

        # > 设置计算器
        calc = Calculator(basename=mol_name, working_dir=calc_dir)
        calc.structure = structure
        # > 中性结构
        calc.structure.charge = 0
        # > 自旋多重度为 1（闭壳层）
        calc.structure.multiplicity = 1
        # r²SCAN-3c 等方法有各自预定义的基组
        if basis:
            calc.input.add_simple_keywords(functional, basis, Scf.NOAUTOSTART)
        else:
            calc.input.add_simple_keywords(functional, Scf.NOAUTOSTART)

        # > 使用 4 个核心运行计算
        calc.input.ncores = 4

        # > 写出 ORCA 输入文件
        calc.write_input()

        # > 运行计算
        calc.run()

        # > 获取并检查输出
        output = calc.get_output()
        # > 检查 ORCA 是否正常终止
        status = output.terminated_normally()
        if not status:
            # > ORCA 未正常终止
            raise RuntimeError(f"ORCA did not terminate normally, see output file: {output.get_outfile()}")
        else:
            # > ORCA 已正常终止，可以解析输出
            output.parse()

        # > 现在检查 SCF 是否收敛
        if not output.results_properties.geometries[0].single_point_data.converged:
            raise RuntimeError("SCF DID NOT CONVERGE")

        # > 获取 DFT 计算的能量
        e_dft = output.results_properties.geometries[0].single_point_data.finalenergy

        dft_energies.append((e_dft, cid))
        comparison.append({
            'Conformer': cid,
            f'{method} (kcal/mol)': delta_e,                
            'DFT (kcal/mol)': 0
        })
        print(f"Conformer {cid}: DFT Energy = {e_dft:.6f} Eh")

    # 计算相对 DFT 能量
    dft_energies.sort()
    min_dft = dft_energies[0][0]

    for entry in comparison:
        cid = entry['Conformer']
        for e_dft, dft_id in dft_energies:
            if cid == dft_id:
                entry['DFT (kcal/mol)'] = (e_dft - min_dft) * unit_conversion
                break

    df_dft = pd.DataFrame(comparison).sort_values(f'{method} (kcal/mol)')

    # 绘制能量差
    x_vals = range(len(df_dft))
    conformer_ids = df_dft['Conformer']

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, df_dft[f'{method} (kcal/mol)'], 'bo', label=f'{method}')
    plt.plot(x_vals, df_dft['DFT (kcal/mol)'], 'ro', label='DFT')
    plt.xticks(x_vals, conformer_ids)
    plt.xlabel('Conformer ID')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title(f'{method} vs DFT Relative Energies')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # 比较按力场能量和 DFT 能量得到的构象排序
    print("\nEnergy Comparison Table:")
    method_rank = df_dft.sort_values(f'{method} (kcal/mol)')['Conformer'].tolist()
    dft_rank = df_dft.sort_values('DFT (kcal/mol)')['Conformer'].tolist()

    comparison_table = []
    for _, row in df_dft.iterrows():
        cid = row['Conformer']
        method_r = method_rank.index(cid) + 1
        dft_r = dft_rank.index(cid) + 1
        delta = method_r - dft_r
        mark = f"{'↑' if delta > 0 else '↓' if delta < 0 else '='}{abs(delta)}" if delta else "-"

        comparison_table.append({
            "ConfID": int(cid),
            f"{method} (kcal/mol)": f"{row[f'{method} (kcal/mol)']:.2f}",
            "DFT (kcal/mol)": f"{row['DFT (kcal/mol)']:.2f}",
            "Rank Change": mark
        })

    df_rank = pd.DataFrame(comparison_table)
    display(df_rank.style.hide(axis="index"))

    # 根据 DFT 能量阈值进行最终筛选
    filtered_mol_dft = copy.deepcopy(mol)
    filtered_mol_dft.RemoveAllConformers()
    final_ids = []

    data_dft = []
    for idx, (e_dft, conf_id) in enumerate(dft_energies, start=1):
        delta_dft = (e_dft - min_dft) * 627.509
        status = "below" if delta_dft <= e_thresh else "above"
        data_dft.append({
            "ConfID": conf_id,
            "DFT Energy (Eh)": f"{e_dft:.6f}",
            "ΔE (kcal/mol)": f"{delta_dft:.2f}",
            "Status": status
        })
        if delta_dft <= e_thresh:
            filtered_mol_dft.AddConformer(mol.GetConformer(conf_id))
            final_ids.append(conf_id)

    df_dft_filtered = pd.DataFrame(data_dft)
    print(f"\nNumber of conformers below {e_thresh} kcal/mol (DFT): {len(final_ids)}")
    display(df_dft_filtered.style.hide(axis="index"))
    return filtered_mol_dft, final_ids

# > 通过 OPI 使用 r²SCAN-3c 执行 DFT 单点能计算
filtered_mol_dft_rdkit, final_ids_rdkit = dft_calculations(smiles=input_smiles, e_thresh=1.0, functional= Dft.R2SCAN_3C, basis=None, rel_e=filtered_rel_e_rdkit, mol=filtered_mol_rdkit, method=method_rdkit, working_dir=working_dir)


# ## 步骤 6：使用 py3Dmol 可视化选定构象
# 
# 最终构象集合可通过 py3Dmol 可视化：

# In[6]:


def visualization(align: list[int], highlight: list[int] | None, mol: Mol, ids: list[int]) -> None:
    """使用 py3Dmol 可视化选定构象，并可选择高亮显示"""
    if not isinstance(highlight, list):
        highlight = [highlight]
    # 使用 py3Dmol 可视化最终构象
    view = py3Dmol.view(width=500, height=500)
    mol = Chem.RemoveHs(mol)

    rdMolAlign.AlignMolConformers(mol, atomIds=align)

    for idx, conf_id in enumerate(ids, start=1):
        view.addModel(Chem.MolToMolBlock(mol, confId=conf_id))
        model = view.getModel()
        if not highlight or idx in highlight:
            model.setStyle({'stick': {'opacity': 1, 'radius': 0.2}})
        else:
            model.setStyle({'stick': {'opacity': 0.8, 'radius': 0.1}})

    view.setBackgroundColor('white')
    view.zoomTo()
    view.show()

visualization(align=[0,1,2,3], highlight=1, mol=filtered_mol_dft_rdkit, ids=final_ids_rdkit)


# ## 步骤 7：使用 GOAT 更充分地生成构象
# [GOAT](https://doi.org/10.1002/anie.202500393) 能更充分地生成构象集合，但计算成本也更高。这里将 GOAT 与 GFN-FF 配合使用，原则上也可使用 ORCA 中的任何方法。之后使用 r²SCAN-3c 重新排序，并按前述方式可视化结果。

# In[7]:


def confGenerator_GOAT(smiles: str, e_thresh: float, method: ForceField, working_dir: Path = Path("RUN"))-> tuple[list[tuple[float, int]], Mol, str]:
    """使用 ORCA GOAT 生成并优化构象"""

    # 生成初始三维结构
    smiles = re.sub(r'[^A-Za-z0-9]', '_', smiles)
    smiles_dir = working_dir / Path(f"{smiles}_{method}")
    smiles_dir.mkdir(exist_ok=True)

    xyz_file = smiles_dir / f"{smiles}.xyz"
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    xyz_block = Chem.MolToXYZBlock(mol)

    with open(xyz_file, 'w') as f:
        f.write(xyz_block)

    # > 设置 GOAT 计算

    # > 准备结构
    xyz_file = Path(xyz_file)
    mol_name = xyz_file.stem
    structure = Structure.from_xyz(xyz_file)

    # > 设置计算器
    calc_dir = xyz_file.parent / mol_name
    calc_dir.mkdir(exist_ok=True)
    calc = Calculator(basename=mol_name, working_dir=calc_dir)
    calc.structure = structure
    # > 中性分子
    calc.structure.charge = 0
    # > 自旋多重度为 1（闭壳层）
    calc.structure.multiplicity = 1

    calc.input.add_simple_keywords(
        method,
        Goat.GOAT
    )

    # > 计算使用 8 个核心
    calc.input.ncores = 8

    # > 写出 ORCA 输入文件
    calc.write_input()

    # > 运行计算
    calc.run()

    # 构象搜索结果位于 basename.finalensemble.xyz 文件中
    result_file = calc_dir / f"{mol_name}.finalensemble.xyz"

    with open(result_file, 'r') as f:
        lines = f.readlines()

    i = 0
    conformer_energies = []
    mol.RemoveAllConformers()
    conf_counter = 0

    while i < len(lines):
        atom_count = int(lines[i].strip())
        energy_line = lines[i + 1].strip()
        atoms_block = lines[i + 2:i + 2 + atom_count]

        match = re.match(r"([-+]?[0-9]*\.?[0-9]+)", energy_line)
        energy = float(match.group(1)) * unit_conversion

        conf = Chem.Conformer(atom_count)
        for idx, line in enumerate(atoms_block):
            _, x, y, z = line.split()
            conf.SetAtomPosition(idx, (float(x), float(y), float(z)))

        conf_id = mol.AddConformer(conf, assignId=True)
        conformer_energies.append((energy, conf_id))

        i += 2 + atom_count
        conf_counter += 1
    print('Number of raw conformers:', len(conformer_energies))

    # 按能量排序，并保留低于能量阈值的构象
    conformer_energies.sort()
    min_energy = conformer_energies[0][0]

    filtered_mol = copy.deepcopy(mol)
    filtered_mol.RemoveAllConformers()
    selected_conf_ids = []
    filtered_rel_e = []

    data = []
    for idx, (energy_value, conf_id) in enumerate(conformer_energies, start=1):
        delta_e = energy_value - min_energy
        status = "below" if delta_e <= e_thresh else "above"
        data.append({
            "ConfID": conf_id,
            "Energy (kcal/mol)": f"{energy_value:.2f}",
            "ΔE (kcal/mol)": f"{delta_e:.2f}",
            "Status": status
        })
        if delta_e <= e_thresh:
            conf = mol.GetConformer(conf_id)
            filtered_mol.AddConformer(conf)
            selected_conf_ids.append(conf_id)
            filtered_rel_e.append((delta_e, conf_id))

    df = pd.DataFrame(data)
    print(f"\nNumber of conformers below {e_thresh} kcal/mol ({method}): {len(selected_conf_ids)}")
    display(df.style.hide(axis="index"))
    return filtered_rel_e, filtered_mol, method

# > 使用 GOAT 和 GFN-FF 生成构象，能量窗口设为 1.5 kcal/mol
filtered_rel_e_goat, filtered_mol_goat, method_goat = confGenerator_GOAT(smiles=input_smiles, e_thresh=1.5, method=ForceField.GFN_FF, working_dir=working_dir)
# > 执行 DFT 单点能计算
filtered_mol_dft_goat, final_ids_goat = dft_calculations(smiles=input_smiles, e_thresh=1.0, functional= Dft.R2SCAN_3C, basis= None, rel_e=filtered_rel_e_goat, mol=filtered_mol_goat, method=method_goat, working_dir=working_dir)
# > 可视化结果
visualization(align=[0,1,2,3], highlight=1, mol=filtered_mol_dft_goat, ids=final_ids_goat)


# ## 小结
# 
# 本 Notebook 演示了构象处理。使用 RDKit 和 GOAT 生成构象，再通过 DFT 重新排序。这两种操作均可直接通过 OPI 完成，使常见构象工作流程更便捷，并在 Notebook 内可视化结果。

