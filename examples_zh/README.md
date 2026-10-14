# ORCA Python 接口（OPI）示例

本目录是仓库 `examples/` 的简体中文副本，包含 58 个示例的脚本与输入文件。文件名及目录层次与原版一致；教程说明、代码注释和模块文档字符串已翻译，计算代码、运行时字符串及配套数据保持原样。

原版示例有以下**两个用途**：

- **学习资源**：向新用户演示如何使用 OPI 设置和运行典型计算。
- **集成测试**：测试套件（见 [`tests/examples`](../tests/examples)）会执行原版 `examples/` 下的示例，以检查 OPI 基本功能与 ORCA 接口是否正常。中文副本供阅读和手动运行。

## 用户须知

手动运行示例时，在对应示例文件夹中执行：

```bash
python3 job.py
```

运行所需的 OPI、ORCA 及可选依赖与原版相同。运行时以当前示例文件夹为工作位置，使脚本能够找到相对路径中的结构文件。

## 贡献者须知

- 修改原版示例后，应确保相应测试仍然通过，并同步中文副本。
- 添加新示例时，请将其加入下方列表，并尽可能添加测试。
- 可在仓库根目录运行以下完整性检查，验证中文脚本的语法、计算代码及配套文件：

```bash
python3 docs_zh/_tools/check_examples_translation.py
```

检查允许修改注释和文档字符串，并比较其余 Python 语法树；同时核对公式、行内代码、引用链接及配套文件的二进制内容。检查本身不执行 ORCA 计算。

## 示例列表

| 目录 | 内容 |
| --- | --- |
| [exmp001_scf](exmp001_scf/) | 执行 HF/def2-SVP 单点计算 |
| [exmp002_scf_ccsdt](exmp002_scf_ccsdt/) | 执行 CCSD(T)/def2-SVP 单点计算 |
| [exmp003_opt](exmp003_opt/) | 执行 wB97X-3c 几何优化 |
| [exmp004_freq](exmp004_freq/) | 执行 TPSS/def2-SVP 频率计算 |
| [exmp005_dft](exmp005_dft/) | 使用自定义 D3 参数执行 B3LYP-D3/def2-SVP 单点计算 |
| [exmp006_mp2](exmp006_mp2/) | 执行 MP2/def2-TZVP 单点计算 |
| [exmp007_cpcm](exmp007_cpcm/) | 采用水的 CPCM 隐式溶剂化模型，执行 BP86/def2-SVP 单点计算 |
| [exmp008_relativ_corr](exmp008_relativ_corr/) | 采用 ZORA 近似处理标量相对论效应，执行 BP86/SARC-ZORCA-SVP 计算 |
| [exmp009_rama](exmp009_rama/) | 使用 PBE0/def2-SVP 计算拉曼光谱 |
| [exmp010_uvvis](exmp010_uvvis/) | 使用 B3LYP/def2-TZVP + CPCM(Hexane) 执行 TD-DFT 计算 |
| [exmp011_epr](exmp011_epr/) | 使用 B3LYP/EPR-II 计算 EPR 性质 |
| [exmp012_nmr](exmp012_nmr/) | 使用 RI-MP2/pcSseg-2 计算 NMR 性质 |
| [exmp013_bs](exmp013_bs/) | 使用 B3LYP/SVP（B3LYP/G）执行破缺对称性计算 |
| [exmp014_led](exmp014_led/) | 使用 DLPNO-CCSD(T)/cc-pVDZ 执行 LED 分解，包括 ADLD 和 ADEX |
| [exmp015_pop_analysis](exmp015_pop_analysis/) | 获取 DSD-PBEP86/def2-SVP 的多种布居分析结果 |
| [exmp016_autoci](exmp016_autoci/) | 使用 AUTO-CI 执行 CISD/def2-SVP 单点计算 |
| [exmp017_roci](exmp017_roci/) | 执行 ROCIS 计算 |
| [exmp018_cipsi](exmp018_cipsi/) | 执行 ICE-CI 计算 |
| [exmp019_engrad](exmp019_engrad/) | 使用 r²SCAN-3c 计算能量和梯度 |
| [exmp020_smd](exmp020_smd/) | 执行 r²SCAN-3c+SMD(Water) 计算 |
| [exmp021_basis](exmp021_basis/) | 为氧添加弥散函数，执行 BP86/def2-SVP 能量计算 |
| [exmp022_scf_block](exmp022_scf_block/) | 执行 PBE0/def2-SVP 能量计算，并旋转初始 SCF 猜测轨道 |
| [exmp023_thermo](exmp023_thermo/) | 从频率计算中获取统计热力学校正 |
| [exmp024_blocks](exmp024_blocks/) | 演示输入块的使用 |
| [exmp025_goat](exmp025_goat/) | 执行 GOAT 构象搜索 |
| [exmp026_scan](exmp026_scan/) | 扫描键长、键角和二面角 |
| [exmp027_constraints](exmp027_constraints/) | 执行带约束的几何优化 |
| [exmp028_nevpt2](exmp028_nevpt2/) | 执行 NEVPT2 计算 |
| [exmp029_oomp2](exmp029_oomp2/) | 执行 OO-MP2 计算 |
| [exmp030_eom](exmp030_eom/) | 执行 EOM-CCSD 计算 |
| [exmp031_dlpno_ccsdt](exmp031_dlpno_ccsdt/) | 执行 DLPNO-CCSD(T) 计算 |
| [exmp032_fragbasis](exmp032_fragbasis/) | 为不同片段指定不同基组 |
| [exmp033_mo_getters](exmp033_mo_getters/) | 访问分子轨道数据 |
| [exmp034_strucfile](exmp034_strucfile/) | 通过外部文件向 ORCA 传递结构 |
| [exmp035_moplot](exmp035_moplot/) | 根据密度或分子轨道生成用于绘图的 Cube 文件 |
| [exmp036_solvator](exmp036_solvator/) | 运行 ORCA 自动溶剂化工作流 solvator |
| [exmp037_s2](exmp037_s2/) | 从非限制性计算中获取 S² 期望值 |
| [exmp038_integrals](exmp038_integrals/) | 从计算中获取积分 |
| [exmp039_neb](exmp039_neb/) | 使用 NEB 搜索过渡态 |
| [exmp040_xzyfraglib](exmp040_xzyfraglib/) | 使用外部片段库 frag_lib.xyz 分配片段 |
| [exmp041_graph](exmp041_graph/) | 输出可用计算结果的概览 |
| [exmp042_element](exmp042_element/) | 访问结构中各元素的原子序数 |
| [exmp043_dummy_atom](exmp043_dummy_atom/) | 执行包含虚原子的计算 |
| [exmp044_from_ase](exmp044_from_ase/) | 从 ASE Atoms 对象获取结构 |
| [exmp045_existing_calc](exmp045_existing_calc/) | 获取已有计算的结果 |
| [exmp046_server](exmp046_server/) | 演示如何通过 OpiServer 使用 ExtOpt 方法 |
| [exmp047_trj](exmp047_trj/) | 对多结构 XYZ 轨迹文件中的所有结构执行计算 |
| [exmp048_loc_gbw](exmp048_loc_gbw/) | 读取局域化分子轨道系数，并将积分变换到 LMO 基下 |
| [exmp049_arbitrary_block_variable](exmp049_arbitrary_block_variable/) | 在 SCF 块中使用任意块变量 |
| [exmp050_docker](exmp050_docker/) | 使用 DOCKER 将水与水进行对接 |
| [exmp051_libxc](exmp051_libxc/) | 通过 method 块修改 DFT LibXC 参数 |
| [exmp052_densities](exmp052_densities/) | 执行 UHF 计算，并获取 AO 基下的密度和自旋密度 |
| [exmp053_opencosmors](exmp053_opencosmors/) | 运行 OpenCOSMO-RS 任务 |
| [exmp054_casscf_fcidump](exmp054_casscf_fcidump/) | 执行 CASSCF 计算，并将活性空间积分导出为 FCIDUMP 文件 |
| [exmp055_gfnff_fallback](exmp055_gfnff_fallback/) | 执行 GFN-FF 优化，并提取能量、梯度和结构 |
| [exmp056_optts_freq](exmp056_optts_freq/) | 执行 GFN2-xTB 过渡态优化和频率计算，并提取虚频 |
| [exmp057_to_ase](exmp057_to_ase/) | 将优化后的 Structure 转换为 ASE Atoms 对象 |
| [exmp058_write_gbw](exmp058_write_gbw/) | 根据已解析的输出数据写出 GBW 文件，并将其用作 SCF 初猜 |
