# OPI 中文手册源码

本目录是 [faccts/opi](https://github.com/faccts/opi) 中 `docs/` 的简体中文译本，依据上游提交 `ab7cad88803cde630ddfee07dfcc269d2d5920f4`。

包括 6 个 Markdown 文档、19 个 Jupyter 教程、14 个中文注释脚本以及自动生成的 API 参考。文件名、目录层次、图片、样式和模板沿用英文版；Notebook 的代码单元、计算输出与元数据保持原样。脚本中的教程说明、注释和函数文档字符串已翻译，计算代码与运行时字符串保持原样。说明文字中的公式、命令、参数名、单位和引用链接予以保留。中文首页将 API 与入门教程入口指向本目录生成的页面。

## 目录

| 路径 | 内容 |
| --- | --- |
| `index.md` | 中文首页及导航 |
| `contents/*.md` | 安装、版本、引用方式、贡献指南和开发指南 |
| `contents/notebooks/*.ipynb` | 19 个完整中文教程 |
| `contents/scripts/` | 14 个示例脚本，含中文教程说明、注释和函数文档字符串 |
| `img/`、`_static/`、`_templates/` | 原版图片、静态资源和模板 |
| `_translations/api_lines.json` | API 说明的中英对照表 |
| `_translations/source.json` | 上游来源与翻译范围 |
| `_tools/check_translation.py` | 源码、代码单元和资源完整性检查 |
| `conf.py`、`Makefile` | 中文 Sphinx 配置与构建命令 |

`contents/api/` 在构建时由 `src/opi/` 生成。Sphinx 在解析 NumPy 风格说明之前应用 API 译文表，并规范上游少量列表、参数段及代码块的 reStructuredText 格式，保留参数名称、类型、类名和可执行示例。维护译文时修改对照表即可。上游增加或修改说明后，应更新相应译文键值。

## 构建 HTML

要求 Python 3.11 或更高版本。在仓库根目录运行：

```bash
cd docs_zh
make html
```

Makefile 会创建 `.venv` 并使用仓库的 `pyproject.toml`、`uv.lock` 安装文档依赖。HTML 入口为 `_build/html/index.html`。Notebook 使用已有输出构建，构建过程不会执行其中的 ORCA 计算。

已安装文档依赖时，也可以从仓库根目录直接运行：

```bash
python -m sphinx -b html docs_zh docs_zh/_build/html
```

## 其他格式

同一套源码支持 Sphinx 的其他输出格式：

```bash
cd docs_zh
make singlehtml
make epub
make latex
```

对应输出位于 `_build/singlehtml/`、`_build/epub/`、`_build/latex/`。`latex` 生成 LaTeX 源码；编译中文 PDF 还需要安装支持中文的 TeX 环境与字体，例如使用 XeLaTeX。

## 术语约定

| 英文 | 中文 |
| --- | --- |
| population analysis | 布居分析 |
| multiplicity | 多重度 |
| local energy decomposition (LED) | 局域能量分解 |
| London dispersion | 伦敦色散 |
| conical intersection | 圆锥交叉 |
| broken symmetry | 破缺对称性 |
| ghost atom | 幽灵原子 |
| core Hamiltonian | 芯哈密顿量 |
| Hessian | Hessian 矩阵 |
| free energy | 自由能 |

软件名称、方法缩写和 API 标识符保留原文。示例脚本的说明文字已翻译；运行时字符串（含控制台输出、异常信息和图中文字）、文件名及数据键名保持原样，便于对照实际运行结果。Notebook 代码单元中的注释继续保留原文。原文中明确的术语笔误在译文中按上下文翻译，必要时注明原文。

## 完整性检查与上游更新

从仓库根目录运行：

```bash
python docs_zh/_tools/check_translation.py
```

检查 Notebook 的代码、输出、元数据、公式和代码块，以及原版资源文件。对中文脚本进行语法编译检查，并比较去除文档字符串后的 Python 抽象语法树，同时核对注释中的公式、代码标识符及引用链接，以确保翻译未改变计算代码。新增上游文档后，按 `docs/` 的相对路径在 `docs_zh/` 中增加译文，并同步 API 译文表。脚本与 Notebook 的上游内容可能不同，翻译脚本时应以脚本自身的步骤编号、路径和说明为准。

中文构建沿用上游 AutoAPI 设置。上游存在重复对象描述、交叉引用歧义及社区示例中缺失的 `gwqm9.csv` 链接等构建警告；示例数据链接与已有计算输出予以保留。译本遵循仓库现有的 `LICENSE` 和 `CITATION.cff`，作者及文献引用保持不变。

