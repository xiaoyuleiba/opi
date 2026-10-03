# 如何贡献

ORCA Python 接口 OPI 被设计为开源社区项目，欢迎外部开发者贡献。为保持项目的高质量，请所有贡献者遵守相应标准，并提前阅读这里提供的信息。

## GitHub

OPI 的开发基于 [GitHub](https://github.com/faccts/opi)。
如需贡献，建议 fork 项目，并通过拉取请求提交修改与改进。
请使用拉取请求模板，在其中注明该请求关闭或影响的 issue，简要说明修改，并记录应加入发行说明的内容。
项目的发行说明遵循 [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) 约定；条目应以大写字母开头，并以拉取请求编号结尾（例如 `(#42)`）。
问题和错误通过 [GitHub issues](https://github.com/faccts/opi/issues) 跟踪。

## 许可证与 CLA

### 开源许可证

本开源项目按照 `GPL-3.0` 开源许可证公开发布。该许可证适用于所有公开发布的代码，允许任何人在遵守其条款的前提下自由使用、修改和分发本项目。

### 专有许可证

程序及全部贡献也可能按照商业许可证纳入我们的专有软件产品。由此，我们可以：

- 将开源组件和闭源组件组合为同一产品；
- 为有特定商业需求的客户提供其他许可条款；
- 确保所有公开部分符合开源要求，同时简化私有或嵌入式分发的许可义务。

### 贡献者许可协议（CLA）

为维持这一许可模式，所有贡献者必须签署贡献者许可协议（CLA）。该协议由行业通用的 Apache CLA 略加修改而成。签署后，您：

- 保留对自己贡献的所有权；
- 授予我们非独占许可，允许我们按照开源和专有许可条款使用、再许可、重新许可和分发您的贡献。

### 我们采用两部分 CLA 体系：

- [个人贡献者许可协议（ICLA）](https://github.com/faccts/opi/blob/main/CLA.md)，用于个人贡献；
- 企业贡献者许可协议（CCLA），用于代表雇主提交的贡献，可向 info@faccts.de 索取。

(sec:howto_code)=
## 代码与开发指南

代码及开发规范请参阅[开发指南](dev_guide.md)。

(sec:howto_examples)=
### 示例 Notebook 与脚本

社区示例可以采用 Jupyter Notebook 或 Python 脚本形式提交。无论采用何种形式，都应提供包含相应 Markdown 说明的 Jupyter Notebook 文档。

- Jupyter Notebook 必须存放在 [docs/notebooks](https://github.com/faccts/opi/tree/main/docs/contents/notebooks) 目录中，扩展名为 `.ipynb`。
- Python 脚本应存放在 [docs/scripts](https://github.com/faccts/opi/tree/main/docs/contents/scripts) 目录中。

  :::{tip}
  可以使用 `jupyter` 和 `jupytext` 包在 Jupyter Notebook 与 Python 脚本之间互相转换：

  ```
  jupyter nbconvert --to script example.ipynb
  ```

  ```
  jupytext example.py --to notebook
  ```
  :::

- 文件名请采用 snake_case，用下划线分隔单词（例如 `this_is_my_file.py`）。
- 请确保两种格式的脚本都有充分的[文档说明](sec:howto_docs)，便于他人理解。
- 请只提交功能完整、可正常运行的代码。
- 欢迎社区长期维护所贡献的示例 Notebook。FACCTs 团队只能维护官方教程示例。

  :::{note}
  社区 Notebook 及其中链接的外部软件由第三方贡献，未经 FACCTs 审查、认可或控制。
  对使用这些材料所产生的内容、功能或后果，我们不承担责任。
  使用前，用户应独立核实外部软件或资源的合法性、安全性及适用性。
  :::

(sec:howto_docs)=
## 文档页面

请尽可能完整、清楚地说明每个贡献的示例。文档页面主要采用 [Sphinx](https://www.sphinx-doc.org/en/master/) 与 [MyST Markdown](https://mystmd.org/) 组合。MyST Markdown 是一种易学、宽容的 Markdown 方言，便于编写文本式文档。

:::{tip}
如果更熟悉 ReStructuredText（ReST），可以使用[转换工具](https://rst-to-myst.readthedocs.io/en/latest/index.html) 将 ReST 转换为 MyST。
:::

- 如需将 Jupyter Notebook 直接纳入文本式文档，只需把它加入 [index.md](https://github.com/faccts/opi/blob/main/docs/index.md) 文件（例如 `docs/contents/notebooks/example.ipynb`）。`myst_nb` Sphinx 扩展会将 Notebook 及其中的 Markdown 文档渲染为独立文档页面。
- 如需说明某个 Python 脚本或代码，可以创建 Markdown 页面（例如 `example.md`），并通过以下方式纳入完整脚本：

  ````
  ```{literalinclude} ../docs/scripts/example.py
  :language: python
  ```
  ````

- 文本文档请尽量采用 Markdown 语法。Markdown 比 ReST 更宽容，也更便于新贡献者学习。
