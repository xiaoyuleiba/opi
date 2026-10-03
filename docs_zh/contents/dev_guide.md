# 开发指南

本指南面向参与本项目的开发者。

## Python

OPI 完全使用 Python 编写，因此只应提交能正常运行的 Python 代码。

## 依赖管理：uv

OPI 项目使用 [uv](https://github.com/astral-sh/uv) 管理依赖和安装。uv 是使用 Rust 编写的高速 Python 包及项目管理器。
uv 有多种安装方式，我们倾向于使用 [pipx](https://github.com/pypa/pipx)：

```
pipx install uv
```

## 构建后端：Hatchling

项目采用 [Hatchling](https://pypi.org/project/hatchling/) 作为构建后端。在项目根目录中执行：

```
uv build [--wheel] [--sdist]
```

构建产物位于 `dist/`。

## 开发者安装

在根目录中运行以下命令，以开发模式安装此包：

```
uv sync
```

这会创建 Python 版本匹配的虚拟环境（venv），并在其中安装该包：

- `--only-dev` 还会安装全部开发依赖。
- `--group <groupname>` 可包含指定依赖组中的依赖。
- `--all-groups` 将所有组的全部依赖加入虚拟环境。

虚拟环境通常位于项目根目录的 `.venv/`。
**安装前应删除已有的 `.venv` 目录！**

## 依赖管理

uv 采用版本锁定，把每个必需包及其间接依赖的准确版本存入中央锁文件 `uv.lock`，类似 Rust 的 `cargo.lock`。

**不应手动编辑 `uv.lock`！**

直接依赖是项目直接依赖的外部包，在 `pyproject.toml` 中指定。
使用 uv 的 `add` 命令添加依赖，例如添加 Numpy：

```
uv add numpy
```

Numpy 会加入 `pyproject.toml` 和锁文件；锁文件也会包含 Numpy 的全部间接依赖。

可以按照 [PEP 508](https://peps.python.org/pep-0508/#grammar) 设置版本约束。通常无需这样做，但新版本存在尚未修复的错误，或两个主版本之间的 API 变化破坏了项目时，可能需要限制版本。
例如，使用最新的 1.x 版本而不升级到 2.0 主版本：

```
uv add 'numpy>=1,<2'
```

使用版本约束时，通常应将包名放在引号中，避免 shell 转义问题。

删除依赖：

```
uv remove numpy
```

### 依赖分组

多数依赖管理器允许将依赖分组，以实现模块化，或区分开发要求和运行要求。
未指定组时，使用默认组（主组），其中包含运行此包所需的全部依赖。
使用 `--group <groupname>` 将依赖加入某个组：

```
uv add --group math 'numpy>=1,<2'
```

普通依赖和开发依赖（`--dev`）均适用：

```
uv add --dev 'numpy>=1,<2'
```

依赖可以同时加入多个组。可分别为每个组中的依赖设置版本约束，但约束之间不能互相排斥。最终只会安装一个满足所有约束的版本。

## Nox

[Nox](https://nox.thea.codes/en/stable/) 是自动化测试框架，通过中央配置文件 `noxfile.py` 配置。
**只有明确理解其作用时才应修改此文件！**

我们使用 Nox 执行测试，并驱动用于确保代码质量和安全性的其他工具。后文会逐一介绍。
在项目根目录中执行：

```
nox <additional nox arguments>
```

不提供其他参数时，Nox 执行 `noxfile.py` 中定义的全部默认 session。在 Nox 中，session 是执行代码并返回状态的预定义 Python 函数。
每个 session 均在独立的 Python 环境中执行，环境由 Nox 按配置建立和管理，从而确保不同环境中的结果一致。

列出可用 session：

```
nox -l
```

运行指定 session：

```
nox -l <SESSION NAME> ...
```

我们将 Nox 配置为任一 session 失败后立即停止，后续 session 均跳过。添加 `--no-stop-on-first-error` 可关闭这一行为。
**如果需要检查某个具体方面，并且必须运行全部期望的 session，就应关闭这一行为！**

Nox 也允许给 session 添加标签：

```
@nox.session(tags=["static_check"])
def type_check(session):
    session.run_always("uv", "sync", "--group", "type_check")
    session.run("mypy")
```

这样可以把多个 session 分组。按标签运行：

```
nox -t <TAG> ...
```

所有 session 均依次独立执行，互不影响。

## 我们的 Nox session

下面简要介绍与 Nox 配合使用的工具。

### 静态类型检查：mypy

项目使用 [mypy](https://mypy-lang.org/) 确保类型安全。Mypy 基于 Python 类型注解执行静态类型检查。“静态”表示不执行待检查的代码。

```
nox -s type_check
```

输出示例：

```
src/opi/output/core.py:1: error: Cannot find implementation or library stub for module named "calc_status"  [import-not-found]
    import calc_status
    ^
src/opi/output/core.py:1: note: See https://mypy.readthedocs.io/en/stable/running_mypy.html#missing-imports
src/opi/output/core.py:2: error: Cannot find implementation or library stub for module named "geometry"  [import-not-found]
[...]
src/opi/output/geometry.py:1: error: Cannot find implementation or library stub for module named "coord"  [import-not-found]
    import coord
    ^
Found 30 errors in 3 files (checked 15 source files)
```

在暂存用于 PR 的修改之前，应解决所有 Mypy 错误。

### 代码检查与格式化：Ruff

[Ruff](https://docs.astral.sh/ruff/) 是目前最快的 Python 代码检查和格式化工具。本项目将其用于：

- 代码检查：

```
nox -s lint
```

- 导入排序，替代 [isort](https://pypi.org/project/isort/)：

```
nox -s imports
```

- 自动执行严格的代码格式化，替代 [Black](https://pypi.org/project/black/)：

```
nox -s format_code
```

### 拼写检查：Codespell

我们使用 [Codespell](https://github.com/codespell-project/codespell) 检查语法与拼写。`codespell` 按照美式英语检查。

```
nox -s spell_check
```

输出示例：

```
./examples/exmp003_opt/job.out:1004: contruction ==> construction
./examples/exmp003_opt/job.out:1480: contruction ==> construction
./examples/exmp003_opt/job.out:1954: contruction ==> construction
./examples/exmp003_opt/job.out:2437: contruction ==> construction
./src/opi/input/models/calculator.py:264: inout ==> input, in out
./build/lib/opi/input/models/calculator.py:264: inout ==> input, in out
6
```

将需要排除检查的单词或短语加入 `.codespellignore`；每个单词或短语必须以换行结尾。

### API 文档：AutoAPI

项目使用 [AutoAPI](https://sphinx-autoapi.readthedocs.io/en/latest/#) 从文档字符串自动生成 API 文档。
文档字符串遵循 [numpydoc 风格约定](https://numpydoc.readthedocs.io/en/latest/format.html)，由 Sphinx 的 Napoleon 扩展处理，以正确格式化并整合到文档中。

### 测试：Pytest

[Pytest](https://docs.pytest.org/en/stable/) 是 Python 单元测试领域事实上的行业标准。所有测试及相关文件位于 `tests/`。
`tests/conftest.py` 是 Pytest 配置文件。

我们使用以下插件：

- [pytest-randomly](https://pypi.org/project/pytest-randomly/)，使单元测试始终按随机顺序执行，避免执行顺序导致偏差。
- [Hypothesis](https://hypothesis.readthedocs.io/en/latest)，根据预定义约束生成随机输入，避免测试始终使用相同输入而导致偏差。

[Fixture](https://docs.pytest.org/en/stable/explanation/fixtures.html) 位于 `tests/fixtures`，应按用途分组到不同模块中。

执行 pytest：

```
nox -s tests [--] [<pytest arguments>] 
```

`--` 后的参数直接传给 Pytest。

## 风格指南

我们严格遵循 [Black](https://black.readthedocs.io/en/stable/) 的规则。出于性能考虑，使用可直接替代 Black 的 Ruff。
Ruff 与 Black 并非完全兼容，详见 [Ruff 格式化差异](https://docs.astral.sh/ruff/formatter/#intentional-deviations)。

### 代码约定

变量、函数及其他标识符的命名遵循 PEP 8，以确保代码的一致性和可读性。
每级缩进统一为 4 个空格，禁止使用制表符。
此外，最大行长扩展为 100 个字符，以提高清晰度并更好地利用现代屏幕宽度。

### 文件名

请使用 snake_case 文件名，以下划线分隔单词（例如 `this_is_my_file.py`）。
