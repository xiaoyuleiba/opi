#!/usr/bin/env python3
"""核对 examples_zh 的计算代码与配套文件，不执行量子化学计算。"""

import ast
from collections import Counter
import io
from pathlib import Path
import re
import subprocess
import tokenize


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'examples'
TARGET = ROOT / 'examples_zh'


class _WithoutDocstrings(ast.NodeTransformer):
    def _strip(self, node):
        self.generic_visit(node)
        if ast.get_docstring(node) is not None:
            node.body.pop(0)
        return node

    visit_Module = _strip
    visit_ClassDef = _strip
    visit_FunctionDef = _strip
    visit_AsyncFunctionDef = _strip


def documentation(text):
    """提取可翻译的注释和文档字符串。"""
    parts = [t.string for t in tokenize.generate_tokens(io.StringIO(text).readline)
             if t.type == tokenize.COMMENT]
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node)
            if doc is not None:
                parts.append(doc)
    return '\n'.join(parts)


def protected_parts(text):
    patterns = {
        '公式': r'(?<!\\)\$\$.*?(?<!\\)\$\$|(?<!\\)\$(?!\$).*?(?<!\\)\$',
        '行内代码': r'(`+)([^`]+)\1',
        '引用链接': r'https?://[^\s)\]>]+',
    }
    return {name: Counter(re.findall(pattern, text, re.S)) for name, pattern in patterns.items()}


def check_script(source, target):
    original = source.read_text(encoding='utf-8')
    translated = target.read_text(encoding='utf-8')
    compile(translated, str(target), 'exec')
    trees = [_WithoutDocstrings().visit(ast.parse(text, type_comments=True))
             for text in (original, translated)]
    assert ast.dump(trees[0]) == ast.dump(trees[1]), f'示例计算代码变化：{target}'
    assert protected_parts(documentation(original)) == protected_parts(documentation(translated)), (
        f'示例公式、行内代码或引用链接变化：{target}'
    )


def check():
    tracked = subprocess.check_output(
        ['git', 'ls-tree', '-r', '--name-only', 'HEAD', '--', 'examples'],
        cwd=ROOT, text=True,
    ).splitlines()
    assert tracked, '未找到 examples 下的已提交文件'
    scripts = resources = 0
    directories = set()
    for name in tracked:
        relative = Path(name).relative_to('examples')
        source, target = SOURCE / relative, TARGET / relative
        assert target.is_file(), f'缺少中文副本文件：{relative}'
        if relative.suffix == '.py':
            check_script(source, target)
            directories.add(relative.parent.as_posix())
            scripts += 1
        elif relative.as_posix() != 'README.md':
            assert source.read_bytes() == target.read_bytes(), f'配套文件变化：{relative}'
            resources += 1
    readme = (TARGET / 'README.md').read_text(encoding='utf-8')
    for directory in directories:
        assert f']({directory}/)' in readme, f'示例目录未列入中文 README：{directory}'
    print(f'检查通过：{scripts} 个中文示例脚本；{resources} 个原版配套文件；中文目录链接完整。')


if __name__ == '__main__':
    check()
