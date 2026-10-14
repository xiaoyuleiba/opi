#!/usr/bin/env python3
"""检查中文文档与英文源码的结构和计算内容是否一致。"""
from collections import Counter
import ast
import io
import json
from pathlib import Path
import re
import subprocess
import tokenize

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'docs'
TARGET = ROOT / 'docs_zh'


class _WithoutDocstrings(ast.NodeTransformer):
    """仅移除模块、函数和类的文档字符串，以比较实际代码。"""

    def _strip(self, node):
        self.generic_visit(node)
        if ast.get_docstring(node) is not None:
            node.body.pop(0)
        return node

    visit_Module = _strip
    visit_FunctionDef = _strip
    visit_AsyncFunctionDef = _strip
    visit_ClassDef = _strip


def python_code(text):
    tree = ast.parse(text, type_comments=True)
    return ast.dump(_WithoutDocstrings().visit(tree), include_attributes=False)


def script_comments(text):
    tokens = tokenize.generate_tokens(io.StringIO(text).readline)
    return '\n'.join(token.string for token in tokens if token.type == tokenize.COMMENT)


def protected_parts(text):
    """提取不应由翻译改变的代码块、公式与行内标识符。"""
    patterns = {
        'code blocks': r'```.*?```',
        'mathematics': r'(?<!\\)\$\$.*?(?<!\\)\$\$|(?<!\\)\$(?!\$).*?(?<!\\)\$',
        'inline code': r'(?<!`)`(?!`)(.*?)`(?!`)',
    }
    return {name: Counter(re.findall(pattern, text, re.S)) for name, pattern in patterns.items()}


def check():
    tracked = subprocess.check_output(
        ['git', 'ls-tree', '-r', '--name-only', 'HEAD', '--', 'docs'], cwd=ROOT, text=True
    ).splitlines()
    notebooks = 0
    resources = 0
    scripts = 0
    for name in tracked:
        relative = Path(name).relative_to('docs')
        source = SOURCE / relative
        target = TARGET / relative
        assert target.is_file(), f'缺少文件：{relative}'
        if relative.suffix == '.ipynb':
            original = json.loads(source.read_text(encoding='utf-8'))
            translated = json.loads(target.read_text(encoding='utf-8'))
            assert len(original['cells']) == len(translated['cells']), relative
            assert {k: v for k, v in original.items() if k != 'cells'} == {
                k: v for k, v in translated.items() if k != 'cells'
            }, f'Notebook 元数据变化：{relative}'
            for index, (a, b) in enumerate(zip(original['cells'], translated['cells'])):
                label = f'{relative}，单元 {index}'
                if a['cell_type'] != 'markdown':
                    assert a == b, f'代码或输出变化：{label}'
                else:
                    assert {k: v for k, v in a.items() if k != 'source'} == {
                        k: v for k, v in b.items() if k != 'source'
                    }, f'Markdown 单元元数据变化：{label}'
                    assert protected_parts(''.join(a['source'])) == protected_parts(''.join(b['source'])), (
                        f'公式、代码块或标识符变化：{label}'
                    )
            notebooks += 1
        elif relative.parent == Path('contents/scripts') and relative.suffix == '.py':
            original = source.read_text(encoding='utf-8')
            translated = target.read_text(encoding='utf-8')
            compile(translated, str(target), 'exec')
            assert python_code(original) == python_code(translated), f'脚本代码变化：{relative}'
            a, b = script_comments(original), script_comments(translated)
            assert protected_parts(a) == protected_parts(b), f'脚本公式或代码标识符变化：{relative}'
            links = r'https?://[^\s)\]>]+'
            assert Counter(re.findall(links, a)) == Counter(re.findall(links, b)), (
                f'脚本引用链接变化：{relative}'
            )
            scripts += 1
        elif relative.suffix != '.md' and relative.name not in {'conf.py', 'Makefile'}:
            assert source.read_bytes() == target.read_bytes(), f'资源或脚本变化：{relative}'
            resources += 1
    catalog = json.loads((TARGET / '_translations/api_lines.json').read_text(encoding='utf-8'))
    assert catalog and all(isinstance(k, str) and isinstance(v, str) for k, v in catalog.items())
    print(f'检查通过：{notebooks} 个 Notebook；{scripts} 个中文注释脚本；'
          f'{resources} 个原版资源；{len(catalog)} 条 API 译文。')


if __name__ == '__main__':
    check()

