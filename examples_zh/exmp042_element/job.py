#!/usr/bin/env python3

import shutil
from pathlib import Path

from opi.core import Calculator
from opi.input.structures import Structure
from opi.output.core import Output
from opi.utils.element import Element


def run_exmp042(
    structure: Structure | None = None, working_dir: Path | None = Path("RUN")
) -> Output:
    # > 重新创建工作目录
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    # > 如未给定结构，则从 inp.xyz 读取
    if structure is None:
        structure = Structure.from_xyz("inp.xyz")

    # > 禁用版本检查，以免依赖 ORCA
    calc = Calculator(basename="job", working_dir=working_dir, version_check=False)
    calc.structure = structure

    # > 输出输入文件中各原子的原子序数
    for atom in calc.structure.atoms:
        print(atom.element.atomic_number)

    # > 输出其他元素的原子序数
    he_element = Element("he")
    print(he_element.atomic_number)

    return calc.get_output()


if __name__ == "__main__":
    run_exmp042()

