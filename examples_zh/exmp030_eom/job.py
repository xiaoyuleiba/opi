#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockMdci
from opi.input.simple_keywords import BasisSet, Wft
from opi.input.structures import Structure
from opi.output.core import Output


# > 进行 EOM-CCSD/def2-SVP 计算，求前五个激发态
def run_exmp030(
    structure: Structure | None = None, working_dir: Path | None = Path("RUN")
) -> Output:
    # > 重新创建工作目录
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    # > 如未给定结构，则从 inp.xyz 读取
    if structure is None:
        structure = Structure.from_xyz("inp.xyz")

    calc = Calculator(basename="job", working_dir=working_dir)
    calc.structure = structure
    calc.input.add_simple_keywords(Wft.EOM_CCSD, BasisSet.DEF2_SVP)

    calc.input.add_blocks(BlockMdci(nroots=5))

    calc.write_input()
    calc.run()

    output = calc.get_output()
    if not output.terminated_normally():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    # << IF 语句结束

    # > 解析 JSON 文件
    output.parse()

    return output


if __name__ == "__main__":
    run_exmp030()

