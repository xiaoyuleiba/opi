#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockSolvator
from opi.input.simple_keywords.solvent import Solvent
from opi.input.structures.structure import Structure
from opi.output.core import Output


def run_exmp036(
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
    calc.structure.charge = 1
    calc.structure.multiplicity = 1
    calc.input.ncores = 16
    calc.input.add_simple_keywords("XTB")
    calc.input.add_simple_keywords(f"ALPB({Solvent.WATER})")
    calc.input.add_blocks(BlockSolvator(nsolv=1))

    calc.write_input()
    calc.run()

    output = calc.get_output()
    if not output.terminated_normally():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)

    output.parse(read_gbw_json=False)

    print("Solvated xyz file:")
    print(f"\t{output.basename}.solvator.xyz")

    # > 检查输出是否正常
    # > solvator 运行可提供的输出信息较少，
    # > 主要是 solvator.xyz、solvator.solventbuild.xyz 或最终单点能

    return output


if __name__ == "__main__":
    run_exmp036()

