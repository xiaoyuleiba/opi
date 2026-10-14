#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockTddft
from opi.input.simple_keywords import BasisSet, Dft, SolvationModel, Solvent
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp010(
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
    calc.input.add_simple_keywords(
        Dft.B3LYP, BasisSet.DEF2_TZVP, SolvationModel.CPCM(Solvent.HEXANE)
    )

    calc.input.add_blocks(BlockTddft(nroots=31))

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

    # > 现在输出本次计算中可用的能量类型
    energies = output.get_energies()
    print("The following energies can be obtained from this calculation:")
    for key in energies:
        print(f"Energy of type {energies[key].method}")

    return output


if __name__ == "__main__":
    run_exmp010()

