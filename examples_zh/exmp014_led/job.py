#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockMdci
from opi.input.simple_keywords import Approximation, AuxBasisSet, BasisSet, Dlpno, Scf, Wft
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp014(
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
        Wft.DLPNO_CCSD_T,
        BasisSet.CC_PVDZ,
        AuxBasisSet.CC_PVDZ_C,
        AuxBasisSet.CC_PVTZ_JK,
        Dlpno.TIGHTPNO,
        Scf.TIGHTSCF,
        Approximation.RIJK,
        Dlpno.LED,
        Dlpno.ADLD,  # 伦敦色散能的原子分解
        Dlpno.ADEX,  # 交换能的原子分解
    )

    # > 使用不同的电荷方案进行原子分解
    calc.input.add_blocks(BlockMdci(ad_loewdin=True, ad_mulliken=True))

    calc.write_input()
    calc.run()

    output = calc.get_output()
    if not output.terminated_normally():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    if not output.scf_converged():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    if not output.cc_converged():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)

    # > 解析 JSON 文件
    output.parse()

    # > 获取结构，并使用相同的片段 ID 写出新的输入文件
    new_structure = output.get_structure()
    new_calc = Calculator(basename="new_job", working_dir=working_dir)
    new_calc.structure = new_structure
    new_calc.input.add_simple_keywords(
        Wft.DLPNO_CCSD_T,
        BasisSet.CC_PVDZ,
        AuxBasisSet.CC_PVDZ_C,
        AuxBasisSet.CC_PVTZ_JK,
        Dlpno.TIGHTPNO,
        Scf.TIGHTSCF,
        Approximation.RIJK,
        Dlpno.LED,
    )

    new_calc.write_input()
    return output


if __name__ == "__main__":
    run_exmp014()

