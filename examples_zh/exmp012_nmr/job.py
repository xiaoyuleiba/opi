#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import (
    BlockEprnmr,
    Nuclei,
    NucleiFlag,
)
from opi.input.simple_keywords import (
    AuxBasisSet,
    BasisSet,
    Property,
    SolvationModel,
    Solvent,
    Wft,
)
from opi.input.structures import Structure
from opi.output.core import Output
from opi.utils.element import Element


def run_exmp012(
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
        Wft.RIMP2,
        BasisSet.PCSSEG_2,
        AuxBasisSet.AUTOAUX,
        Property.NMR,
        SolvationModel.CPCM(Solvent.CHLOROFORM),
    )

    calc.input.add_blocks(
        BlockEprnmr(
            nuclei=Nuclei(atom=Element.HYDROGEN, flags=(NucleiFlag(shift=True, ssall=True)))
        )
    )

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
    run_exmp012()

