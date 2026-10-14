#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Dft, Scf, SolvationModel, Solvent
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp007(
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
        Scf.NOAUTOSTART, Dft.BP86, BasisSet.DEF2_SVP, SolvationModel.CPCM(Solvent.WATER)
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

    ngeoms = len(output.results_properties.geometries)
    print("N GEOMETRIES")
    print(ngeoms)
    print("DFT ENERGY")
    print(output.results_properties.geometries[0].dft_energy.finalen)
    print("Solvent")
    print(output.results_properties.geometries[0].solvation_details.solvent)
    print("Epsilon")
    print(output.results_properties.geometries[0].solvation_details.epsilon)
    print("Surface points")
    print(output.results_properties.geometries[0].solvation_details.npoints)
    print("CPCM ENERGY")
    print(output.results_properties.geometries[0].solvation_details.cpcmdielenergy)

    return output


if __name__ == "__main__":
    run_exmp007()

