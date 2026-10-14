#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockNeb
from opi.input.simple_keywords import Neb, Scf, Sqm
from opi.input.structures import Properties, Structure
from opi.output.core import Output
from opi.utils.units import AU_TO_KCAL


def run_exmp039(working_dir: Path | None = Path("RUN")) -> Output:
    example_folder = Path(__file__).parent
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    calc = Calculator(basename="job", working_dir=working_dir)
    shutil.copy(example_folder / "prod.xyz", working_dir / "prod.xyz")
    calc.structure = Structure.from_xyz(example_folder / "reac.xyz")
    calc.input.add_simple_keywords(Scf.NOAUTOSTART, Sqm.NATIVE_GFN2_XTB, Neb.NEB_TS)

    calc.input.add_blocks(BlockNeb(neb_end_xyzfile="prod.xyz"))

    calc.write_input()
    calc.run()

    output = calc.get_output()
    if not output.terminated_normally() and output.neb_converged():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    # << IF 语句结束

    # > 解析 JSON 文件
    output.parse()

    N = len(output.results_gbw)
    print(f"N: {N}")
    # > 输出扫描过程中的 HOMO-LUMO 能隙
    print("Printing HOMO-LUMO gap along the minimal energy pathway.")
    for index, gbw in enumerate(output.results_gbw[1:], start=1):
        print(index, output.get_hl_gap(index))

    # > 输出能量
    print("Printing the final energies along the final TS-optimization trajectory.")
    for index, gbw in enumerate(output.results_properties.geometries[1:], start=1):
        print(index, output.get_final_energy(index=index))

    # > 读取最小能量路径（MEP）上的构型及其能量
    print("Printing the energies along the minimal energy pathway.")
    mep_file = working_dir / f"{calc.basename}_MEP_trj.xyz"
    structures = Structure.from_trj_xyz(mep_file)
    properties_list = Properties.from_trj_xyz(mep_file, mode="neb")

    # > 也可以改为读取 _MEP.allxyz 文件
    mep_file = output.get_file("_MEP.allxyz")
    structures = Structure.from_trj_xyz(mep_file, comment_symbols=">")
    properties_list = Properties.from_trj_xyz(mep_file, mode="neb", comment_symbols=">")

    # > 输出相对于首个构型的 MEP 能量
    energy_first = properties_list[0].energy_total
    for image, (structure, properties) in enumerate(zip(structures, properties_list)):
        energy_relative = (properties.energy_total - energy_first) * AU_TO_KCAL
        print(f"IMAGE {image + 1}: {properties.energy_total} Eh ({energy_relative:.2f} kcal/mol)")

    return output


if __name__ == "__main__":
    run_exmp039()

