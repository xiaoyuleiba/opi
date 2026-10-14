#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Method, Scf, Task
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp047(working_dir: Path | None = Path("RUN")) -> list[Output]:
    example_folder = Path(__file__).parent
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    # > 从 inp.xyz 读取结构
    structures = Structure.from_trj_xyz(example_folder / "inp_trj.xyz")
    print(f"Number of structures in inp_trj.xyz: {len(structures)}")
    # > 从 other.xyz 读取结构（中间以空行分隔）
    structures_blank_lines = Structure.from_trj_xyz(
        example_folder / "with_blank_lines_trj.xyz", comment_symbols=tuple("\n")
    )
    print(f"Number of structures in with_blank_lines_trj.xyz: {len(structures_blank_lines)}")
    # > 从 another.xyz 读取结构（XYZ 块之间用 > 行分隔）
    structures_comments = Structure.from_trj_xyz(
        example_folder / "with_comments_trj.xyz", comment_symbols=tuple(">")
    )
    print(f"Number of structures in with_comments_trj.xyz: {len(structures_comments)}")

    output_list = []
    for index, structure in enumerate(structures):
        calc = Calculator(basename="job", working_dir=working_dir)
        calc.structure = structure
        calc.input.add_simple_keywords(Scf.NOAUTOSTART, Method.HF, BasisSet.DEF2_SVP, Task.SP)

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

        # 检查 SCF 是否收敛
        if not output.scf_converged():
            print("SCF DID NOT CONVERGE")
            sys.exit(1)

        print(index, output.get_final_energy())

        output_list.append(output)

    return output_list


if __name__ == "__main__":
    run_exmp047()

