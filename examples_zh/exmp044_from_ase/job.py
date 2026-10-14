#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from ase import Atoms

from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Method, Scf, Task
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp044(working_dir: Path | None = Path("RUN")) -> Output:
    # 创建水分子（H2O）
    # 坐标以 Ångström 为单位
    positions = [
        [0.000, 0.000, 0.000],  # O
        [0.757, 0.586, 0.000],  # H
        [-0.757, 0.586, 0.000],  # H
    ]
    symbols = ["O", "H", "H"]

    # 创建 ASE Atoms 对象
    water = Atoms(symbols=symbols, positions=positions)

    # 设置初始电荷和磁矩
    # 示例：O 的电荷为 -1，H 为中性，因此总电荷为 -1
    water.set_initial_charges([-1.0, 0.0, 0.0])

    # 示例：指定 1 个未成对电子，对应二重态（S=1/2，2S=1）
    water.set_initial_magnetic_moments([1.0, 0.0, 0.0])

    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    calc = Calculator(basename="job", working_dir=working_dir)
    calc.structure = Structure.from_ase(water)
    calc.input.add_simple_keywords(Scf.NOAUTOSTART, Method.HF, BasisSet.DEF2_SVP, Task.SP)

    # > 输出电子数
    print(f"Number of electrons is {calc.structure.nelectrons}")
    print(f"Number of electrons is {'even' if calc.structure.nelec_is_even else 'odd'}")

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
    if output.results_properties.geometries[0].single_point_data.converged:
        print("SCF CONVERGED")
    else:
        print("SCF DID NOT CONVERGE")
        sys.exit(1)

    print("FINAL SINGLE POINT ENERGY")
    print(output.get_final_energy())

    return output


if __name__ == "__main__":
    run_exmp044()

