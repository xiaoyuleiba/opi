#!/usr/bin/env python3

import shutil
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockGeom
from opi.input.simple_keywords import BasisSet, Method, Scf, Task
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp026(
    structure: Structure | None = None, working_dir: Path | None = Path("RUN")
) -> Output:
    # > 重新创建工作目录
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    # > 如未给定结构，则从 inp.xyz 读取
    if structure is None:
        structure = Structure.from_xyz("inp.xyz")

    calc_bond = Calculator(basename="job", working_dir=working_dir)
    calc_bond.structure = structure
    calc_bond.input.add_simple_keywords(
        Scf.NOAUTOSTART,
        Method.HF,
        BasisSet.MINIX,
        Task.OPT,
    )

    calc_bond.input.add_blocks(BlockGeom(scan="B  0  1 = 1.1, 2.1, 5"))
    calc_bond.write_input()
    calc_bond.run()
    output_bond = calc_bond.get_output()
    output_bond.parse()

    # > 输出扫描过程中的 HOMO-LUMO 能隙
    for index, gbw in enumerate(output_bond.results_gbw[1:], start=1):
        print(index, output_bond.get_hl_gap(index))

    # > 绘制扫描过程中的分子轨道
    for index, gbw in enumerate(output_bond.results_gbw[1:], start=1):
        homo_id = output_bond.get_homo().index
        cube_output = output_bond.plot_mo(homo_id, gbw_index=index)
        print(index, cube_output)

    # > 绘制扫描过程中的密度
    for index, gbw in enumerate(output_bond.results_gbw[1:], start=1):
        cube_output = output_bond.plot_density(gbw_index=index)
        print(index, cube_output)

    # > 访问扫描过程中的能量
    for index, gbw in enumerate(output_bond.results_properties.geometries[1:], start=1):
        print(index, output_bond.get_final_energy(index=index))

    # > 通过柔性势能面扫描改变键角
    calc_angle = Calculator(basename="scan_angle", working_dir=working_dir)
    calc_angle.structure = structure
    calc_angle.input.add_simple_keywords(
        Scf.NOAUTOSTART,
        Method.HF,
        BasisSet.DEF2_SVP,
        Task.OPT,
    )

    calc_angle.input.add_blocks(BlockGeom(scan="A  1 0 2 = 110, 170, 6"))
    calc_angle.write_input()
    calc_angle.run()

    # > 通过柔性势能面扫描改变二面角
    calc_dihedral = Calculator(basename="scan_dihedral", working_dir=working_dir)
    calc_dihedral.structure = structure
    calc_dihedral.input.add_simple_keywords(
        Scf.NOAUTOSTART,
        Method.HF,
        BasisSet.DEF2_SVP,
        Task.OPT,
    )

    calc_dihedral.input.add_blocks(BlockGeom(scan="D 6 1 0 2 = 60, 120, 6"))
    calc_dihedral.write_input()
    calc_dihedral.run()

    return output_bond

