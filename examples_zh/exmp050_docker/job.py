#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockDocker
from opi.input.simple_keywords import Sqm
from opi.input.structures import Properties, Structure


def run_exmp050(
    structure: Structure | None = None, working_dir: Path | None = Path("RUN")
) -> tuple[list[Structure], list[Properties]]:
    # > 重新创建工作目录
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    # > 客体结构位于示例文件夹中
    example_folder = Path(__file__).parent

    # > 如未给定结构，则从 inp.xyz 读取
    if structure is None:
        structure = Structure.from_xyz("inp.xyz")

    # > 设置计算器
    calc = Calculator(basename="job", working_dir=working_dir)
    calc.structure = structure
    calc.input.add_simple_keywords(Sqm.GFN2_XTB)

    calc.input.add_blocks(BlockDocker(guest=example_folder / "inp.xyz"))

    # > 写出输入文件并运行计算
    calc.write_input()
    calc.run()

    # > 获取输出并检查部分结果
    output = calc.get_output()
    if not output.terminated_normally():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    # << IF 语句结束

    # > 解析 JSON 文件
    output.parse()

    # > 输出结果结构图，以查看可用结果
    output.print_graph()

    structures = Structure.from_trj_xyz(
        working_dir / f"{calc.basename}.docker.struc1.all.optimized.xyz"
    )

    properties_list = Properties.from_trj_xyz(
        working_dir / f"{calc.basename}.docker.struc1.all.optimized.xyz", mode="docker"
    )

    # > 输出读取的结构
    for structure, properties in zip(structures, properties_list):
        print(f"FINAL ENERGY: {properties.energy_total}")
        print(structure.to_xyz_block())

    return structures, properties_list


if __name__ == "__main__":
    output = run_exmp050()

