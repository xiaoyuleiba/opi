#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import BlockMethod
from opi.input.simple_keywords import (
    BasisSet,
    DispersionCorrection,
    Method,
    Scf,
    SolvationModel,
    Solvent,
    Task,
)
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp024(
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
        Scf.NOAUTOSTART,
        Method.HF,
        BasisSet.DEF2_SVP,
        Task.SP,
        SolvationModel.CPCM(Solvent.WATER),
        DispersionCorrection.D3,
    )

    # 检查是否已添加指定的简单关键词
    keywordlist = [Scf.NOAUTOSTART, Task.FREQ, "GOAT", "noAutoStart"]
    checks = calc.input.has_simple_keywords(*keywordlist)
    for keyword, check in zip(keywordlist, checks):
        print(f"keyword {keyword} has {'been added' if check else 'not been added'}")

    # 添加输入块，获取并修改后，再添加回 Calculator
    calc.input.add_blocks(BlockMethod(d3s6=0.64, d3a1=0.3065, d3s8=0.9147, d3a2=5.0570))
    method_block = calc.input.get_block("method")
    # 修改 %method 块，该块以其 ORCA 块名称作为键
    method_block.d3s6 = 0.4
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

    print("FINAL SINGLE POINT ENERGY")
    print(output.results_properties.geometries[0].single_point_data.finalenergy)
    return output


if __name__ == "__main__":
    run_exmp024()

