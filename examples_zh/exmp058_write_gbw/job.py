#!/usr/bin/env python3
"""
示例：根据已解析的输出数据写出 GBW 文件，并将其用作 SCF 初猜。

先执行 Hartree-Fock 计算，再将其 GBW JSON 文件解析为 `GbwResults`。
`GbwResults.to_gbw_file()` 通过运行
`orca_2json <json-file> -gbw`，将这些数据重新转换为 GBW 文件。
写出的文件包含几何结构、基组和分子轨道，
因此可用作第二次计算的 SCF 初猜（`MORead`）。

注意，ORCA 此时只创建一个基础 GBW 文件，其中的信息少于
计算本身写出的 GBW 文件。
"""

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Method, Scf, Task
from opi.input.structures import Structure
from opi.output.core import Output
from opi.output.models.json.gbw.gbw_results import GbwResults


def run_exmp058(
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
    calc.input.add_simple_keywords(Scf.NOAUTOSTART, Method.HF, BasisSet.DEF2_SVP, Task.SP)

    calc.write_input()
    calc.run()

    output = calc.get_output()
    if not output.terminated_normally() and output.scf_converged():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    # << IF 语句结束

    # > 解析 JSON 文件
    output.parse()

    print("FINAL SINGLE POINT ENERGY")
    print(output.get_final_energy())

    # > 将本次计算的 GBW JSON 文件解析为 `GbwResults`
    results_gbw = output.results_gbw[0]
    print("NUMBER OF MOLECULAR ORBITALS")
    print(len(results_gbw.molecule.molecularorbitals.mos))

    # > 根据已解析的数据写出 GBW 文件。该文件必须创建在工作
    # > 目录中，因为 ORCA 仅接受其自身工作目录下的 `%moinp` 文件。
    gbw_file = results_gbw.to_gbw_file(working_dir / "from_results.gbw")
    print("GBW FILE WRITTEN FROM THE PARSED DATA")
    print(f"{gbw_file.name} ({gbw_file.stat().st_size} bytes)")

    # > 将写出的 GBW 文件用作第二次计算的 SCF 初猜
    calc_moread = Calculator(basename="job_moread", working_dir=working_dir)
    # > 从 GBW 文件获取结构
    calc_moread.structure = GbwResults.from_gbw_file(gbw_file).get_structure()
    calc_moread.input.add_simple_keywords(
        Scf.NOAUTOSTART, Scf.MOREAD, Method.HF, BasisSet.DEF2_SVP, Task.SP
    )
    calc_moread.input.moinp = gbw_file

    calc_moread.write_input()
    calc_moread.run()

    output_moread = calc_moread.get_output()
    if not output_moread.terminated_normally() and output.scf_converged():
        print(f"ORCA calculation failed, see output file: {output_moread.get_outfile()}")
        print(output_moread.error_message())
        sys.exit(1)
    # << IF 语句结束

    output_moread.parse()

    print("FINAL SINGLE POINT ENERGY WITH THE WRITTEN GBW FILE AS GUESS")
    print(output_moread.get_final_energy())

    return output_moread


if __name__ == "__main__":
    run_exmp058()

