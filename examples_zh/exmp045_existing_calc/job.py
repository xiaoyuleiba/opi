#!/usr/bin/env python3

import sys
from pathlib import Path

from opi.core import Calculator
from opi.output.core import Output


def run_exmp045() -> Output:
    # > 获取并解析已有计算的输出
    example_folder = Path(__file__).parent
    working_dir = example_folder / "RUN"

    # > 与 'RUN' 目录中的已有计算使用相同的基本名称
    calc = Calculator(basename="job", working_dir=working_dir)
    output = calc.get_output()
    if not output.terminated_normally():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)
    # << IF 语句结束

    # > 解析 JSON 文件（parse 会根据 GBW 和 property.txt 文件创建这些文件）
    output.parse()

    # 检查 SCF 是否收敛
    if output.results_properties.geometries[0].single_point_data.converged:
        print("SCF CONVERGED")
    else:
        print("SCF DID NOT CONVERGE")
        sys.exit(1)

    print("FINAL SINGLE POINT ENERGY")
    print(output.get_final_energy())
    print(output.get_homo())
    return output


if __name__ == "__main__":
    run_exmp045()

