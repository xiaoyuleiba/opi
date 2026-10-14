#!/usr/bin/env python3
"""
示例：执行 GFN2-xTB 过渡态优化，随后进行频率计算。

从接近过渡态的几何结构出发，OPTTS 定位一阶鞍点，
随后 FREQ 计算 Hessian 矩阵。真正的过渡态恰有一个虚频（负频率），
可通过 ``get_imaginary_frequencies()`` 获取。
"""

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import (
    Opt,
    Scf,
    Sqm,
    Task,
)
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp056(
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
    calc.input.add_simple_keywords(Scf.NOAUTOSTART, Sqm.NATIVE_GFN2_XTB, Opt.OPTTS, Task.NUMFREQ)

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

    print("ALL VIBRATIONAL FREQUENCIES [cm-1]")
    frequencies = output.get_frequencies()
    for mode, freq in frequencies.items():
        print(f"  Mode {mode:3d}: {freq:12.4f}")

    print("IMAGINARY FREQUENCIES [cm-1]")
    imaginary = output.get_imaginary_frequencies()
    if imaginary:
        for mode, freq in imaginary.items():
            print(f"  Mode {mode:3d}: {freq:12.4f}")
    else:
        print("  None")

    print(f"Number of imaginary frequencies: {len(imaginary)}")
    print("Is PES transition state")
    print(output.is_pes_transition_state())
    print("FINAL SINGLE POINT ENERGY")
    print(output.get_final_energy())

    return output


if __name__ == "__main__":
    run_exmp056()

