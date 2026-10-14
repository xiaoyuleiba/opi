#!/usr/bin/env python3
"""
示例：使用 `Grepper` 回退机制进行 GFN-FF 几何优化。

GFN-FF 是 ORCA 中的外部力场方法（内部调用 xtb）。
它会写出性质 JSON 文件，但其中不包含几何结构条目，
因此常规的 JSON 输出解析会对能量、梯度和结构返回 None。
默认启用的 `Grepper` 回退机制会直接从纯文本 .out 文件中
提取这些量。
"""

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import ForceField, Task
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp055(
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
    # > GFN-FF 为外部方法；ORCA 写出 .out 文件，但不写出性质 JSON
    calc.input.add_simple_keywords(ForceField.GFN_FF, Task.OPT)

    calc.write_input()
    calc.run()

    # > 获取输出并禁用版本检查，因为外部方法无法进行该检查
    output = calc.get_output(version_check=False)
    if not output.terminated_normally():
        print(f"ORCA calculation failed, see output file: {output.get_outfile()}")
        print(output.error_message())
        sys.exit(1)

    # > 对于 GFN-FF，parse() 加载的性质 JSON 不包含几何结构条目
    output.parse()

    # > 当没有可用的 JSON 数据时，get_final_energy / get_structure / get_gradient 会自动回退到从 .out
    # > 文件中检索数据（fallback=True 为默认值）
    final_energy = output.get_final_energy()
    print("FINAL SINGLE POINT ENERGY (from .out fallback)")
    print(final_energy)

    optimized = output.get_structure()
    print("OPTIMIZED STRUCTURE (from .out fallback)")
    print(optimized.to_xyz_block())

    # > 几何优化不会为最终几何结构输出梯度。
    # > 使用 index=-2 获取倒数第二步的梯度。
    gradient = output.get_gradient(index=-2)
    print("GRADIENT AT SECOND-TO-LAST STEP in Eh/Bohr (from .out fallback)")
    print(gradient)

    # > 回退机制支持逐步遍历优化过程。
    # > 索引超出优化步数时，get_final_energy 返回 None。
    step = 0
    print("ENERGIES ALONG OPTIMISATION TRAJECTORY")
    while True:
        e = output.get_final_energy(index=step)
        if e is None:
            break
        s = output.get_structure(index=step)
        print(f"  step {step:3d}  E = {e:.10f}  ({len(s.atoms)} atoms)")
        step += 1

    return output


if __name__ == "__main__":
    run_exmp055()

