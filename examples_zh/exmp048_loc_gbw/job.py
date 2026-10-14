#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

import numpy as np

from opi.core import Calculator
from opi.input.simple_keywords import BasisSet, Dft, Scf, Task
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp048(
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
    calc.charge = 0
    calc.multiplicity = 1
    calc.input.add_simple_keywords(
        Scf.NOAUTOSTART,
        Dft.TPSS,
        BasisSet.DEF2_SVP,
        Task.SP,
    )

    # > 添加用于局域化的任意输入字符串
    calc.input.add_arbitrary_string("%loc\nLocMet PM\nend\n")

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
    if output.scf_converged():
        print("SCF CONVERGED")
    else:
        print("SCF DID NOT CONVERGE")
        sys.exit(1)

    # > 保存分子轨道（MO）
    mos = output.get_mos()["mo"]

    # > 获取积分
    cmo = np.array([mo.mocoefficients for mo in mos])
    s = output.get_int_overlap(recreate_json=True)
    h = output.get_int_hcore(recreate_json=True)
    f = output.get_int_f(recreate_json=True)
    fao = h + f

    smo = cmo @ s @ cmo.T

    fmo = cmo @ fao @ cmo.T

    # > 创建包含局域分子轨道（LMO）的 GBW JSON 文件并解析
    output.create_gbw_json(force=True, suffix=".loc")
    output.parse()

    # > 保存局域分子轨道（LMO）
    lmos = output.get_mos()["mo"]

    # > 将 Fock 矩阵变换到 LMO 基下
    clmo = np.array([mo.mocoefficients for mo in lmos])
    flmo = clmo @ fao @ clmo.T

    # > 将 MO 基下的重叠矩阵保存到文件
    np.savetxt(working_dir / "smo.txt", smo, fmt="%.4f")
    # > 将 AO 基下的 Fock 矩阵保存到文件
    np.savetxt(working_dir / "fao.txt", fao, fmt="%.4f")
    # > 将 MO 基下的 Fock 矩阵保存到文件
    np.savetxt(working_dir / "fmo.txt", fmo, fmt="%.4f")
    # > 将 LMO 基下的 Fock 矩阵保存到文件
    np.savetxt(working_dir / "flmo.txt", flmo, fmt="%.4f")


if __name__ == "__main__":
    run_exmp048()

