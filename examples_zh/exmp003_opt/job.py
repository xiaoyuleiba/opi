#!/usr/bin/env python3

import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import Dft, Scf, Task
from opi.input.structures import Structure
from opi.output.core import Output
from opi.output.models.json.gbw.gbw_results import GbwResults


def run_exmp003(
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
    calc.input.add_simple_keywords(Scf.NOAUTOSTART, Dft.WB97X3C, Task.OPT)
    calc.input.ncores = 4

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

    # > 检查 SCF 是否收敛
    if not output.scf_converged():
        print(f"ORCA SCF failed to converge, see output file: {output.get_outfile()}")
        sys.exit(1)

    # > 检查几何优化是否收敛
    if not output.geometry_optimization_converged():
        print(
            f"ORCA geometry optimization failed to converge, see output file: {output.get_outfile()}"
        )
        sys.exit(1)

    ngeoms = len(output.results_properties.geometries)
    print("N GEOMETRIES")
    print(ngeoms)
    print("FINAL SINGLE POINT ENERGY")
    print(output.results_properties.geometries[-1].single_point_data.finalenergy)
    print("SCF Energy along trajectory")
    # > 几何结构索引从 1 到 *ngeom*
    for igeom in range(0, ngeoms):
        print(
            f"{igeom})", output.results_properties.geometries[igeom].single_point_data.finalenergy
        )
    print("Mulliken charges along trajectory")
    # > 几何结构索引从 1 到 *ngeom*
    for igeom in range(0, ngeoms):
        try:
            charges = (
                output.results_properties.geometries[igeom]
                .mulliken_population_analysis[0]
                .atomiccharges
            )
        except TypeError:
            charges = "n/a"
        print(f"{igeom})", charges)

    # > 现在将最终结构输出为 XYZ 文件
    optimized = output.get_structure()
    print(optimized.to_xyz_block())

    # > 最终结构也保存在 GBW 文件中，因此可从该文件获取
    # > 而非从性质 JSON 中获取。
    from_gbw = output.get_structure_from_gbw()
    print("RMSD between the structure from the property JSON and from the gbw JSON")
    print(optimized.rmsd(from_gbw))

    # > 即使仅持有 GBW 文件而没有 `Output` 对象，也可以这样做。
    # > 使用 `orca_2json` 即时将 GBW 文件转换为 JSON。
    standalone = GbwResults.from_gbw_file(working_dir / "job.gbw").get_structure()
    print("STRUCTURE FROM GBW FILE")
    print(standalone.to_xyz_block())

    # > 现在输出最后一次计算的梯度，它对应于
    # > 最终结构之前一步的结构
    print(output.get_gradient(index=-2))

    return output


if __name__ == "__main__":
    run_exmp003()

