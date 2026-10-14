import shutil
import sys
from pathlib import Path

from opi.core import Calculator
from opi.input.simple_keywords import SolvationModel, Solvent
from opi.input.structures import Structure
from opi.output.core import Output


def run_exmp053(
    structure: Structure | None = None, working_dir: Path | None = Path("RUN")
) -> Output:
    # > 重新创建工作目录
    shutil.rmtree(working_dir, ignore_errors=True)
    working_dir.mkdir()

    # > 如未给定结构，则从 inp.xyz 读取
    if structure is None:
        structure = Structure.from_xyz("inp.xyz")

    # > 设置计算器
    calc = Calculator(basename="job", working_dir=working_dir)
    calc.structure = structure
    calc.input.add_simple_keywords(SolvationModel.COSMORS(Solvent.ACETONITRILE))

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

    # > 解析 JSON 文件
    output.parse()

    return calc.get_output()


if __name__ == "__main__":
    output = run_exmp053()

