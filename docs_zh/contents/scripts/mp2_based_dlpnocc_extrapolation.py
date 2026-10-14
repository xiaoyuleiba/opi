from pathlib import Path
from opi.core import Calculator
import shutil
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import (
    Wft,
    BasisSet,
    Dlpno,
    AuxBasisSet,
    RelativisticCorrection,
    ShellType,
    Scf,
)
from opi.input.blocks import BlockMdci, BlockMp2


def run_mp2_based_dlpnocc_extrapolation(structure: Structure, wd: Path = Path("RUN")):
    """
    采用基于 MP2 的校正方案，使 DLPNO-CCSD(T) 计算趋近完整对自然轨道空间的极限，
    依据文献 J. Chem. Theory Comput. 2023, 19, 13, 4023–4032。
    """

    # > HF 和 CC 类型
    HF_TYPE = ShellType.UHF
    CC_TYPE = Wft.DLPNO_CCSD_T1 # 或 Wft.DLPNO_CCSD_T

    # > 按需手动修改这些设置
    BASIS_SET = BasisSet.AUG_CC_PVDZ_DK
    REL = RelativisticCorrection.DKH2
    AUX_BASIS = AuxBasisSet.AUTOAUX
    DLPNO_THRESH = Dlpno.TIGHTPNO  # 或 Dlpno.NORMALPNO

    # > 用于 DLPNO-CCSD(T) 计算的 MDCI 块
    cc_mdci_block = BlockMdci()
    cc_mdci_block.usefulllmp2guess = False
    cc_mdci_block.maxiter = 100

    # > 用于 DLPNO-MP2 计算的 MP2 块
    mp2_block = BlockMp2()

    # > 使 DLPNO-MP2 设置与 DLPNO-CCSD(T) 设置一致
    # > Normal PNO 设置
    if DLPNO_THRESH == Dlpno.NORMALPNO:
        mp2_block.tcutpno = 3.33 * 10e-7
        mp2_block.tcutdo = 1 * 10e-2
        mp2_block.tcutmkn = 1 * 10e-3
    # > Tight PNO 设置
    elif DLPNO_THRESH == Dlpno.TIGHTPNO:
        mp2_block.tcutpno = 1 * 10e-7
        mp2_block.tcutdo = 5 * 10e-3
        mp2_block.tcutmkn = 1 * 10e-3
    else:
        raise ValueError("Unkown DLPNO_THRESH")

    # > DLPNO-CCSD(T1)
    dlpno_cc_calc = Calculator(basename="dlpno_ccsdt", working_dir=wd)
    dlpno_cc_calc.structure = structure
    dlpno_cc_calc.input.add_simple_keywords(
        HF_TYPE, CC_TYPE, REL, BASIS_SET, AUX_BASIS, DLPNO_THRESH
    )
    dlpno_cc_calc.input.add_blocks(cc_mdci_block)

    # > 写出输入文件并运行计算
    status = dlpno_cc_calc.write_and_run()
    if not status:
        raise RuntimeError("DLPNO-CCSD(T) calculation did not terminate normally!")

    # > 获取输出
    dlpno_cc_out = dlpno_cc_calc.get_output()
    dlpno_cc_out.parse()

    # > DLPNO-MP2 计算
    dlpno_mp2_calc = Calculator(basename="dlpno_mp2", working_dir=wd)
    dlpno_mp2_calc.structure = structure
    dlpno_mp2_calc.input.add_simple_keywords(
        HF_TYPE, Wft.DLPNO_MP2, REL, BASIS_SET, AUX_BASIS, Scf.MOREAD, Scf.NOITER
    )
    dlpno_mp2_calc.input.moinp = wd / "dlpno_ccsdt.gbw"
    dlpno_mp2_calc.input.add_blocks(mp2_block)

    # > 写出输入文件并运行计算
    status = dlpno_mp2_calc.write_and_run()
    if not status:
        raise RuntimeError("DLPNO-MP2 calculation did not terminate normally!")

    # > 获取输出
    dlpno_mp2_out = dlpno_mp2_calc.get_output()
    dlpno_mp2_out.parse()

    # > RI-MP2 计算
    ri_mp2_calc = Calculator(basename="ri_mp2", working_dir=wd)
    ri_mp2_calc.structure = structure
    ri_mp2_calc.input.add_simple_keywords(
        HF_TYPE, Wft.RIMP2, REL, BASIS_SET, AUX_BASIS, Scf.MOREAD, Scf.NOITER
    )
    ri_mp2_calc.input.moinp = wd / "dlpno_ccsdt.gbw"

    # > 写出输入文件并运行计算
    ri_mp2_calc.write_and_run()
    if not status:
        raise RuntimeError("RI-MP2 calculation did not terminate normally!")

    # > 获取输出
    ri_mp2_out = ri_mp2_calc.get_output()
    ri_mp2_out.parse()

    energy_dlpno_cc = dlpno_cc_out.get_final_energy()
    corr_energy_dlpno_mp2 = dlpno_mp2_out.get_energies()["MP2"].correnergy[0][0]
    corr_energy_ri_mp2 = ri_mp2_out.get_energies()["MP2"].correnergy[0][0]

    extrapolated_energy = energy_dlpno_cc + (corr_energy_ri_mp2 - corr_energy_dlpno_mp2)

    # > 返回最终外推能量
    return extrapolated_energy

if __name__ == "__main__":
    xyz_string = """3

    O         -3.56626        1.77639        0.00000
    H         -2.59626        1.77639        0.00000
    H         -3.88959        1.36040       -0.81444"""
    structure = Structure.from_xyz_block(xyz_string)
    wd = Path("RUN")
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir()
    energy = run_mp2_based_dlpnocc_extrapolation(structure, wd)
    print(f"FINAL EXTRAPOLATED ENERGY: {energy}")

