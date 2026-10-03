#!/usr/bin/env python
# coding: utf-8

# # Properties for ML models
# 
# This tutorial demonstrates how to use ORCA together with the ORCA python interface (OPI) to calculate electronic properties relevant for, e.g., training machine-learned interatomic potentials (MLIPs). As an example, we will use the reference level &omega;B97M-V/def2-TZVPD used for generating reliable training data, e.g., for [OMol25](https://arxiv.org/abs/2505.08762). As an exemplary molecule, we use 4-methylpyrazole.
# 
# In this notebook we will:
# 1.  Import the required python dependencies.
# 2.  Define a working directory.
# 3.  Set up the input structure.
# 4.  Define the calculation settings.
# 5.  Perform the ORCA calculation.
# 6.  Check for proper termination and convergence.
# 7.  Process the results.

# ## Step 1: Import Dependencies
# 
# We start by importing the modules needed for:
# - Interfacing with ORCA input/output
# - Numerical calculations and data handling
# - Plotting results
# 

# In[1]:


from pathlib import Path
import shutil

# > pandas and numpy for data handling
import pandas as pd
import numpy as np

# > OPI imports for performing ORCA calculations and reading output
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, \
    Scf, AuxBasisSet, Approximation, Dft, Task, OutputControl, Grid
from opi.input.structures.structure import Structure

# > for visualization of molecules
import py3Dmol


# ## Step 2: Working directory
# 
# We define a subfolder `RUN` in which the actual ORCA calculations will take place.

# In[ ]:


# > Calculation is performed in `RUN`
working_dir = Path("RUN")
# > The `working_dir`is automatically (re-)created
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()


# ## Step 3: Setup the Input Structure
# 
# We define the structure with Cartesian coordinates in angstrom and visualize it:
# 

# In[3]:


# > define cartesian coordinates in angstrom as python string 
xyz_data = """\
12

N  1.3237  0.6845 -0.0002
N  1.4406 -0.6544 -0.0001
C -0.7492 -0.0043  0.0001
C  0.0327  1.1249  0.0002
C -2.2189 -0.0810 -0.0002
C  0.1711 -1.0697  0.0001
H -0.2057  2.1791  0.0003
H -2.6522  0.6920 -0.6433
H -2.5652 -1.0527 -0.3670
H -2.6091  0.0584  1.0131
H -0.0275 -2.1321  0.0003
H  2.1699  1.2401 -0.0003\n
"""
# > Visualize the input structure
view = py3Dmol.view(width=400, height=400)
view.addModel(xyz_data, 'xyz')
view.setStyle({}, {'stick': {}, 'sphere': {'scale': 0.3}})
view.zoomTo()
view.show()

# > Write the input structure to a file
with open(working_dir / "struc.xyz","w") as f:
    f.write(xyz_data)
# > Read structure into object
structure = Structure.from_xyz(working_dir / "struc.xyz")


# ## Step 4: Define the Calculation Settings
# Next, we have to define the settings used in the ORCA calculation. This can be done by defining a calculator and adding simple keywords to it:

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 4) -> Calculator:
    # > Set up a Calculator object, the basename for the ORCA calculation is set to LED
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > Assign structure to calculator
    calc.structure = structure

    # define a level of theory: wB97M-V/def2-TZVPPD with simple
    sk_list = [
        Dft.WB97M_V, # > wB97M-V method
        BasisSet.DEF2_TZVPD, # > def2-TZVPD basis set
        AuxBasisSet.DEF2_J, # > auxiliary basis set
        Scf.TIGHTSCF, # > tight SCF
        Approximation.RIJCOSX, # > Use RIJ and COSX
        Grid.DEFGRID3, # > Large grid for exchange-correlation and COSX
        Task.ENGRAD, # > Calculate energy and gradient
        OutputControl.PRINTGAP, # > Print HOMO-LUMO gap
    ]

    # use simple keywords in calculator
    calc.input.add_simple_keywords(*sk_list)

    # Define number of CPUs for the calcualtion
    calc.input.ncores = ncores # > 4 CPUs for this ORCA run

    return calc

calc = setup_calc("single_point", working_dir=working_dir, structure=structure)


# ## Step 5: Perform the ORCA calculation
# 
# To perform the calculation, we have to write the input in ORCA format first, and then run the job by using the `run()` function of the calculator class. The output of ORCA can be written to the calculator with its `get_output()` function.

# In[5]:


def run_calc(calc: Calculator) -> Output:
    # > Write the ORCA input file
    calc.write_input()
    # > Run the ORCA calculation
    print("Running ORCA calculation ...", end="")
    calc.run()
    print("   Done")

    # > Get the output object
    output = calc.get_output()

    return output

output = run_calc(calc)


# ## Step 6: Check for Proper Termination and Convergence
# When performing automated computations, a termination and convergence check should always be included and, when computing larger data sets, failed convergence should be handled. A simple check could look like:

# In[6]:


def check_and_parse_output(output: Output):
    # > Check for proper termination of ORCA
    status = output.terminated_normally()
    if not status:
        # > ORCA did not terminate normally
        raise RuntimeError("ORCA did not terminate normally. Please check RUN/single_point.out")
    else:
        # > ORCA did terminate normally so we can parse the output
        output.parse()

    # Now check for convergence of the SCF
    if output.results_properties.geometries[0].single_point_data.converged:
        print("SCF CONVERGED")
    else:
        raise RuntimeError("SCF DID NOT CONVERGE")

check_and_parse_output(output)


# ## Step 7: Process the Results
# After a successful calculation, we can access the results and store them. For example:

# In[7]:


def process_output(output: Output) -> dict:
    # > Conversion factor for Hartree -> eV
    h_to_eV = 27.2114
    # > Conversion factor for Bohr -> Angstrom
    bohr_to_A = 0.529177

    # > Energy and gradient information
    # > Total energy in eV
    Etot = output.results_properties.geometries[0].single_point_data.finalenergy*h_to_eV
    # > VV10 energy part in eV
    # > Note that VV10 is self-consistent in gradient computations and thus ecnl is used
    E_VV10 = output.results_properties.geometries[0].dft_energy.ecnl*h_to_eV
    # > Gradient in eV/Bohr (list of size n_atoms containing one list per atom containing three entries)
    Grad = output.results_properties.geometries[0].nuclear_gradient[0].grad
    Grad = [[x * (h_to_eV/bohr_to_A) for x in nuc_grad] for nuc_grad in Grad]

    # > General information
    # > Charge of molecule (float value)
    tot_chrg = output.results_properties.calculation_info.charge
    # > Spin of molecule (float value)
    tot_spin = output.results_properties.calculation_info.mult
    # > Number of atoms (float value)
    n_at = output.results_properties.calculation_info.numofatoms
    # > Number of electrons (float value)
    n_el = output.results_properties.calculation_info.numofelectrons
    # > Number of basis functions (float value)
    n_bf = output.results_properties.calculation_info.numofbasisfuncts
    # > Number of alpha electrons (float value)
    n_el_alp = output.results_properties.geometries[0].dft_energy.nalphael
    # > Number of beta electrons (float value)
    n_el_beta = output.results_properties.geometries[0].dft_energy.nbetael

    # > Population analysis
    # > Loewdin charges (list of size n_atoms containing lists with 3 entries)
    chrgs_lowedin = output.results_properties.geometries[0].loewdin_population_analysis[0].atomiccharges
    # > Mulliken charges (list of size n_atoms containing lists with 3 entries)
    chrgs_mulliken = output.results_properties.geometries[0].mulliken_population_analysis[0].atomiccharges
    # > Mayer bond order (list of size of the number of bonds)
    bond_order = output.results_properties.geometries[0].mayer_population_analysis[0].bondorders
    # > Mayer's total valence (list of size n_atoms)
    va_mayer = output.results_properties.geometries[0].mayer_population_analysis[0].va

    # > Dipole moments
    # > Total dipole moments (list with 3 entries: x, y, z)
    dip_tot = output.results_properties.geometries[0].dipole_moment[0].dipoletotal
    # > Magnitude of dipole (float value)
    dip_mag = output.results_properties.geometries[0].dipole_moment[0].dipolemagnitude
    # > Electronic dipole contributions (list with 3 entries: x, y, z)
    dip_el = output.results_properties.geometries[0].dipole_moment[0].dipoleeleccontrib
    # > Nuclear dipole contributions (list with 3 entries: x, y, z)
    dip_nu = output.results_properties.geometries[0].dipole_moment[0].dipolenuccontrib

    # > MO handling
    mos = output.results_gbw.molecule.molecularorbitals.mos # > saves a list of MOs with size n_atoms
    mo_energies = []
    mo_occupancy = []
    for mo in mos:
        mo_energies.append(mo.orbitalenergy) # > Extract the MO energies
        mo_occupancy.append(mo.occupancy) # > Extract the occupation numbers
    mo_energies = np.array(mo_energies) # > Make float list
    mo_occupancy = np.array(mo_occupancy) # > Make float list
    if np.all(np.isin(mo_occupancy, [1.0, 0.0])):
        print("Unrestricted calculation detected.")
    elif np.all(np.isin(mo_occupancy, [2.0, 0.0])):
        print("Restricted calculation detected.")

    HOMO_energy = np.max(mo_energies[(mo_occupancy == 1.0) | (mo_occupancy == 2.0)]) # > Detect HOMO energy
    HOMO_energy = HOMO_energy * h_to_eV
    LUMO_energy = np.min(mo_energies[mo_occupancy == 0.0]) # > Detect LUMO energy
    LUMO_energy = LUMO_energy * h_to_eV

    summary_data = {
    "Total Energy [eV]": [Etot],
    "Gradient [ev/A]": [str(Grad[0][0])[:10]+', '+str(Grad[1][0])[:10]+', '+str(Grad[2][0])[:10]],
    "VV10 Energy [eV]": [E_VV10],
    "Dipole Magnitude": [dip_mag],
    "Total Charge": [tot_chrg],
    "Spin Multiplicity": [tot_spin],
    "Number of Atoms": [n_at],
    "Number of Electrons": [n_el],
    "Alpha Electrons": [n_el_alp],
    "Beta Electrons": [n_el_beta],
    "Basis Functions": [n_bf],
    "HOMO Energy [eV]": [HOMO_energy],
    "LUMO Energy [eV]": [LUMO_energy],
    }

    return summary_data

summary_data = process_output(output)


# In the following, the resulting values are shown with pandas:

# In[8]:


df_summary = pd.DataFrame(summary_data, index=[""])
df_summary.T

