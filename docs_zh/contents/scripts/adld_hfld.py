#!/usr/bin/env python
# coding: utf-8

# # ADLD - HFLD
# 
# This tutorial demonstrates how to use ORCA together with the ORCA python interface (OPI) to perform the **[Atomic Decomposition of London Dispersion energy (ADLD)](https://doi.org/10.1021/acscentsci.5c00356)** analysis to decompose the London dispersion energy of the **[Hartree-Fock plus London Dispersion (HFLD)](https://doi.org/10.1021/acs.jctc.9b00425)** interaction energy within a dimer into atomic contributions. Decomposing the HFLD energy allows a non-empirical estimate of atomic London dispersion contributions. For the quicker but more empirical ADLD(D4) decomposition, see the corresponding notebook. In comparison to ADLD(D4) only one calculation of the whole system is required, since the subsystems consisiting of only one fragment have no contribution.
# 
# In this notebook we will:
# 1.  Import the required python dependencies.
# 2.  Define a working directory.
# 3.  Set up the input structure.
# 4.  Run a HFLD calculation for the supersystem.
# 7.  Process the ADLD results.
# 8.  Visualize the results.
# 9.  Generating the Cube File.

# ## Step 1: Import Dependencies
# 
# We start by importing the modules needed for:
# - Interfacing with ORCA input/output
# - Numerical calculations and data handling
# - Plotting results
# 
# > **Note:** We additionally import modules for visualization/plotting like `py3Dmol`. For this, it might be necessary to install `py3Dmol` into your OPI `venv` (e.g., by activating the `.venv` and using `uv pip install py3Dmol`).

# In[1]:


from pathlib import Path
import shutil

# > pandas and numpy for data handling
import numpy as np

# > OPI imports for performing ORCA calculations and reading output
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, Dlpno, \
    Scf, AuxBasisSet, Approximation
from opi.input.structures.structure import Structure
from opi.input.blocks import BlockFrag

# > For visualization of molecules
import py3Dmol


# ## Step 2: Working Directory and Conversion Factor
# 
# We define a subfolder `RUN` in which the actual ORCA calculations will take place. Also, we define a conversion factor, since we want the resulting interaction energies in kcal/mol for better interpretability.

# In[ ]:


# > Calculation is performed in `RUN`
working_dir = Path("RUN")
# > The `working_dir`is automatically (re-)created
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()
# > Conversion factor for atomic units to kcal/mol
unit_conversion = 627.509 


# ## Step 3: Setup the Input Structure
# 
# As an example we will decompose the LD interaction in an cubane dimer. The 3D structure in Cartesian coordinates is defined and visualized:
# 

# In[3]:


# > Define cartesian coordinates in Angstroem as python string 
# > In principle, this can be replaced with other neutral NCI dimers
xyz_data = """\
32
cubane dimer
C -2.37438447 -0.01391040 0.09884370
H -1.28799640 -0.03512861 0.15212364
C -4.14778752 -0.93512832 -0.83041188
C -4.12750489 1.22660964 -0.39560995
C -3.20643343 0.25290551 -1.19176459
C -4.23643748 -0.22901780 1.25693113
C -3.31541388 -1.20282062 0.46084844
C -3.29511100 0.95905424 0.89568028
C -5.06903241 0.03775498 -0.03358762
H -4.49182155 -1.69883200 -1.52468842
H -4.45488733 2.20468479 -0.74189213
H -2.79034514 0.44620851 -2.17829474
H -4.65224637 -0.42236391 2.24356374
H -2.98771769 -2.18079631 0.80710684
H -2.95068094 1.72266937 1.58984614
H -6.15558078 0.05880427 -0.08645048
C 2.05734830 0.00538088 -0.06151765
H 0.97096030 0.02659937 -0.11479860
C 2.97807526 -0.96758450 -0.85835285
C 2.99837839 1.19429062 -0.42352221
C 3.91940243 0.22048707 -1.21960352
C 2.88939603 -0.26143458 1.22909153
C 3.81046792 -1.23513946 0.43293829
C 3.83075080 0.92659878 0.86773902
C 4.75199612 -0.04628527 0.07091613
H 2.63364560 -1.73119990 -1.55251858
H 2.67068285 2.17226622 -0.76978149
H 4.33521233 0.41383252 -2.20623583
H 2.47330675 -0.45473693 2.21562140
H 4.13784970 -2.21321450 0.77922136
H 4.17478444 1.69030273 1.56201542
H 5.83854442 -0.06733485 0.12378001

\n
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


# ## Step 4: Run HFLD for the Supramolecular System
# 
# We run a HFLD/def2-TZVP calculation for the complete system. We choose the moderate def2-TZVP basis set, because the HFLD London dispersion energy is dependent on the basis and this is a good compromise for computational costs vs. accuracy. The fragments in the calculation will automatically be defined by ORCA (in this example two fragments, one for each monomer). We define and run three helper functions: `setup_calc`, `run_calc`, and `check_and_parse_output`. 
# 
# First, we set up the calculation with `setup_calc`:

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 4) -> Calculator:
    # > Set up a Calculator object
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > Assign structure to calculator
    calc.structure = structure

    # > Define a simple keyword list for the ORCA calculation.
    sk_list = [
    Dlpno.HFLD,
    BasisSet.DEF2_TZVP,
    AuxBasisSet.DEF2_TZVP_C,
    Scf.NORMALSCF,
    Approximation.RIJCOSX,
    Dlpno.ADLD,
    Dlpno.LED
    ]

    # > Use simple keywords in calculator
    calc.input.add_simple_keywords(*sk_list)

    # > Automatically assign fragments by setting the frag block
    calc.input.add_blocks(BlockFrag(storefrags=True))

    # Define number of CPUs for the calcualtion
    calc.input.ncores = ncores # > 4 CPUs for this ORCA run

    return calc

calc_supersystem = setup_calc("supersystem", working_dir=working_dir, structure=structure)


# Then, we run the calculation with `run_calc` and obtain the output:

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

output_supersystem = run_calc(calc_supersystem)


# The successful calculation has to be checked and parsed, which is done in `check_and_parse_output`:

# In[6]:


def check_and_parse_output(output: Output):
    # > Check for proper termination of ORCA
    status = output.terminated_normally()
    if not status:
        # > ORCA did not terminate normally
        raise RuntimeError(f"ORCA did not terminate normally, see output file: {output.get_outfile()}")
    else:
        # > ORCA did terminate normally so we can parse the output
        output.parse()

    # Now check for convergence of the SCF
    if output.results_properties.geometries[0].single_point_data.converged:
        print("ORCA terminated normally and the SCF converged")
    else:
        raise RuntimeError("SCF DID NOT CONVERGE")

check_and_parse_output(output_supersystem)


# ## Step 5: Processing of ADLD Energies
# 
# Now we perform the post-processing for the ADLD. First, we gather atomic contributions of the LD. 

# In[7]:


# > Path to atomic contributions to the LD in the dimer (product) [Eh]
product_vdw = output_supersystem.results_properties.geometries[0].mdci_adld.adlddispatomic_loewdin

# > Energies are already relative since reactant contributions are zero
# > Apply conversion factor for result in kcal/mol
atomic_disp: list[float] = []
for i in range(len(product_vdw)):
    product_x = float(product_vdw[i][0])
    product_x_kcal = (product_x) * unit_conversion
    atomic_disp.append(product_x_kcal)


# Now we have an array of atomic contributions to the relative London dispersion interaction energies. A London dispersion density ($\rho_{disp}$) can be defined for easily visualizaing and analyzing the atomic contributions:

# In[8]:


def rho_disp(r, alpha, epsilons, positions):
    """
    Calculate the dispersion density rho_disp at position r.

    Parameters:
    - r: (3,) array-like, the 3D point to evaluate at.
    - alpha: float, parameter alpha (typically set to 0.3-0.5).
    - epsilons: (N,) array-like, dispersion strengths for each atom A.
    - positions: (N, 3) array-like, positions R_A of atoms.

    Returns:
    - rho: float, the dispersion density at r.
    """
    # Pre-factor
    prefactor = (np.pi / alpha)**(-1.5)

    # Sum over atoms
    r = np.array(r)
    positions = np.array(positions)
    deltas = r - positions  # shape (N, 3)
    squared_distances = np.sum(deltas**2, axis=1)  # shape (N,)

    contributions = epsilons * np.exp(-alpha * squared_distances)
    total_sum = np.sum(contributions)

    return prefactor * total_sum


# ## Step 6: Crude Visualization of ADLD Contributions
# 
# In the original ADLD publication, the `rho_disp` function was evaluated on a grid and the results was written to a `.cube` file for visualization. For the sake of a simple visualization within this notebook, we will first evaluate this function only on the atomic positions and color the atoms by the values at their centers. These plots offer a qualitative analysis of the data and may differ in detail from those published. Further below, we also generate the cube file for the full visualization. We refer the interested reader to the [LDDsuite](https://github.com/bistonigroup/LDDsuite) on GitHub.
# 
# We start by calculating $\rho_{disp}$ at the atomic positions:

# In[9]:


# > Path to the Cartesian coordinates of the dimer
elements = output_supersystem.results_properties.geometries[0].geometry.coordinates.cartesians

# > List of element symbols
elem_list: list[str] = []
# > List of atomic positions
positions: list[tuple[float,float,float]] = []

# > Fill the lists
for elem in elements:
    elem_list.append(elem[0])
    positions.append([elem[1],elem[2],elem[3]])

# > Calculate the London dispersion density (rho_disp) at atomic positions
rho_disp_atomic: list[float] = []
for i in range(len(positions)):
    # > alpha value of 0.3 is selected
    rho_disp_atomic.append(rho_disp(positions[i], 0.3, atomic_disp, positions))


# Next, we define a small helper function for mapping the values of the dispersion density to a color gradient:

# In[10]:


def value_to_color(value, min_val):
    """
    Maps a value (min_val to 0) to a white-to-red gradient.
    Numbers near 0 -> white; near min_val -> red.
    """
    if value >= 0:
        return '#FFFFFF'

    # Ensure min_val is negative and avoid division by zero
    if min_val == 0:
        min_val = -1e-6  # small epsilon to avoid division by zero

    # Normalize: 0 (white) to min_val (red)
    norm = min(abs(value / min_val), 1.0)

    r = 255
    g = int(255 * (1 - norm))
    b = int(255 * (1 - norm))

    return f'#{r:02X}{g:02X}{b:02X}'


# Finally, we plot $\rho_{disp}$ by adding colored spheres to the py3Dmol visualization. We require custom vdw radii for the sizes of the spheres to add:

# In[11]:


# > Custom vdw radii for visualization
vdw_radii = {
    'H': 1.20, 'He': 1.40, 'Li': 1.82, 'Be': 1.53,
    'B': 1.92, 'C': 1.70,  'N': 1.55,  'O': 1.52,
    'F': 1.47
}

# > Definition of maximum and minimum values for color gradient
min_val = -0.01 # > kcal/mol Bohr³

# > Setup the viewer
view = py3Dmol.view(width=400, height=400)
view.addModel(xyz_data)

# > Initial style: white stick representation for all atoms
view.setStyle({}, {'stick': {'color' : 'white'}})

# Apply spheres as overlays for color-coding
for i, val in enumerate(rho_disp_atomic):
    color = value_to_color(val, min_val)
    radius = vdw_radii[elem_list[i]] * 0.35 # > arbitrary scaling factor for nice spheres
    selection = {'serial': i}
    view.addStyle(selection, {'sphere': {'color': color, 'radius': radius}})

view.zoomTo()
view.show()


# ## Step 7: Generating the Cube File
# 
# For the correct visualization, also the cube file can be generated that might be plotted separately:

# In[12]:


# > Helper dict to obtain atomic numbers
element_dict = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9
}

# > Helper function for defining bounding box
def get_bounding_box(coords, padding):
    coords = np.array(coords)
    min_coords = coords.min(axis=0) - padding
    max_coords = coords.max(axis=0) + padding
    return min_coords, max_coords

# > Helper function for writing the actual cube file
def write_cube_file(filename, coords, atomic_numbers, values, min_coords, voxel_size, n_points):
    with open(filename, 'w') as f:
        f.write("Cube file generated by adld notebook\n")
        f.write("rho_disp\n")
        f.write(f"{len(coords):5d} {min_coords[0]:12.6f} {min_coords[1]:12.6f} {min_coords[2]:12.6f}\n")
        f.write(f"{n_points[0]:5d} {voxel_size[0]:12.6f} 0.0 0.0\n")
        f.write(f"{n_points[1]:5d} 0.0 {voxel_size[1]:12.6f} 0.0\n")
        f.write(f"{n_points[2]:5d} 0.0 0.0 {voxel_size[2]:12.6f}\n")
        for atom_num, (x, y, z) in zip(atomic_numbers, coords):
            f.write(f"{atom_num:5d} 0.0 {x:12.6f} {y:12.6f} {z:12.6f}\n")
        flat_values = values.flatten(order='F')  # Fortran order for cube files
        for i in range(0, len(flat_values), 6):
            line = "".join(f"{v:13.5e}" for v in flat_values[i:i+6])
            f.write(line + "\n")


# > Define the box size
padding = 3.0  # 3 Å padding around the molecule
min_coords, max_coords = get_bounding_box(positions, padding)

# > Define the grid
n_points = (40, 40, 40)
x_lin = np.linspace(min_coords[0], max_coords[0], n_points[0])
y_lin = np.linspace(min_coords[1], max_coords[1], n_points[1])
z_lin = np.linspace(min_coords[2], max_coords[2], n_points[2])
voxel_size = (max_coords - min_coords) / (np.array(n_points) - 1)

# Evaluate rho_disp on the grid
values = np.zeros(n_points)
for i, x in enumerate(x_lin):
    for j, y in enumerate(y_lin):
        for k, z in enumerate(z_lin):
            r = (x, y, z)
            # > we set alpha to 0.3
            values[i, j, k] = rho_disp(r, 0.3, atomic_disp, positions)

# > obtain atomic numbers
atomic_numbers = [element_dict[symbol] for symbol in elem_list]

# Write to cube file
write_cube_file(working_dir / "rho_disp.cube", positions, atomic_numbers, values, min_coords, voxel_size, n_points)
print("Cube file written: rho_disp.cube")


# ## Summary
# 
# In this notebook we demonstrated how the ORCA Python interface can be employed to decompose the London dispersion energy from HFLD into atomic contributions. We performed a crude visualization within this notebook and generated a cube file for the full analysis.
