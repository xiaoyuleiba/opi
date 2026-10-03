#!/usr/bin/env python
# coding: utf-8

# # ADLD - DFTD4
# 
# This tutorial demonstrates how to use ORCA together with the ORCA python interface (OPI) to perform the **[Atomic Decomposition of London Dispersion energy (ADLD)](https://doi.org/10.1021/acscentsci.5c00356)** analysis to decompose the **[DFT-D4](https://doi.org/10.1063/1.5090222)** London dispersion interaction energy within molecular dimers into atomic contributions. Decomposing the DFT-D4 energy allows a quick estimate of atomic London dispersion contributions. For the more elaborate decomposition within the LED scheme see the corresponding notebook.
# 
# In this notebook we will:
# 1.  Import the required python dependencies.
# 2.  Define a working directory.
# 3.  Set up the input structure.
# 4.  Run a DFT-D4 calculation for the supersystem.
# 5.  Assign the non-covalent subsystems automatically.
# 6.  Run calculations for the separate dimers.
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
    Scf, Approximation, Dft, DispersionCorrection
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
# As an example we will decompose the LD interaction in an adamantane dimer. The 3D structure in Cartesian coordinates is defined and visualized:
# 

# In[3]:


# > Define cartesian coordinates in Angstroem as python string 
# > In principle, this can be replaced with other neutral NCI dimers
xyz_data = """\
52
adamantane dimer
C -1.2512580945 0.0014447688 -4.5331335448
C -0.0001235489 1.2445251646 -2.7612423124
C 1.2481335316 -0.0013367292 -4.5311413071
C -0.0029283957 -1.2565268787 -2.7691903511
C -1.2503012639 -1.2486921591 -3.6528854863
C -0.0006594466 0.0020950040 -5.4138613643
C -1.2467375929 1.2467239390 -3.6461299866
C 1.2473963684 1.2436834160 -3.6447126890
C 1.2453783869 -1.2511955542 -3.6514603233
C -0.0026505360 -0.0089965776 -1.8844115499
H -2.1450794262 0.0043889516 -5.1648857415
H 0.0009877220 2.1356012499 -2.1257000212
H 2.1418366957 -0.0004543894 -5.1618030462
H -0.0038482437 -2.1525328115 -2.1403679585
H 2.1455059233 -1.2737565650 -3.0272449767
H 1.2692678550 -2.1503564250 -4.2769208835
H -0.0009339824 -0.8777868157 -6.0657312229
H -1.2754300701 -2.1478464803 -4.2782930880
H -2.1504389072 -1.2698005946 -3.0287598493
H -0.8848964964 -0.0102264960 -1.2357477040
H 0.0010895729 0.8848621166 -6.0618708738
H 0.8778642353 -0.0123636823 -1.2336945592
H -2.1479497683 1.2694405004 -3.0236737797
H 1.2698733945 2.1461993692 -4.2654872765
H -1.2663652760 2.1493731136 -4.2667886778
H 2.1487325387 1.2637698441 -3.0222209102
H 3.8360448211 0.0002950370 -6.3613914539
C 4.7298358432 0.0010844825 -6.9919263799
C 4.7303085767 1.2481697816 -7.8756345239
C 4.7330471126 -1.2464687062 -7.8748194344
C 5.9783164308 0.0029481456 -6.1088865644
C 5.9773801566 1.2523639250 -8.7598945071
H 4.7076765449 2.1490254113 -7.2526031839
H 3.8288296955 1.2686817563 -8.4978754295
C 5.9807806158 -1.2475248312 -8.7579867258
H 3.8320676507 -1.2689344965 -8.4978189406
H 4.7112191662 -2.1476227034 -7.2521085138
C 7.2287109083 0.0062218143 -6.9896685254
H 5.9796434423 -0.8795548254 -5.4604383280
H 5.9753753912 0.8830178962 -5.4572372534
C 7.2245631611 1.2521204877 -7.8754596256
C 5.9782475715 0.0016408752 -9.6407077164
H 5.9761534374 2.1461192284 -9.3917461569
C 7.2279643357 -1.2408795855 -7.8736092341
H 5.9837006737 -2.1426289557 -9.3878126810
H 8.1226008273 0.0075858975 -6.3578267144
H 7.2461857452 2.1548254344 -7.2550714782
H 8.1256650763 1.2722286312 -8.4982189724
H 6.8586688519 0.0016709201 -10.2918095343
H 5.0963929708 0.0001414723 -10.2900312493
H 8.1290236854 -1.2603983485 -8.4964381158
H 7.2516713723 -2.1413794735 -7.2504037938
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


# ## Step 4: Run DFTD4 for the Supramolecular System
# 
# First, we run a B3LYP-D4/MINIX calculation for the complete system. We choose the minimal basis set, because the D4 London dispersion energy is not dependent on the basis. The fragments in the calculation will automatically be defined by ORCA (in this example two fragments, one for each molecule). We define and run three helper functions: `setup_calc`, `run_calc`, and `check_and_parse_output`. 
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
    Dft.B3LYP,
    DispersionCorrection.D4,
    BasisSet.MINIX,
    Scf.NORMALSCF,
    Approximation.RIJCOSX,
    Dlpno.ADLD
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


# ## Step 5: Definition of Subsystems
# Next, we extract the automatic generated fragments and define the subsystems: 

# In[7]:


# > Shortcut to the fragments data
frags = output_supersystem.results_properties.geometries[0].geometry.fragments
# > Obtain the number of fragments in the led
nfrags = max(item[0] for item in frags)

# > Extract Fragment data from autofragmentation to a dict 
# > Assumes continuous fragment numbering

# > Generate a list of atoms in the fragments
atoms_in_frag = [[] for _ in range(nfrags)]
for index, value in enumerate(frags):
    frag_id = int(value[0])
    atoms_in_frag[frag_id-1].append(index)
    structure.atoms[index].fragment_id = frag_id

# > Define a list of atoms in the subsystems
atoms_in_sub = [] 
atoms_in_sub.append(atoms_in_frag[0])   # > First subsystem 
atoms_in_sub.append([atom for frag in atoms_in_frag[1:] for atom in frag])  # > Second subsystem 
nsub = len(atoms_in_sub)

# > Define a list of fragments in the subsystems
frag_in_sub = [[] for _ in range(nsub)]
for ifrag in range(nfrags):
    if ifrag == 0: frag_in_sub[0].append(ifrag) # > First fragment is the first subsystem
    if ifrag > 0: frag_in_sub[1].append(ifrag)  # > Second fragment is the second subsystem

# > Visualize the fragments 
print("Here you can view the fragment IDs:")
for atom in structure.atoms:
    x, y, z = atom.coordinates.coordinates
    view.addLabel(f"{atom.fragment_id}", {
        'position': {'x': x, 'y': y, 'z': z},
        'backgroundColor': 'black',
        'fontColor': 'white',
        'showBackground': True,
        'fontSize': 12
    })
view.show()

# > Set IDs to zero 
for atom in structure.atoms:
    atom.fragment_id = None


# ## Step 6: Run the Subsystem Calculations
# 
# The subsystem structures are extracted and calculations are performed for each subsystem. We utilize the previously defined helper functions:

# In[8]:


sub_outputs: list[Output] = []
# Loop over subsystems (there are two), obtain substructures, and perform calculations
for i in range(len(atoms_in_sub)):

    # > Prepare calculation for the subsystem 
    subsystem = i + 1
    sub_structure = calc_supersystem.structure.extract_substructure(atoms_in_sub[i])
    calc_subsystem = setup_calc(f"subsystem_{subsystem}", working_dir, sub_structure)

    # > Visualize the subsystem for which the calculation will run
    # > First remove fragment data for visualization
    for atom in calc_subsystem.structure.atoms: atom.fragment_id = None
    xyz_data_sub = calc_subsystem.structure.to_xyz_block()
    # > Visualize the subsystem 
    view_sub = py3Dmol.view(width=400, height=400)
    view_sub.addModel(xyz_data_sub, 'xyz')
    view_sub.setStyle({}, {'stick': {}, 'sphere': {'scale': 0.3}})
    view_sub.zoomTo()
    print(f"Visualizing subsystem {subsystem}")
    view_sub.show()

    # > Run the calculation
    output_sub = run_calc(calc_subsystem)

    # > Check and parse the output
    check_and_parse_output(output_sub)

    sub_outputs.append(output_sub)


# ## Step 7: Processing of ADLD Energies
# 
# Now, we perform the post-processing for the ADLD. Therefore, we gather the total atomic contributions of the LD and calculate relative atomic contributions by subtracting the contributions from the separate monomers (reactants) from the contributions of the combined dimer (product): 

# In[9]:


# > Path to total atomic contributions to the LD in the dimer (product) [Eh]
product_vdw = output_supersystem.results_properties.geometries[0].vdw_correction.vdw_atomic
# > Path to total atomic contributions to the LD in the monomers (reactants) [Eh]
reactant_vdw = sub_outputs[0].results_properties.geometries[0].vdw_correction.vdw_atomic + sub_outputs[1].results_properties.geometries[0].vdw_correction.vdw_atomic

# > Calculate relative atomic dispersion contributions
atomic_disp: list[float] = []
for i in range(len(product_vdw)):
    product_x = float(product_vdw[i][0])
    reactant_x = float(reactant_vdw[i][0])
    difference = (product_x - reactant_x) * unit_conversion
    atomic_disp.append(difference)


# Now, we have an array of atomic contributions to the relative London dispersion interaction energies. A London dispersion density ($\rho_{disp}$) can be defined for easily visualizaing and analyzing the atomic contributions:

# In[10]:


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


# ## Step 8: Crude Visualization of ADLD Contributions
# 
# In the original ADLD publication, the `rho_disp` function was evaluated on a grid and the results were written to a `.cube` file for visualization. For the sake of a simple visualization within this notebook, we will first evaluate this function only on the atomic positions and color the atoms by the values at their centers. These plots offer a qualitative analysis of the data and may differ in detail from those published. Further below, we also generate the cube file for the full visualization. We refer the interested reader to the [LDDsuite](https://github.com/bistonigroup/LDDsuite) on GitHub.
# 
# We start by calculating $\rho_{disp}$ at the atomic positions:

# In[11]:


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

# In[12]:


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

# In[13]:


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


# ## Step 9: Generating the Cube File
# 
# For the correct visualization, also the cube file can be generated that might be plotted separately:

# In[14]:


# > helper dict to obtain atomic numbers
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

# > helper function for defining bounding box
def get_bounding_box(coords, padding):
    coords = np.array(coords)
    min_coords = coords.min(axis=0) - padding
    max_coords = coords.max(axis=0) + padding
    return min_coords, max_coords

# > helper function for writing the actual cube file
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

# > Obtain atomic numbers
atomic_numbers = [element_dict[symbol] for symbol in elem_list]

# Write to cube file
write_cube_file(working_dir / "rho_disp.cube", positions, atomic_numbers, values, min_coords, voxel_size, n_points)
print("Cube file written: rho_disp.cube")


# ## Summary
# 
# In this notebook we demonstrated how the ORCA Python interface can be employed to decompose the London dispersion energy from DFT-D4 into atomic contributions. We performed a crude visualization within this notebook and generated a cube file for the full analysis.
