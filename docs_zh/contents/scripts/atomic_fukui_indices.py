#!/usr/bin/env python
# coding: utf-8

# # Atomic Fukui Indices
# 
# This tutorial demonstrates how atomic Fukui indices are calculated using the ORCA python interfacte (OPI).
# 
# In this notebook we will:
# 1. Import the Required Dependencies
# 2. Define a Working Directory and Prepare Structures
# 3. Visualize the Initial Molecular Structure
# 4. Run ORCA Calculations
# 5. Charge Extraction and Fukui Function Computation
# 6. 3D Visualization of Fukui Functions

# ## Step 1: Import Dependencies
# 
# We start by importing the modules needed for:
# - Interfacing with ORCA input/output
# - Numerical calculations and data handling
# - Plotting results
# - Handling for directory
# 
# > **Note:** We additionally import modules for visualization/plotting like `py3Dmol`. For this, it might be necessary to install `py3Dmol` into your OPI `venv` (e.g., by activating the `.venv` and using `uv pip install py3Dmol`).

# In[ ]:


# > Import pathlib for directory handling
from pathlib import Path
import shutil

# > Import numpy for data handling and numerical operations
import numpy as np

# > OPI imports for performing ORCA calculations and reading output
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import BasisSet, Dft, DispersionCorrection
from opi.input.structures.structure import Structure

# > Import libraries for visualization
import py3Dmol
from IPython.display import display, HTML


# ## Step 2: Define Working directory
# 
# All actual calculations will be performed in a subfolder `RUN`.

# In[ ]:


# > Calculation is performed in `RUN`
working_dir = Path("RUN")
# > The `working_dir`is automatically (re-)created
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()


# ## Step 3: Setup an Input Structure
# 
# We use formaldehyde as our example molecule. The 3D structure in Cartesian coordinates is defined in XYZ format and visualized.

# In[3]:


# > Define the molecule's Cartesian coordinates in Angstroem as python string 
xyz_data = """\
4

C      0.59761    0.92710    0.00015
O      1.78965    0.92659    0.00044
H      0.00000    1.85448    0.00000
H      0.00000    0.00000    0.00000\n
"""

# > Visualize the molecular structure using py3Dmol
view = py3Dmol.view(width=500, height=500)
view.addModel(xyz_data, "xyz")
view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})
view.setBackgroundColor('white')
view.zoomTo()
view.show()

# > Save the XYZ structure to a file
with open(working_dir / "struc.xyz","w") as f:
    f.write(xyz_data)


# ## Step 4: ORCA Calculations
# 
# We perform single-point energy calculations for three electronic states of the molecule:
# - Neutral (default)
# - Anion (additional electron, charge = -1, multiplicity = 2)
# - Cation (electron removed, charge = +1, multiplicity = 2)

# The successful calculation has to be checked and parsed, which is done in `check_and_parse_output`:

# In[4]:


def check_and_parse_output(output: Output)-> None:
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


# We start with the neutral calculation:

# In[5]:


# > Path to structure
xyz_file = working_dir / "struc.xyz"

# > Create a Calculator object for ORCA input generation and execution
calc = Calculator(basename="neutral", working_dir=working_dir)

# > Load the molecular structure from XYZ file
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 0
calc.structure.multiplicity = 1

# > Add calculation keywords
calc.input.add_simple_keywords(
    Dft.B3LYP,
    BasisSet.DEF2_SVP,
    DispersionCorrection.D3
)

# > Write ORCA input file and run the calculation
calc.write_input()
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")

# > Get the output object
output_neutral = calc.get_output()
check_and_parse_output(output_neutral)


# Next, the calculation for the anion is performed:

# In[6]:


calc = Calculator(basename="anion", working_dir=working_dir)
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = -1
calc.structure.multiplicity = 2
calc.input.add_simple_keywords(
    Dft.B3LYP,
    BasisSet.DEF2_SVP,
    DispersionCorrection.D3
)
calc.write_input()
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")
output_anion = calc.get_output()
check_and_parse_output(output_anion)


# Finally, the calculation of the cation is performed:

# In[7]:


calc = Calculator(basename="cation", working_dir=working_dir)
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 1
calc.structure.multiplicity = 2
calc.input.add_simple_keywords(
    Dft.B3LYP,
    BasisSet.DEF2_SVP,
    DispersionCorrection.D3
)
calc.write_input()
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")
output_cation = calc.get_output()
check_and_parse_output(output_cation)


# ## Step 5: Charge Extraction and Fukui Function Computation

# The Fukui function is computed in three common forms:
# 
# - f⁺ (nucleophilic attack): ρ(N+1) - ρ(N)
# - f⁻ (electrophilic attack): ρ(N) - ρ(N-1)
# - f⁰ (radical attack): ½[ρ(N+1) - ρ(N-1)]
# 
# Here, we employ Mulliken charges for calculating atomic Fukui indices (see [J. Am. Chem. Soc. 1986, 108, 19, 5708–5711](https://doi.org/10.1021/ja00279a008) for details). Note that employing Mulliken charges instead of populations/densities introduces a minus sign in all equations above.

# In[8]:


def compute_fukui_charge(charges1: list[float], charges2: list[float], mode: str) -> list[float]:
    """Compute the Fukui function from atomic charges."""
    q1 = np.array(charges1)
    q2 = np.array(charges2)

    # > Reverted order due to using charges instead of populations/densities
    if mode == 'plus/minus':
        diff_charge = (q2 - q1)
    elif mode == 'zero':
        diff_charge = (q2 - q1) / 2
    return diff_charge.tolist()

# > Extract atomic charges from Mulliken population analysis
charges_neutral = [c[0] for c in output_neutral.results_properties.geometries[-1].mulliken_population_analysis[-1].atomiccharges]
charges_anion = [c[0] for c in output_anion.results_properties.geometries[-1].mulliken_population_analysis[-1].atomiccharges]
charges_cation = [c[0] for c in output_cation.results_properties.geometries[-1].mulliken_population_analysis[-1].atomiccharges]

# > Compute atomic Fukui indices based on Mulliken charges, note the negative sign
# > since atomic charges are used instead of the density
f_plus_charges = compute_fukui_charge(charges1=charges_anion, charges2=charges_neutral, mode='plus/minus')
f_minus_charges = compute_fukui_charge(charges1=charges_neutral, charges2=charges_cation, mode='plus/minus')
f_zero_charges = compute_fukui_charge(charges1=charges_anion, charges2=charges_cation, mode='zero')


# ## Step 6: 3D Visualization of Fukui Functions

# In[9]:


def read_xyz(xyz_path: str) -> tuple[list[tuple[str, float, float, float]], str]:
    """Read an XYZ file and return atom info and raw string."""
    with open(xyz_path, "r") as f:
        lines = f.readlines()
    num_atoms = int(lines[0])
    atom_lines = lines[2:2 + num_atoms]
    coords = []
    for line in atom_lines:
        parts = line.strip().split()
        atom = parts[0]
        x, y, z = map(float, parts[1:])
        coords.append((atom, x, y, z))
    return coords, "".join(lines)

def create_molecule_view(xyz_file: str, title: str, charges: list[float] = None) -> str:
    """Generate an interactive 3D molecular visualization using py3Dmol."""
    atoms, xyz_str = read_xyz(xyz_file)
    view = py3Dmol.view(width=300, height=300)
    view.addModel(xyz_str, "xyz")
    view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})

    for (atom, x, y, z), charge in zip(atoms, charges):
        view.addLabel(f"{charge:.3f}", {
            'position': {'x': x, 'y': y, 'z': z},
            'backgroundColor': 'black',
            'fontColor': 'white',
            'showBackground': True,
            'fontSize': 12
        })
    title += " - Label"

    view.zoomTo()
    view.zoom(1.8)
    view.setBackgroundColor('white')
    html = view._make_html()

    return f"""
    <div style="display:inline-block; text-align:center; margin:10px;">
        <h4>{title}</h4>
        {html}
    </div>
        """

# Define XYZ path
xyz_path = working_dir / "struc.xyz"

# Generate visualizations
label_views = [
    create_molecule_view(xyz_path, title="Mulliken charges", charges=charges_neutral),
    create_molecule_view(xyz_path, title="F+", charges=f_plus_charges),
    create_molecule_view(xyz_path, title="F-", charges=f_minus_charges),
    create_molecule_view(xyz_path, title="F0", charges=f_zero_charges),
]

# Display label views in 2x2 grid
label_html = f"""
<h3>Atomic Fukui Indices (Labels)</h3>
<div style='display: flex; flex-wrap: wrap; justify-content: left; gap: 10px;'>
    {label_views[0]} {label_views[1]}
    {label_views[2]} {label_views[3]}
</div>
"""
display(HTML(label_html))


# The Fukui indices suggest that a nucleophilic attack would occur at the carbon atom, an electrophilic attack would occur at the oxygen atom, and a radical attack would also occur at the oxygen atom.

# ## Summary
# 
# In this notebook, we calculated and visualized atomic fukui indices based on Mulliken charges with the ORCA Python Interface (OPI). Other charge schemes are available (Loewdin, Mayer) and could be used as well.
