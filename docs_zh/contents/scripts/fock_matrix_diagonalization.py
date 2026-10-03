#!/usr/bin/env python
# coding: utf-8

# # Fock Matrix Diagonalization
# 
# This tutorial demonstrates how to extract the Fock Matrix from an ORCA calculation using the ORCA python interfacte (OPI). After extraction the Fock matrix is diagonalized and Koopmans' theorem is used to estimate the ionization energy.
# 
# In this notebook we will:
# 1. Import the Required Dependencies
# 2. Define a Working Directory and Prepare Structures
# 3. Download and Prepare Molecular Structure from PubChem
# 4. Run ORCA Calculations
# 5. Extract Matrices and Perform Fock Matrix Diagonalization

# ## Step 1: Import Dependencies
# 
# We begin by importing all required Python modules for this tutorial. These include:
# 
# - **File and directory handling**: for managing input/output files.
# - **JSON parsing**: for reading ORCA-generated GBW JSON files.
# - **Numerical libraries**: for matrix operations, diagonalization, and eigenvalue problems.
# - **PubChem data retrieval**: for downloading molecular structures directly from PubChem using REST APIs.
# - **RDKit**: for molecule handling, conversion, and generating 3D geometries.
# - **ORCA Python Interface (OPI)**: for setting up, running, and parsing ORCA quantum chemistry calculations.
# - **py3Dmol**: for interactive 3D visualization of molecular structures directly in the notebook.
# 
# > **Note:** We additionally import modules for visualization/plotting like `py3Dmol`. For this, it might be necessary to install `py3Dmol` into your OPI `venv` (e.g., by activating the `.venv` and using `uv pip install py3Dmol`).

# In[1]:


# > Import pathlib for directory handling
from pathlib import Path
import shutil

# > Import library for parsing and handling JSON data
import json

# > Import libraries for data handling and numerical operations
import numpy as np
from scipy.linalg import eigh

# > Import libraries for retrieving molecular data from PubChem
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

# > OPI imports for performing ORCA calculations and reading output
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import BasisSet, Scf
from opi.input.simple_keywords.method import Method

# > Import py3Dmol for 3D molecular visualization
import py3Dmol


# ## Step 2: Define Working Directory
# 
# All actual calculations will be performed in a subfolder `RUN`.

# In[ ]:


# > Calculation is performed in `RUN`
working_dir = Path("RUN")
# > The `working_dir`is automatically (re-)created
shutil.rmtree(working_dir, ignore_errors=True)
working_dir.mkdir()
# > Conversion factor for atomic units to eV
unit_conversion = 27.211386


# ## Step 3: Download and Prepare Molecular Structure from PubChem
# 
# In this step, we retrieve the molecular structure of the target molecule directly from PubChem using its REST API.
# 

# In[3]:


# > Download molecule data from PubChem
mol_name = "Benzene"
cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{mol_name}/cids/JSON"
cid_response = requests.get(cid_url)
cid_response.raise_for_status()
cid = cid_response.json()['IdentifierList']['CID'][0]

sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
sdf_response = requests.get(sdf_url)
sdf_response.raise_for_status()
sdf_data = sdf_response.text

# > Parse SDF and generate 3D conformer with RDKit
mol = Chem.MolFromMolBlock(sdf_data)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()

# > Convert to XYZ format
num_atoms = mol.GetNumAtoms()
lines = [str(num_atoms), f"{mol_name} from PubChem"]
for atom in mol.GetAtoms():
    pos = conf.GetAtomPosition(atom.GetIdx())
    symbol = atom.GetSymbol()
    lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

xyz_data = "\n".join(lines)

# > Visualize the molecular structure using py3Dmol
view = py3Dmol.view(width=500, height=500)
view.addModel(xyz_data, "xyz")
view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.2}})
view.setBackgroundColor('white')
view.zoomTo()
view.show()

# Save the XYZ structure to a file
with open(working_dir / "struc.xyz","w") as f:
    f.write(xyz_data)


# ## Step 4: ORCA Calculations
# 
# Now, we perform a HF/STO-3G calculation. We set up a config_dict and use it to request integrals.  

# In[4]:


# > Set up path for the structure
xyz_file = working_dir / "struc.xyz"

# > Define GBW JSON configuration
config_dict = {
    "MOCoefficients": True,
    "1elPropertyIntegrals": ["dipole"],
    "1elIntegrals": ["H", "S"],
    "FockMatrix": ["F", "J", "K"],
    "JSONFormats": ["json"]
}

# > Create a Calculator object for ORCA input generation and execution
calc = Calculator(basename="benzene", working_dir=working_dir)

# > Load the molecular structure from XYZ file
structure = Structure.from_xyz(xyz_file)
calc.structure = structure
calc.structure.charge = 0
calc.structure.multiplicity = 1

# > Add calculation keywords
calc.input.add_simple_keywords(
            Method.HF,
            BasisSet.STO_3G,
            Scf.TIGHTSCF
)

# > Write the ORCA input file
calc.write_input()
# > Run the ORCA calculation
print("Running ORCA calculation ...", end="")
calc.run()
print("   Done")

# > Get output and use it to create the gbw json output with config
output = calc.get_output()
status = output.terminated_normally()
if status: 
    output.parse()
    output.create_gbw_json(force=True,config=config_dict)
else:
    raise RuntimeError("ORCA did not terminate normally.")


# ## Step 5: Extract Matrices and Perform Fock Matrix Diagonalization

# In[5]:


# Load the GBW JSON file
with open(working_dir / "benzene.json") as f:
    json_data = json.load(f)

# Extract integrals from the JSON file:
h = np.array(json_data['Molecule']['H-Matrix'])
f = np.array(json_data['Molecule']['F-Matrix'][0])
s = np.array(json_data['Molecule']['S-Matrix'])

# Build the full AO Fock matrix by adding core Hamiltonian and Fock correction
fao = h + f
# Perform generalized eigenvalue problem
eigenvalues, _ = eigh(fao, s)
# Convert orbital energies from Hartree to eV
eigenvalues_ev = eigenvalues * unit_conversion

# Estimate ionization energy using Koopmans' theorem
print("Ionization energy (eV):", -max([e for e in eigenvalues_ev if e < 0]))


# ## Summary
# 
# In this notebook we demonstrated how to utilize the ORCA Python Interface (OPI) to obtain Integrals from an ORCA calculation. We used this functionality to obtain the Fock matrix of a benzene molecule, diagonalized it with `eigh`, and used Koopmans' theorem to estimate the ionization energy.
