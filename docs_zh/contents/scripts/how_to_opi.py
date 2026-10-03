#!/usr/bin/env python
# coding: utf-8

# # How to OPI &ndash; A First Example
# 
# This tutorial demonstrates the basic use of ORCA together with the ORCA python interface (OPI). As an example, we will simply compute the energy of a water molecule. Therefore, we will write a script that performs the OPI calculation and access the results. Both can be done with OPI.
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
# The first thing to do for using OPI is to import the modules needed. We additionally import also the modules for visualization like `py3Dmol`. For this, it might be necessary to install `py3Dmol` into your OPI `venv` (e.g., by activating the `.venv` and using `uv pip install py3Dmol`).

# In[1]:


from pathlib import Path
import shutil

# > OPI imports for performing ORCA calculations and reading output
from opi.core import Calculator
from opi.output.core import Output
from opi.input.simple_keywords import Dft, Task
from opi.input.structures.structure import Structure

# > for visualization of molecules
import py3Dmol


# ## Step 2: Working Directory
# 
# When working with OPI, it is always a good idea to define a subfolder in which the actual ORCA calculations will take place. In this case, we set up a `RUN` directory

# In[2]:


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


# > Define cartesian coordinates in angstrom as python string 
xyz_data = """\
3

O      0.00000   -0.00000    0.00000
H      0.00000    0.96899    0.00000
H      0.93966   -0.23409    0.03434\n
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
# Next, we have to define the settings used in the ORCA optimization. This can be done by defining a calculator and adding simple keywords to it:

# In[4]:


def setup_calc(basename : str, working_dir: Path, structure: Structure, ncores: int = 1) -> Calculator:
    # > Set up a Calculator object, the basename for the ORCA calculation is also set
    calc = Calculator(basename=basename, working_dir=working_dir)
    # > Assign structure to calculator
    calc.structure = structure

    # Define a level of theory: we use r2SCAN-3c for geometry optimization
    sk_list = [
        Dft.R2SCAN_3C, # > r2SCAN-3c method (Comes with a predefined basis set)
        Task.SP # > Perform the singlepoint
    ]

    # > Use simple keywords in calculator
    calc.input.add_simple_keywords(*sk_list)

    # > Define number of CPUs for the calcualtion
    calc.input.ncores = ncores # > CPUs for this ORCA run (default: 1)

    return calc

calc = setup_calc("single_point", working_dir=working_dir, structure=structure)


# ## Step 5: Perform the Calculation
# 
# To perform the single-point calculation, we have to write the input in ORCA format first, and then run the job by using the `run()` function of the calculator class. The output of ORCA can be written to the calculator with its `get_output()` function.

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
# When performing automated computations, a termination and convergence check should always be included. A simple check could look like:

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


# ## Step 7: Access the Results
# After a successful calculation, we can access the resulting properties. They are generally available via the `output` object:

# In[7]:


# > Total energy in Hartee
Etot = output.results_properties.geometries[0].single_point_data.finalenergy

# > Print total energy
print("Final energy in Hartee: ",Etot)


# ## Summary
# 
# In this notebook we demonstrated how to perform a DFT single point energy calculation with the ORCA Python Interface (OPI). 
