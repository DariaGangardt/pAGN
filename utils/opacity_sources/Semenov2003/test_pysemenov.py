"""
test_pysemenov.py

Authors: Alessandro A. Trani, Daria Gangardt

This file performs a simple test of the python interface for the Semenov2003 opacity fortran code

Check pysemenov_module.sh for instructions how to compile the interface,
and check generate_combined_tables.py for how to generate a table of opacities for tabulation
"""
import pysemenov

# This file follows the flow of the Fortran main
# except that input and ouput is via Python

# See opacitypy.f for description of the options

ross = True  # False for Planck opacities
model = 'nrm'  # Normal iron abundances
top = 'c'  # composite grains
shap = 'a'  # aggregate shape

rho = 1e-11  # gas density, g/cm^3
T = 1e2  # temperature, kelvin

kappa = pysemenov.compute_kappa(top=top,ross=ross,model=model,shap=shap,rho=rho,T=T)

print("kappa in cm^2/g", kappa) # cm^2/g
print("should be      ", 2.941995025939745)
