# Reduced-Precision-IMC
---
Half and double floating-point precision Implict Monte Carlo implementations of the Su & Olson volume source benchmark problem written in Julia. 
Plotted IMC results of radiation and material energy densities for cases with and without isotropic scattering versus benchmark solutions are included.

# Dependencies
* Julia
## Julia Libraries
 * Random
 * Plots
 * DelimitedFiles

# How to run:
Choose either the half or double floating-point precision script.  
Change the inputs line at the top of the script to include the path to your desired input: 

`inputs = readdlm(raw"calc1SO.in")`

The script can then be executed. 
The output is plots of radiation and material energy densities as well as csv files tabulating those quantities at each timestep.

# Inputs
* calc1SO.in - Isotropic Scatter, tau_max = 10
* calc2SO.in - Isotropic Scatter, tau_max = 100
* calc3SO.in - No Scatter, tau_max = 10
* calc4SO.in - No Scatter, tau_max = 100
