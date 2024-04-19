# Description
ElectroMagnetic SIMulator (emsim): Package for a multitude of electromagnetic simulations.

**Capabilities**
- Absorptance, reflectance and transmittance through 1D multilayered stacks
- Any arbitrary combination of thin (coherent) and thick (incoherent) layers
- EQE, Jsc, Vsc and Pmax of adjacent semiconductor solar cell
- Lab values reflected/transmitted light

**Features**
- Wavelengths and thicknesses are all given is nanometers

# Benchmark
## T-matrix in thin films
For checking the validity of our code, one can use simple [online calculator](http://people.ece.umn.edu/users/taylo589/tmm.shtml) and our `ex_pln_bnch.mxl` script

# Examples

- `ex_pln_01.mlx` calculates ATR spectra according to Figure 2 in Scalora et al, see folder [papers](https://github.com/iliasundensity/emsim/tree/main/papers) for pdf;
- `ex_pln_02.mlx` calculates reflectance spectra of double-sided PSC (Bragg-glass-Bragg) stack in different situations;
- `ex_pln_03.mlx` calculates ATR and Lab colors for low-E slide;
- `ex_pln_04.mlx` calculates transmittance and electrical characteristics of industry standard (DSM/Glass/EVA/PV) and SunDensity (PSC/EVA/PV) design.
