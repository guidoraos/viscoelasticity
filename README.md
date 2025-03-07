# viscoelasticity

This directory contains two C++ program that can be used to exract viscoelastic properties from atomistic MD simulationf of liquids:
1) "autocorr.cpp" computes the autocorrelation functions (e.g. of the components of the stress tensor) by recursive coarse-graining in time. 
2) "G1G2.cpp" computes the in-phase and out-of-phase components of the shear modulus [G'(omega), G"(omega)] and complex viscosity [eta'(omega), eta"(omega)] from the autocorrelation function of the stress tensor (shear relaxation modulus).

THESE PROGRAMS ARE PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED. IT IS THE USER'S RESPONSABILITY TO ENSURE THAT CALCULATIONS BASED ON THEM ARE EXECUTED, VALIDATED AND INTERPRETED CORRECTLY.

If you publish work based on this code, please cite:
David, A., De Nicola, A., Tartaglino, U., Milano, G. and Raos, G., 2019.
Viscoelasticity of Short Polymer Liquids from Atomistic Simulations.
Journal of The Electrochemical Society, 166(9), p.B3246.
https://iopscience.iop.org/article/10.1149/2.0371909jes
