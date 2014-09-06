antenna-array-optimization
==========================

An optimization framework for constructing linear dipole antenna arrays written in MATLAB. 
Cost function for the optimization problem is defined as the combined reflection coefficient of the array.
Array structure can be illustrated as follows;

==TL==D1==TL==D2== . . . ==TL==Dn
TL: Transmission line   D: Dipole

Cost function is defined in antenna_cost.m. Parameters can be changed in this file. Input parameters of this function are
dipole lengths, TL lengths and TL impedances.
impfits.mat file is required for operation as it consists of a curve fit data which replaces complex integral calculations. Fit data can be created by using impedance_curvefit.m.
A supplementary script, newton_method.m, is included to perform the optimization process. It is slower compared to native
toolbox function of MATLAB, however, process flow can be tracked.
