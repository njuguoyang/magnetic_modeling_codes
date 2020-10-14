# Magnetic Modeling Codes
FORTRAN and IDL codes to analyze solar magnetic field observations and construct magnetic models.
Some readme files can be found under different pathes.

These codes fit into the following topics, which have been described in [this paper](http://adsabs.harvard.edu/abs/2017ScChE..60.1408G).

1. Analysis of solar vector magnetic field observations
- Stokes profiles inversion
- Removing the 180° ambiguity of the transverse components of vector magnetic field
- Projection effect correction
- Optical flow techniques to derive velocity from time series of magnetic field observations

2. Magnetic field models
- Force-free field models including potential, linear force-free, and nonlinear force-free
- Magnetohydrostatic (MHS) models
- Magnetohydrodynamic (MHD) models (dependent on [MPI-AMRVAC](https://github.com/amrvac/amrvac)

3. Mangetic topology analysis
- Null point, spine, and fan structures
- Bald pathces and magnetic dips
- Quasi-separatrix layers (QSLs)

4. Magnetic helicity computation
- Finite volume method, vector potential, etc.
- Twist, writhe, linking numbers
- Helicity flux integration
