CALS10k Model:

Text files:
CALS10l.1b & CALS10k.2 are the data file

Executables:
fieldcoefficients : used to get the field coefficients
fieldpred: used to get the field predictions

Folder configurationsCALS10k:
configurationsCALS10k: *.f are used to compile the executables using gfortran

This folder contains the configurations for gfortran to create the executables for CALS10k model
NEED GFORTRAN to be working
Command window: gfortran fieldcoeffs10k.f -o field_coefficients