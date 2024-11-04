# MHforCvT
A Python script for calculating  C_v(T) curves from the output of canonical simulations across a range of Nt temperatures. This tool processes simulation data and computes thermodynamic properties, providing a detailed analysis of the heat capacity as a function of temperature.  

Based primarily on the method given in overview in:  Citation: R. Poteau, F. Spiegelmann, and P. Labastie, Isomerisation and phase transitions in small sodium clusters, Z. Phys. D 30, 57 (1994).

Depends on a set of outputs from mh.py file, which should be run first to generate:
  - mh_xxx.dat files - individual histograms per constant temperature
  - TotArr_yyyy.dat -- file containing the n_ij counts ... sized NumberOfTemperatures x ConfigurationEnergyBins
  - Gasum_yyyy.dat -- file contiaining the average energy and temperature per simulation temperature
    
Some parts are hard coded and should be modified - these are listed in the comments of the script header.
