# MHforCvT
A Python script for calculating  C_v(T) curves from the output of canonical simulations across a range of Nt temperatures. This tool processes simulation data and computes thermodynamic properties, providing a detailed analysis of the heat capacity as a function of temperature.  

Based primarily on the method given in overview in:  Citation: R. Poteau, F. Spiegelmann, and P. Labastie, Isomerisation and phase transitions in small sodium clusters, Z. Phys. D 30, 57 (1994).

1)  Run mh.py first -- which assumes data is contained in VASP output files.  However, could pretty easily be modified for any other formats.

2)  canonCvtMain.py depends on a set of outputs from that mh.py execution:
  - mh_xxx.dat files - individual histograms per constant temperature
  - TotArr_yyyy.dat -- file containing the n_ij counts ... sized NumberOfTemperatures x ConfigurationEnergyBins
  - Gasum_yyyy.dat -- file contiaining the average energy and temperature per simulation temperature
    
Some parts of both scripts are hard coded and should be modified - these are listed in the comments of the script header.
