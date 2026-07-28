This folder contains the necessary files and notebook for RDF/CDF analysis of MD simulation outputs. It is currently in progress for general overview and to be updated for best fit in code examples folder for lab.

RDF stands for radial distribution function - a 2D form of structural analysis of specified atoms/molecules within a simulation. Specifically, this form of anlaysis takes the pairwise correlations between two atoms some distance, r, from each other, and outputs the probability of the atom located at that distance, g(r).

Files needed:
    Simulation output files (.h5)
    r coordinates (.npy) --> can generate these in notebook

Import modules noted in notebook. Refer to numpy and matplotlib documentation for additional work help
    