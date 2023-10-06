# Alkyne Universal Displacement DVR
This code takes a .xyz file containing a system with a terminal alkyne moeity, and calculates the anharmonic vibrational frequency for the CC stretch using the discrete variable representation (DVR) method. 

The code uses a defined set of localized, universal displacements for the CC stretch normal mode from a localized harmonic vibrational frequency calculation for gas-phase propyne with MP2/aug-cc-pVTZ level of theory. These 'universal displacements' are applied to any molecule containing a terminal alkyne moeity, and the rationalization for this approach is explained in Streu, K et al, 2023 [1].


# How to use the code
1. Copy the files in this directory to a new directory in your workspace
2. Check that you have the following packages installed:
    * numpy
        * The code was developed with version 1.15.1
    * scipy
        * The code was developed with version 1.1.0
3. Update n_proc in run_dvr.py to reflect the number of CPUs to allocate for the single point energy calculations. Typically 8 CPUs is appropriate for the squirtle/bulbasaur/charmander workstations or 20 CPUs for eevee/mew/Palmetto. 
    Note: Increasing beyond 20 CPUs may be appropriate for large systems or expensive calculations, and should be done on mew or Palmetto.
4. In your workspace nagivate to the copied directory and type the following python command into terminal:

python run_dvr.py > run_dvr-out.txt

5. The above code will save the terminal output to `run_dvr-out.txt`, and you will find the CC stretch DVR frequency printed at the bottom of this output file.
6. Other outputs include:
    * /SPEs and /geometries directories which contain the QChem single point energy output files and .xyz geometries for each normal mode coordinate grid point respectively.
    * .npy binary file of the numpy array for the normal mode coordinate grid energies.
    * _elvls.txt and _wfns.txt files containing the discrete energy eigenvalues and wavefunctions.
    * calculation_data.db and __pycache__ containing stored calculation data.


# How to edit the code for a new .xyz file
1. Update the following parameters in the run_dvr.py file:
    * N is the integer total number of atoms in the system that appear in the .xyz file.
    * H, CH, CR, and R are the 0-indexed integer atom numbers corresponding to the H, CH, CR, and R atoms. The H atom is the terminal hydrogen, CH is the carbon in the triple-bond adjacent the terminal hydrogen, CR is the carbon in the triple-bond adjacent the R atom, and R is the atom adjacent the -CCH moeity.
2. Remove the example propyne.xyz file and replace with the .xyz file of interest.
3. If applicable, update the DFT method/basis. 
    Note: by default the calculation includes the following QChem keywords in the $rem section of the input file:
        * DFT_D = D3
        * DFT_GRID = 
    To change/remove the dispersion or DFT grid settings, edit extra_rem_keywords in DVR_3atom.py. 


# How the code works
1. The -CCH atoms are moved in space along the R-CR, CR-CH, and CH-H bond vectors according to the defined universal displacements to generate 20 unique geometries along the normal mode coordinate grid.
2. The energy is calculated for each geometry to create the discrete potential energy surface (PES) along the normal mode coordinate grid.
3. The PES grid minimum energy is set to zero.
4. The vibrational quantum eigenvalues and eigenvectors (wavefuctions and energy levels) are calculated for the given PES.


# Theoretical Background
This code was developed to calculate accurate and efficient CC stretch vibrational frequencies for molecules with terminal alkyne moieties. A brief theoretical background on the code is provided here, and a more detailed version is available in Streu, K. et al, 2023 [1]. 

[@Kristy write a brief theoretical background]


# References
1. Streu K, Hunsberger S, Patel J, Wan X, Daly Jr. CA. Localized-Normal Mode DVR for Terminal Alkyne CC Stretch Vibrations. ChemRxiv. Cambridge: Cambridge Open Engage; 2023; This content is a preprint and has not been peer-reviewed. DOI: https://doi.org/10.26434/chemrxiv-2023-2w6xx.

[@Kristy add references]
