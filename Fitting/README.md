# Overview

This is a guide on calculating and fitting potential energy surfaces with specific attention to obtaining CHARMM-compatable values to be used in OpenMM simulations.

# Required Packages

* psi4
    * Any electronic structure program can be used, but Psi4 is easy to install and use in an otherwise pythonic environment.
    * Some of the optimizations were performed using Q-Chem.
* jupyter
* notebook
* matplotlib
* numpy
* scipy

These can be installed using conda/pip:

`conda install -c conda-forge psi4 jupyter notebook matplotlib numpy scipy`

# Most Useful Links

* scipy.optimize.curve_fit documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
* CHARMM-GUI on CHARMM Force Field functions and parameter files: https://www.charmm-gui.org/?doc=lecture&module=molecules_and_topology&lesson=3

# Included Files

* Jupyter notebooks for fittng bonds, angles, torsions, and out-of-plane angles to CHARMM style potentials.
* The energies and angles for the PSi4 angle scan.
* An example Q-Chem PES scan.


# Notes

This folder contains code used to obtain CHARMM compatable force field parameters for azirimycin. Atomic point charges were separatly obtained using FFParam, and most other interactions were taken from CGenFF where the analogy penalties were low. Quantum chemistry calculations were performed in Q-Chem 5.4 for all interactions except angles. For these, Psi4 was used. The code for Psi4 and the input for an example Q-Chem calculation are included.


# Other Useful Links

* Psi4 documentation: https://psicode.org/psi4manual/master/index.html
* On forcefields: https://computecanada.github.io/molmodsim-md-theory-lesson-novice/01-Force_Fields_and_Interactions/index.html  
* On molecular mechanics: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4026342/

# Citation

1. Daly C, Seebald L, Wolk E. Employing Metadynamics to Predict the Membrane Partitioning of Carboxy-2H-Azirine Natural Products. ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-m8v9m  This content is a preprint and has not been peer-reviewed.
2. Vanommeslaeghe, K.; Hatcher, E.; Acharya, C.; Kundu, S.; Zhong, S.; Shim, J.; Darian, E.; Guvench, O.; Lopes, P.; Vorobyov, I.; Mackerell Jr., A. D. CHARMM General Force Field: A Force Field for Drug-like Molecules Compatible with the CHARMM All-Atom Additive Biological Force Fields. Journal of Computational Chemistry 2010, 31 (4), 671–690. https://doi.org/10.1002/jcc.21367.
3. Kumar, A.; Yoluk, O.; MacKerell Jr., A. D. FFParam: Standalone Package for CHARMM Additive and Drude Polarizable Force Field Parametrization of Small Molecules. Journal of Computational Chemistry 2020, 41 (9), 958–970. https://doi.org/10.1002/jcc.26138.
4. Evgeny Epifanovsky et al. Software for the frontiers of quantum chemistry: An overview of developments in the Q-Chem 5 package. [J. Chem. Phys.. 155, 084801 (2021)]
5. Psi4 1.4: Open-Source Software for High-Throughput Quantum Chemistry”, D. G. A. Smith et al. J. Chem. Phys. (2020). (doi: 10.1063/5.0006002).
