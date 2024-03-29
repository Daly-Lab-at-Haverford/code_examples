DVR Frequency and Isotropic Polarizability Calculations: Code Example


Purpose: 


Overview: 

The files provided contain code to calculate the DVR (Discrete Variable Representation) frequency and isotropic polarizability.  

    <solventName>.xyz: This file contains the coordinates for the atoms in a given system.  In this example, the two molecules present are propyne and water.  

    run_dvr.py: This file contains the variables to be modified for each calculation.  The items to be changed are as follows: 
    
    probe = '<probe name>'
    calc_method = '<method name>'
    n_proc = '<number of processors>'
    N = <number of atoms in the system>
    H = <index of the terminal hydrogen atom>
    CH = <index of the carbon atom bonded to the terminal H atom>
    CR = <index of the carbon atom bonded to the CH atom and to the R group>
    R = <index of the atom in the R group bonded to the CR atom>
    

    DVR_3atom.py: This file contains code for running a 3-atom DVR calculation.  Unless the user wishes to modify which or where outputs are saved or which calculations are performed, there is no need to modify this file.  
    The code first collects the molecular geometry and atom symbols from the xyz file and sets up arrays of coordinates.  Then, equilibrium bond lengths are calculated and the universal reduced mass and universal bond vector displacements (using propyne MP2/aug-cc-pVTZ) are defined.  Next, the potential energy surface along the normal mode coordinate is calculated, and 20 geometries along the normal mode are generated.  Each geometry is saved to the 'geometries' folder.  Then the code to run calculations starts: an SPE (single point energy) calculation is set up, and the polarizability and DVR frequency are calculated.  Outputs (energy levels and wavefunctions) are saved to .txt files.  

    dvr_v3.py:  In the DVR_3atom.py file, the dvr_v3.py file is called upon to calculate the energy levels and wavefunctions.  ("elvls, wfns = dvr.dvr1D(q[1]-q[0], e_pes_0, red_mass)")
    As with the DVR_3atom.py file, the dvr_v3.py file should not have to be modified by the user in order to successfully run a calculation.  More details about the code can be found in the file itself.  

    extract_dipoles.py: This is a file parser that looks through the 'POLs/pes_q_{}_pol.txt' output file (one is created for each of the 20 geometries along the normal mode), and writes the resulting dipole moment to a <probe name>_dm.txt file.  This file does not need to be modified by the user in order to successfully run a calculation.  
    
    
    
Running a Calculation: 
First, the user should create a system of molecules using IqMol.  This code is designed to calculate the frequency of the CC triple bond and the polarizability of an alkyne, so the system should include exactly one alkyne and potentially another solvent molecule; here, propyne and water are used as the two molecules.  
After uploading the xyz file for the system to the same folder in which these files are located, the user should follow the instructions above on how to modify the run_dvr.py file appropriately.  

To run the calculation on a workstation, navigate to the directory through the terminal.  Note that it is recommended to use a 'qchem' conda environment.  Enter the following lines to run the calculation: 
conda activate qchem
 python run_dvr.py > run_dvr-out.txt
 
 After the calculation is complete, verify that the run_dvr-out.txt file contains both the DVR frequency and the transition isotropic polarizability for the CC stretch of the alkyne present.  



    
References: 

Kristina Streu, Sara Hunsberger, Jeanette Patel, Xiang Wan, Clyde A. Daly; Development of a universal method for vibrational analysis of the terminal alkyne C≡C stretch. J. Chem. Phys. 21 February 2024; 160 (7): 074106. https://doi.org/10.1063/5.0185580

Evgeny Epifanovsky, Andrew T. B. Gilbert, Xintian Feng, Joonho Lee, Yuezhi Mao, Narbe Mardirossian, Pavel Pokhilko, Alec F. White, Marc P. Coons, Adrian L. Dempwolff, Zhengting Gan, Diptarka Hait, Paul R. Horn, Leif D. Jacobson, Ilya Kaliman, Jörg Kussmann, Adrian W. Lange, Ka Un Lao, Daniel S. Levine, Jie Liu, Simon C. McKenzie, Adrian F. Morrison, Kaushik D. Nanda, Felix Plasser, Dirk R. Rehn, Marta L. Vidal, Zhi-Qiang You, Ying Zhu, Bushra Alam, Benjamin J. Albrecht, Abdulrahman Aldossary, Ethan Alguire, Josefine H. Andersen, Vishikh Athavale, Dennis Barton, Khadiza Begam, Andrew Behn, Nicole Bellonzi, Yves A. Bernard, Eric J. Berquist, Hugh G. A. Burton, Abel Carreras, Kevin Carter-Fenk, Romit Chakraborty, Alan D. Chien, Kristina D. Closser, Vale Cofer-Shabica, Saswata Dasgupta, Marc de Wergifosse, Jia Deng, Michael Diedenhofen, Hainam Do, Sebastian Ehlert, Po-Tung Fang, Shervin Fatehi, Qingguo Feng, Triet Friedhoff, James Gayvert, Qinghui Ge, Gergely Gidofalvi, Matthew Goldey, Joe Gomes, Cristina E. González-Espinoza, Sahil Gulania, Anastasia O. Gunina, Magnus W. D. Hanson-Heine, Phillip H. P. Harbach, Andreas Hauser, Michael F. Herbst, Mario Hernández Vera, Manuel Hodecker, Zachary C. Holden, Shannon Houck, Xunkun Huang, Kerwin Hui, Bang C. Huynh, Maxim Ivanov, Ádám Jász, Hyunjun Ji, Hanjie Jiang, Benjamin Kaduk, Sven Kähler, Kirill Khistyaev, Jaehoon Kim, Gergely Kis, Phil Klunzinger, Zsuzsanna Koczor-Benda, Joong Hoon Koh, Dimitri Kosenkov, Laura Koulias, Tim Kowalczyk, Caroline M. Krauter, Karl Kue, Alexander Kunitsa, Thomas Kus, István Ladjánszki, Arie Landau, Keith V. Lawler, Daniel Lefrancois, Susi Lehtola, Run R. Li, Yi-Pei Li, Jiashu Liang, Marcus Liebenthal, Hung-Hsuan Lin, You-Sheng Lin, Fenglai Liu, Kuan-Yu Liu, Matthias Loipersberger, Arne Luenser, Aaditya Manjanath, Prashant Manohar, Erum Mansoor, Sam F. Manzer, Shan-Ping Mao, Aleksandr V. Marenich, Thomas Markovich, Stephen Mason, Simon A. Maurer, Peter F. McLaughlin, Maximilian F. S. J. Menger, Jan-Michael Mewes, Stefanie A. Mewes, Pierpaolo Morgante, J. Wayne Mullinax, Katherine J. Oosterbaan, Garrette Paran, Alexander C. Paul, Suranjan K. Paul, Fabijan Pavošević, Zheng Pei, Stefan Prager, Emil I. Proynov, Ádám Rák, Eloy Ramos-Cordoba, Bhaskar Rana, Alan E. Rask, Adam Rettig, Ryan M. Richard, Fazle Rob, Elliot Rossomme, Tarek Scheele, Maximilian Scheurer, Matthias Schneider, Nickolai Sergueev, Shaama M. Sharada, Wojciech Skomorowski, David W. Small, Christopher J. Stein, Yu-Chuan Su, Eric J. Sundstrom, Zhen Tao, Jonathan Thirman, Gábor J. Tornai, Takashi Tsuchimochi, Norm M. Tubman, Srimukh Prasad Veccham, Oleg Vydrov, Jan Wenzel, Jon Witte, Atsushi Yamada, Kun Yao, Sina Yeganeh, Shane R. Yost, Alexander Zech, Igor Ying Zhang, Xing Zhang, Yu Zhang, Dmitry Zuev, Alán Aspuru-Guzik, Alexis T. Bell, Nicholas A. Besley, Ksenia B. Bravaya, Bernard R. Brooks, David Casanova, Jeng-Da Chai, Sonia Coriani, Christopher J. Cramer, György Cserey, A. Eugene DePrince III, Robert A. DiStasio Jr., Andreas Dreuw, Barry D. Dunietz, Thomas R. Furlani, William A. Goddard III, Sharon Hammes-Schiffer, Teresa Head-Gordon, Warren J. Hehre, Chao-Ping Hsu, Thomas-C. Jagau, Yousung Jung, Andreas Klamt, Jing Kong, Daniel S. Lambrecht, WanZhen Liang, Nicholas J. Mayhall, C. William McCurdy, Jeffrey B. Neaton, Christian Ochsenfeld, John A. Parkhill, Roberto Peverati, Vitaly A. Rassolov, Yihan Shao, Lyudmila V. Slipchenko, Tim Stauch, Ryan P. Steele, Joseph E. Subotnik, Alex J. W. Thom, Alexandre Tkatchenko, Donald G. Truhlar, Troy Van Voorhis, Tomasz A. Wesolowski, K. Birgitta Whaley, H. Lee Woodcock III, Paul M. Zimmerman, Shirin Faraji, Peter M. W. Gill, Martin Head-Gordon, John M. Herbert, and Anna I. Krylov. Software for the frontiers of quantum chemistry: An overview of developments in the Q-Chem 5 package. [J. Chem. Phys.. 155, 084801 (2021)]
