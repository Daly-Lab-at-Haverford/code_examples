## QM + MM DVR Frequencies (from MD Snapshots)

## What is this code for?
* This code example can best be implemented on any workstation that uses a Slurm Workload Manager (Expanse, Palmetto)
* 
## What you need
**note** it is helpful to create a conda environment with all the necessary packages 
(add 'source activate [conda env]' to .job files)
* An MD production file of type .h5
* a forcefield field file of type .xml
* python v.3.8 or later
    * NumPy
    * pyQchem - Hauser Github
* QChem v.5.4
* mdTraj

## File structure
* There are two template folders that contains all of the base python scripts and slurm job files necessary
to calculate the DVR frequency for a terminal alkyne stretch (20-point PES)
    * **template_qm-mm** finds the frequency of a structure with solvents in the quantum region
    and a molecular mechanics region. The mm-region treats all other solvents as point charges
        * the qm region can be set to any value, including 0
    * **template_qm** finds the frequency of a structure with solvents in the quantum region
    and with NO molecular mechanics region
* The example h5 file here is from an MD simulation of propargyl acetate in triethylamine
and it contains 500 frames
* getQCInput.py is where you can customize the calculation to the specifications of your system.
    * DFT functional and basis set
    * probe and solvent type
    * center of the trajectory, specified by atom number
        * in this example it is the terminal hydrogen [12,12]
    * reference point for the probe, specified by atom number
        * in this example it is the terminal hydrogen [12]
    * number of solvent molecules in quantum region
        * in this example there are 2
    * the cutoff distance for the molecular region
        * in this example it is 2.8 nm
* Within the getQCInput.py file is the startingGeometry() function
    * This is where you can specify the number of gridpoints in your potential energy surface
        * here a 20-point grid from -0.3-0.5 is used
    * The atoms involved in the frequency stretch are also specified here by atom number
        * important note: atom numbers are 0-indexed here

## How to use
The implementation of files in template_qm-mm and template_qm are the same
1. Make all necessary updates to getQCInput.py (see File Structure section)
2. Use getQCInputs.job to create input files for each frame to calculate the single-point energy across a given PES
    a. '#SBATCH --array=0-1' this is where you specify the frames of interest, 
    to run a single frame set array = [frame number]
    b. sbatch getQCInputs.job, all SPEs will run simultaneously
3. Use runSPEs.job to calculate single-point energy for each point in each frame
4. Use getQCOutputs.job to run getQCOutput.py to extract final SCF energy of each point for each frame
    a. as above, update array as necessary
    b. sbatch getQCOutputs.job, to extract simultaneously
5. getQCOutputs.job also runs calcDVR.py to calculate the DVR frequency of each frame
and saves the frequencies to "snapshotDataPrecise.txt"
    a. "snapshotDataPrecise.txt" is formatted as (frame number, DVR frequency)

## Important Notes
* This code-example uses an implementation of the Universal DVR approach for terminal alkynes
see AlkyneUniversalDVR code-example for more information