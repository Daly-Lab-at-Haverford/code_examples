# This program will perform a DVR calculation for
# a given terminal alkyne containing molecule 
# at a specified level of theory

# Import packages and code
import argparse
import DVR_3atom
import extract_dipoles

parser = argparse.ArgumentParser()
parser.add_argument('solvent', type=int)
solvent = parser.parse_args().solvent
probe = 'PPY_1butanol'

calc_method = 'B3LYP-D3'
calc_basis = 'def2-TZVPD'
n_proc = 8
N = 22
H = 2
CH = 1
CR = 0
R = 3

# Run DVR Calculation
print('Beginning {}/{}_{} DVR Calculation'.format(probe, calc_method, calc_basis))
DVR_3atom.run_3atom_DVR(probe, calc_method, calc_basis, n_proc, N, H, CH, CR, R, solvent)
print('{}/{}_{} DVR Calculation Complete!'.format(probe, calc_method, calc_basis))

print('Extracting Dipoles...')
extract_dipoles.extract_dm(probe)
print('{}/{}_{} Dipole Moments Extracted!'.format(probe, calc_method, calc_basis))
