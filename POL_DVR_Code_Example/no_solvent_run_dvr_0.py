# This program will perform a DVR calculation for
# a given terminal alkyne containing molecule
# at a specified level of theory

# Import packages and code
import argparse
import DVR_3atom
import extract_dipoles

probe = 'TA16'

calc_method = 'TPSS'
calc_basis = '6-311++G**'
n_proc = 80
N = 4
H = 0
CH = 1
CR = 2
R = 3

# Run DVR Calculation
print('Beginning {}/{}_{} DVR Calculation'.format(probe, calc_method, calc_basis))
DVR_3atom.run_3atom_DVR(probe, calc_method, calc_basis, n_proc, N, H, CH, CR, R)
print('{}/{}_{} DVR Calculation Complete!'.format(probe, calc_method, calc_basis))

print('Extracting Dipoles...')
extract_dipoles.extract_dm(probe)
print('{}/{}_{} Dipole Moments Extracted!'.format(probe, calc_method, calc_basis))
