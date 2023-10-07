# This program will perform a DVR calculation for
# a given terminal alkyne containing molecule 
# at a specified level of theory

# Import packages and code
import DVR_3atom

probe = 'PPY'
calc_method = 'TPSS'
calc_basis = '6-311++G**'
n_proc = 20
N = 7
H = 2
CH = 1
CR = 0
R = 3

# Run DVR Calculation
print('Beginning {}/{}_{} DVR Calculation'.format(probe, calc_method, calc_basis))
DVR_3atom.run_3atom_DVR(probe, calc_method, calc_basis, n_proc, N, H, CH, CR, R)
print('{}/{}_{} DVR Calculation Complete!'.format(probe, calc_method, calc_basis))

