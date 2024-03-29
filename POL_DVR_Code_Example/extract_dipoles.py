import numpy as np


#Define code for extracting dipole moments from polarizability calculation output file 
def extract_dm(probe):   
    steps = 20
    qmin  = -0.3
    qmax  = 0.5
    q= np.linspace(qmin, qmax, steps) #return 20 evenly spaced #s btwn -.3 and .5

    dline = 0
    dm   = np.zeros((steps, 3))
    for i, n in enumerate(q):
        name = 'POLs/pes_q_{}_pol.txt'.format(round(n,2))
        f = open(name)
        lines = f.readlines()
        for j, line in enumerate(lines):
            if "Dipole Moment (Debye)" in line:
                dline = j 
        dipolex = lines[dline+1].split()[1]
        dipoley = lines[dline+1].split()[3]
        dipolez = lines[dline+1].split()[5]
        dipole = np.array([dipolex, dipoley, dipolez])
        dm[i] = dipole


    np.savetxt('{}_dm.txt'.format(probe), dm)