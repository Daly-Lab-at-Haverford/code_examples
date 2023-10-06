from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.file_io import write_structure_to_xyz
from pyqchem.errors import OutputError
import glob
import numpy as np
import os
import dvr_v3 as dvr

# Define code for running 3 atom DVR calculation
def run_3atom_DVR(probe, calc_method, calc_basis, n_proc, N, H, CH, CR, R):

    # Collect the molecule geometry and atom symbols
    xyz_file = glob.glob('*.xyz')

    with open(xyz_file[0], 'r') as file:
        xyz = file.read()

    xyz = xyz.split()
    xyz = np.array(xyz)
    xyz = xyz[1:]
    xyz = np.reshape(xyz, (N,4))

    orig_coor = xyz[:,[1,2,3]]
    orig_coor = orig_coor.astype('float64')
    print('Optimized Coordinates for DVR = \n' + str(orig_coor))

    symbols = xyz[:,0]
    symbols = list(symbols)
    print('Atom symbols for DVR = \n' + str(symbols))

    # Calculate the equilibrium bond lengths
    r0_r_cr  = np.linalg.norm(orig_coor[R] - orig_coor[CR])
    r0_cr_ch = np.linalg.norm(orig_coor[CR] - orig_coor[CH])
    r0_ch_h  = np.linalg.norm(orig_coor[CH] - orig_coor[H])
    print('r0_r_cr = ' + str(r0_r_cr))
    print('r0_cr_ch = ' + str(r0_cr_ch))
    print('r0_ch_h = ' + str(r0_ch_h))

    # Define universal reduced mass using propyne MP2/aug-cc-pVTZ
    red_mass = 6.286641266460249
    
    # Define universal bond vector displacements using propyne MP2/aug-cc-pVTZ
    v_r_cr = 0.562764142774114
    v_cr_ch = 0.9671529309507624
    v_ch_h = 0.31655937056996647

    # Create folders for organizing results
    try:
        os.mkdir('geometries')
    except FileExistsError:
        pass
    
    try:    
        os.mkdir('SPEs')
    except FileExistsError:
        pass
    
    ### Calculate the potential energy surface (PES) along the normal mode coordinate
    steps = 20 # Specify the normal mode coordinate grid spacing
    qmin  = -0.3
    qmax  = 0.5
    q     = np.linspace(qmin, qmax, steps) #return 20 evenly spaced values ranging from -0.3 and +0.5

    e_pes = np.zeros(steps)

    print('Calculating 20 single point energies')

    # Generate 20 geometries along the normal mode coordinate grid space
    for i, n in enumerate(q):
        orig_coor_q = np.copy(orig_coor)

        r = orig_coor_q[R]

        cr = orig_coor_q[CR] - r

        cr_len = np.linalg.norm(cr)
        cr_dir = cr/cr_len
        cr_new = cr_dir*(n*(-1)*v_r_cr + r0_r_cr)
        cr_new += r


        ch = orig_coor_q[CH] - cr_new

        ch_len = np.linalg.norm(ch)
        ch_dir = ch/ch_len
        ch_new = ch_dir*(n*v_cr_ch + r0_cr_ch)
        ch_new += cr_new


        h = orig_coor_q[H] - ch_new

        h_len = np.linalg.norm(h)
        h_dir = h/h_len
        h_new = h_dir*(n*v_ch_h + r0_ch_h)
        h_new += ch_new

        orig_coor_q[H]  = h_new[0], h_new[1], h_new[2]
        orig_coor_q[CH] = ch_new[0], ch_new[1], ch_new[2]
        orig_coor_q[CR] = cr_new[0], cr_new[1], cr_new[2]

        # Set up pyqchem/QChem single point energy (SPE) calculation
        molecule = Structure(coordinates=orig_coor_q,
                             symbols=symbols,
                             charge=0,
                             multiplicity=1)

        write_structure_to_xyz(molecule, 'geometries/' + 'geom_q_{}.xyz'.format(round(n,2)))

        qc_input = QchemInput(molecule,
                              jobtype='sp',
                              method=calc_method,
                              basis=calc_basis,
                              scf_convergence=9,
                              thresh=14,
                              max_scf_cycles=200,
                              scf_algorithm='rca_diis',
                              mem_total=100000,
                              mem_static=2000,
                              symmetry=False,
                              sym_ignore=True,
                              extra_rem_keywords={'DFT_D': 'D3',
                                                  'XC_SMART_GRID': True,
                                                  'XC_GRID': '000099000590',
                                                    })

        # Run pyqchem/QChem SPE calculation
        try:
            output = get_output_from_qchem(qc_input,
                                           processors=n_proc)
            
            # Save pyqchem/QChem SPE calculation output file
            outfile = open('SPEs/' + 'pes_q_{}_spe.txt'.format(round(n,2)), 'a')
            outfile.write(output)
            outfile.close()

            # Parse QChem output file for energy and save in a numpy array
            enum = output.find('Total energy in the final basis set =')
            energy = float(output[enum:enum+200].split()[8])
            e_pes[i] = energy

        except OutputError as e:
            print('Something wrong happened!:\n', e.error_lines)
            print('Recovering usable data...')

            try:
                output = e.full_output
            
                outfile = open('SPEs/' + 'pes_q_{}_spe.txt'.format(round(n,2)), 'a')
                outfile.write(output)
                outfile.close()

                enum = output.find('Total energy in the final basis set =')
                energy = float(output[enum:enum+200].split()[8])
                e_pes[i] = energy

            except:
                pass
        print('Finished single point energy calculation for Q = ' + str(round(n,2)))

    print('The PES energy array = \n' + str(e_pes))

    # Save PES energy array to .npy binary file
    np.save('{}_pes.npy'.format(probe), e_pes)


    ### Calculate the DVR frequency
    e_pes_0 = e_pes - e_pes.min() # zero the PES energies

    elvls, wfns = dvr.dvr1D(q[1]-q[0], e_pes_0, red_mass)
    dvr_freq = elvls[1]-elvls[0]
    print('The calculated normal mode DVR frequency = ' + str(dvr_freq))

    # Save energy levels and wavefunctions to .txt files
    np.savetxt('{}_elvls.txt'.format(probe), elvls)
    np.savetxt('{}_wfns.txt'.format(probe), wfns)

