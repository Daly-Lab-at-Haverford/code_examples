import numpy as np

def PES_CCH(orig_coor, steps, qmin, qmax, H, CH, CR, R):

    qm_N = len(orig_coor)
    pes = np.zeros((1,3))

    # Calculate the equilibrium bond lengths
    r0_r_cr  = 1.4540038408329599
    r0_cr_ch = 1.2081221370842437
    r0_ch_h  = 1.0657282305752254
    
    # Define universal bond vector displacements using propyne MP2/aug-cc-pVTZ
    v_r_cr = 0.562764142774114
    v_cr_ch = 0.9671529309507624
    v_ch_h = 0.31655937056996647
    
    ### Calculate the potential energy surface (PES) along the normal mode coordinate

    q = np.linspace(qmin, qmax, steps) #return 20 evenly spaced values ranging from -0.3 and +0.5

    # Generate 20 geometries along the normal mode coordinate grid space
    for n in q:
        orig_coor_q = np.copy(orig_coor)

        # calculate the unit bond vectors
        r_cr = orig_coor_q[CR] - orig_coor_q[R]
        cr_ch = orig_coor_q[CH] - orig_coor_q[CR]
        ch_h = orig_coor_q[H] - orig_coor_q[CH]

        r_cr_len = np.linalg.norm(r_cr)
        r_cr_dir = r_cr/r_cr_len

        cr_ch_len = np.linalg.norm(cr_ch)
        cr_ch_dir = cr_ch/cr_ch_len

        ch_h_len = np.linalg.norm(r_cr)
        ch_h_dir = ch_h/ch_h_len

        # mass-weight and average the unit bond vectors 
        new_dir = (12.*r_cr_dir + 12.*cr_ch_dir + 1.*ch_h_dir)/25.

        # move the atoms to scan q
        r = orig_coor_q[R]
        cr_new = new_dir*(n*(-1)*v_r_cr + r0_r_cr)
        cr_new += r

        ch_new = new_dir*(n*v_cr_ch + r0_cr_ch)
        ch_new += cr_new

        h_new = new_dir*(n*v_ch_h + r0_ch_h)
        h_new += ch_new

        orig_coor_q[H]  = h_new[0], h_new[1], h_new[2]
        orig_coor_q[CH] = ch_new[0], ch_new[1], ch_new[2]
        orig_coor_q[CR] = cr_new[0], cr_new[1], cr_new[2]

        pes = np.vstack((pes, orig_coor_q))
    
    pes = pes[1:]

    return pes, qm_N