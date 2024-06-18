import mdtraj as md
import numpy as np
import xmltodict

def center_traj(traj, res_name, c_idx):
    '''
    Takes a trajectory and the atom indices 
    (at least two), and centers the system with respect 
    to the specified atoms.

    Parameters
    ----------
    traj     : trajectory object
    res_name : a three letter string corresponding to the residue 
               name for the reference molecule(s), (ex: 'PAC')
    c_idx : a three letter string corresponding to the residue 
               name for the reference molecule(s), (ex: 'PAC')

    Returns
    -------
    c_traj : a trajectory object with the reference molecule 
             centered in each frame
    '''

    top = traj.topology
    res_idx = [residue.index for residue in top.residues if residue.name==res_name]
    for j in res_idx:
        anchor_list = list(set(top.residue(j).atoms))

    c_list = []
    for k in c_idx:
        c_list.append(anchor_list[k])

    c_traj = traj.image_molecules(anchor_molecules = [set(c_list)])

    return c_traj


def ref_xyz(traj, frame, ref_idx):
    '''
    Takes a trajectory, frame number, and a list of atom indices used as
    the reference (ex: the terminal alkyne hydrogen).
    Returns the Cartesian coordinates for the center of mass of the reference
    atoms for the specified frame.

    Parameters
    -----------
    traj    : trajectory object
    frame   : an integer corresponding to the trajectory frame of 
              interest, (ex: 1000)
    ref_idx : a list of integers corresponding to the reference atom 
              indices of interest, (ex: [6])

    Returns
    --------
    ref_coord : a numpy array with 3 entries that represent the xyz coordinates 
                of the center of mass for the specified reference atoms 
                for the specified frame.
    '''

    top = traj.topology

    xyz = traj[frame].atom_slice(ref_idx).xyz

    mass = np.zeros(len(ref_idx))
    for j, k in enumerate(ref_idx):
        mass_j = [atom.element.mass for atom in top.atoms
                            if atom.index==k]
        mass[j] = mass_j[0]

    ref_coord = sum(xyz[0, j, :]*mass[j] for j in 
                    range(len(mass)))/np.sum(mass)

    return ref_coord


def mol_xyz(traj, frame, res_name):
    '''
    Takes a trajectory, frame number, and molecule residue name. 
    Returns the Cartesian coordinates of the center of mass (COM) 
    for all molecules in the system.


    Parameters
    -----------
    traj     : trajectory object
    frame    : an integer corresponding to the frame of 
               interest, (ex: 1000)
    res_name : a three letter string corresponding to the molecule 
               residue name, (ex: 'TEA')

    Returns
    --------
    mol_coord : a 2D numpy array of shape (N,4) with N entires 
                corresponding to the number of molecules, and
                4 entires for each molecule that represent the
                molecule residue index number followed by the 
                Cartesian coordinates for the COM of each molecule 
                for the specified frame.
    '''

    top = traj.topology

    # make a list of all molecule indices
    res_idx = [residue.index for residue in top.residues 
               if residue.name==res_name]

    # create an empty numpy array with 4 entries (res_idx and xyz coords)
    # for each molecule
    mol_coord = np.zeros((len(res_idx), 4))

    # fill numpy array with molecule residue id and COM xyz coords
    for n, k in enumerate(res_idx):
        xyz = traj[frame].atom_slice([atom.index for atom in top.atoms 
                                      if atom.residue.index==k]).xyz

        mass = np.array([atom.element.mass for atom in top.atoms 
                         if atom.residue.index==k])

        mol_coord[n, 0] = k

        mol_coord[n, 1:] = np.sum(np.array([xyz[0, j, :]*mass[j] 
                                            for j in range(len(mass))]), 
                                            axis=0)/np.sum(mass)
    return mol_coord

def mol_dist(ref_coord, mol_coord):
    '''
    Takes the Cartesian coordinates for a reference point and the 
    center of masses for a set of molecules, and returns the molecule 
    distances from the reference.

    Parameters
    ----------
    ref_coord : a numpy array with 3 entries, xyz coordinates of a reference point.
    mol_coord : a 2D numpy array of shape (N,4) with N entries corresponding 
                to the number of molecules, and 4 entires for each molecule 
                which represent the molecule residue index number followed by the 
                Cartesian coordinates for the COM of each molecule for the 
                specified frame.

    Returns
    -------
    mol_dist : a numpy array with 2 entries for each molecule:
                    1. molecule residue index number.
                    2. computed molecule distances from atom_xyz.
    '''
    
    mol_dist = np.zeros((len(mol_coord),2))

    for j in np.arange(len(mol_coord)):
        mol_dist[j] = mol_coord[j,0], np.linalg.norm(mol_coord[j,1:]-ref_coord)

    return mol_dist


def create_qm(traj, framePath, frame, mm_cutoff, N, res_name, mol_dist):
    '''
    Takes the xyz Cartesian coordinates of a reference point and identifies 
    the closest N specified number of solvent molecules that will be treated 
    quantum mechanically.
    Writes an .xyz file with the Cartesian coordinates (in Angstroms) for all atoms 
    in the quantum region:
        1. All atoms for reference alkyne molecule
        2. All atoms for N number of solvent molecules
    Returns a numpy array of the solvent residue indices included in the qm region.

    Parameters
    ----------
    traj     : trajectory object.
    frame    : an integer corresponding to the trajectory frame of 
               interest, (ex: 1000).
    N        : an integer number of solvent molecules to be included in the 
               QM region.
    res_name : a three letter string corresponding to the reference residue 
               name, i.e. the alkyne probe molecule (ex: 'PAC').
    mol_dist : a numpy array with 2 entries for each molecule:
                  1. molecule residue index number.
                  2. computed molecule distances from atom_xyz.
    
    Returns
    -------
    qm_N : a 1D numpy array containing the solvent molecule residue indices for all 
           molecules included in the qm region.
    '''

    top = traj.topology

    # sort the list of solvent molecules by distance from reference
    mol_dist_sort = np.argsort(mol_dist[:,1])

    # identify closest N solvent molecules
    qm_N = mol_dist[mol_dist_sort][:N].T[0]

    qm_atom_idx = [atom.index for atom in top.atoms
                   if atom.residue.name==res_name
                   or atom.residue.index in qm_N]

    traj[frame].atom_slice(qm_atom_idx).save_xyz(framePath + 'qmregion_frame' 
                                                 + str(frame) + '.xyz')

    return qm_N


def create_mm(traj, framePath, frame, res_name, qm_N, mol_dist, mm_cutoff, ff_file):
    '''
    Identifies the remaining solvent molecules that will be treated as point charges
    with molecular mechanics. Writes an .xyz file with the Cartesian coordinates 
    (in Angstroms) for all atoms in the molecular mechanics region, and a .txt file 
    with atom charge and Cartesian coordinates (in Angstroms) for each mm region atom.
    Returns the sliced trajectory for the mm region atoms.

    Note: this code does not include a mm cutoff, which could be added later.

    Parameters
    ----------
    traj      : trajectory object.
    frame     : an integer corresponding to the trajectory frame of 
                interest, (ex: 1000).
    res_name  : a three letter string corresponding to the reference residue 
                name, i.e. the alkyne probe molecule (ex: 'PAC').
    qm_N      : a 1D numpy array containing the solvent molecule residue indices 
                for all molecules included in the qm region.
    mol_dist  : a numpy array with 2 entries for each molecule:
                  1. molecule residue index number.
                  2. computed molecule distances from atom_xyz.
    mm_cutoff : the cutoff distance for the mm region in nanometers. All molecules
                outside of this cutoff will be removed from the mm region.
    ff_file   : a file path to the solvent molecule force field .xml file.

    Returns
    -------
    mm_region : a 2D numpy array of shape (N, 4) with N entries for each mm
                atom included in the mm region, and 4 entries for each atom which
                represent the atom charge followed by the atom Cartesian coordinates.
    '''

    top = traj.topology
    N = len(qm_N)

    # use mm cutoff to cut solvent molecules out of the mm region
    cut_mm = mol_dist[:,0][mol_dist[:,1] > mm_cutoff]

    mm_atoms = [atom.index for atom in top.atoms
                if atom.residue.index not in qm_N
                and atom.residue.name!=res_name
                and atom.residue.index not in cut_mm]
    mm_frame = traj[frame].atom_slice(mm_atoms)
    mm_top = mm_frame.topology
    mm_frame.save_xyz(framePath + 'mmregion_frame' 
                      + str(frame) + '.xyz')
    mm_xyz = mm_frame.xyz[0, :, :] #xyz coords in nanometers
    mm_atom_types = [atom.name for atom in mm_top.atoms]

    # load all data from xml
    with open(ff_file) as fd:
        ff = xmltodict.parse(fd.read())

    # parse ff for relevant charge data by atom type
    charge_dict = ff['ForceField']['NonbondedForce']['Atom']
    typetoname = ff['ForceField']['Residues']['Residue']['Atom']

    # create a numpy array with charge and xyz coords for all mm atoms
    mm_region = np.zeros((len(mm_atom_types),4))

    for j, n in enumerate(mm_atom_types):
        for m in charge_dict:
            for k in typetoname:
                if k['@name']==n and m['@type']==k['@type']:
                    mm_region[j] = [float(m['@charge']), 
                                    mm_xyz[j, 0]*10, 
                                    mm_xyz[j, 1]*10, 
                                    mm_xyz[j, 2]*10]

    np.savetxt(framePath + '/mmregion_frame' + str(frame) 
               + '.txt', mm_region, delimiter=' ')

    return mm_region


def mdprune(traj_file, framePath, frame, probe_resname, solv_resname, c_idx, ref_idx, N, mm_cutoff, ff_file):

    traj = md.load(traj_file)
    ctraj = center_traj(traj, probe_resname, c_idx)
    H_coord = ref_xyz(ctraj, frame, ref_idx)
    tea_coord = mol_xyz(ctraj, frame, solv_resname)
    tea_dist = mol_dist(H_coord, tea_coord)
    qm_N_residx = create_qm(ctraj, framePath, frame, mm_cutoff, N, probe_resname, tea_dist)

    create_mm(ctraj, framePath, frame, probe_resname, qm_N_residx, tea_dist, mm_cutoff, ff_file)