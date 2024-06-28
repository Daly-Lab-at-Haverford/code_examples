import numpy as np
import mdtraj as md

def atom_center_traj(trajectory, atom_num, inplace=False):
    """
    center the trajectory based on a specific atom. 
    This atom will be at the center of the box in all frames

    Parameters
    ----------
    trajectory : mdtraj trajectory object
        The trajectory that should be centered

    atom_num : integer
        The index of the atom of that should be in the center of the box

    inplace : bool
        If True, the input trajectory will be centered; 
        if False, the input trajectory will be left alone
        detault is False, so original trajectory will be left as is.

    Returns
    -------
    traj_c : mdtraj trajectory object
        Centered trajectory
    """

    top = trajectory.topology
    
    atoms = [
    [i for i in top.atoms if i.index == atom_num][0],
    [i for i in top.atoms if i.index == atom_num][0]
    ]

    traj_c = trajectory.image_molecules(anchor_molecules = [set(atoms)])

    return traj_c

def bin_vol(r, dr):
    """
    calculates the volume of a bin in a CDF calculation

    Parameters
    ----------
    r : float or numpy array
        r value(s) in CDF calculation

    dr : float
        binsize. assumes dz=dr. 
        
    Returns
    -------
    bin_volume : array or single float of bin volume
    """
    return np.pi*(2*r*dr**2 + dr**3)

def atom_CDF(trajectory, center_atom, direction_atom, binning_atoms, n_bins=100):
    """
    Compute a cylindrical distribution function for a particular set of atoms
    
    Parameters
    ----------
    trajectory : mdtraj trajectory object
        The trajectory that the CDF should be calculated for

    center_atom : list containing one integer
        The index of the atom of that should be at r=0 and z=0

    direction_atom : list containing one integer
        The index of the atom of that should be at r=0 and along the z-axis

    binning_atoms : list of integers
        The indecies of the atoms which should be binned in the CDF

    n_bins : integer
        The number of bins in each direction. Default is 100, is usually sufficient

    Returns
    -------
    z_vals : array of z-axis values
    r_vals : array of r-axis values
    rho : g(r, z), the cylindrical distribution function
    """

    #center trajectory
    traj_c = atom_center_traj(trajectory, center_atom[0])

    #making sure n_bins is even so there aren't r=0 bins
    if n_bins % 2 != 0:
        n_bins += 1
        
    #defining useful variables, like the overall density and number of atoms
    n_atoms = len(binning_atoms)
    ndens = n_atoms/np.mean(traj_c.unitcell_volumes)
    #note that atoms can be outside this, mdtraj's image_molecules only requires the molecule's COMs to be inside the box.
    max_vec = np.linalg.norm(traj_c.unitcell_lengths[0])/2

    #binning settings, distance units are nm    
    rmin = -max_vec-0.5
    rmax = max_vec+0.5
    dr = (rmax-rmin)/(n_bins-1)
    zmin = -max_vec-0.5
    zmax = max_vec+0.5
    dz = (zmax-zmin)/(n_bins-1)

    frames = traj_c.n_frames
    rho = np.zeros((n_bins, n_bins))
    
    z_unit = np.zeros((n_atoms, 3))
    
    for frame in np.arange(frames):
    
        # get xyz coordinates for atoms of interest
        center_point = traj_c.xyz[frame, center_atom]
        dir_vec = traj_c.xyz[frame, direction_atom]-center_point
        binning_atom_vecs = traj_c.xyz[frame, binning_atoms]-center_point
    
        z_unit[:] = (dir_vec)[0]/np.linalg.norm(dir_vec)
    
        z = np.sum((binning_atom_vecs*z_unit), axis=-1)
        r = np.linalg.norm(binning_atom_vecs - np.expand_dims(z, 1)*z_unit, axis=-1)

        rho[np.round((z-zmin)/dz).astype(int), np.round((r-rmin)/dr).astype(int)] += 1/frames/ndens/bin_vol(r, dr) 
        rho[np.round((z-zmin)/dz).astype(int), np.round((rmax-r)/dr).astype(int)] += 1/frames/ndens/bin_vol(r, dr) 

        #make arrays for plotting, divide double counted r=0 line.
        z_array = np.linspace(zmin, zmax, n_bins)
        r_array = np.linspace(rmin, rmax, n_bins)
        #rho[:, r_array == 0] /= 2
        
    return (z_array, r_array, rho)
    