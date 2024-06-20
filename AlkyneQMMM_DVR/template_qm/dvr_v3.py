import numpy as np
from scipy import linalg

def T(i, j, dx, redMass):
    a = ((-1.0)**(i-j)/(2.0*redMass*dx**2))
    if i == j:
        return a*np.pi**2/3.0
    else:
        return a*2.0/(i-j)**2

def d(i, j):
    if i == j:
        return 1.0
    else:
        return 0.0

def p(i, j, dx):
    b = ((-1.0)**(i-j))/(dx)
    if i == j:
        return b*0
    else:
        return b*1/(i-j)

def dvr1D(dx, energies, redMass):
    """
    Calculate vibrational
    quantum eigenvalues and eigenvectors
    (wavefuctions and energy levels)
    given a bonding potential energy surface.
    This is performed by following Eqn 2.1 of
    https://doi.org/10.1063/1.462100

    Parameters
    ----------
    distances : numpy array
        Sampled bond lengths of potential energy surface.
        The points should be equidistant, ordered,
        and in units of angstroms.

    energies : numpy array
        Sampled energies of potential energy surface.
        The points should be in Hartree.

    redMass : float
        Reduced mass of the bond or vibrational motion
        of interest in atomic mass units.

    Returns
    -------
    eigenvalues : numpy array
        ordered energy levels of potential energy surface.

    eigenvectors : numpy array
        ordered matrix of wavefunctions
    """

    #conversion units
    hart2waveNumber = 2.195*10**5
    amu2eMass = 1822.888
    ang2bohr = 1.88973

    #conversions to hartree atomic units
    redMass = redMass*amu2eMass
    dx = dx*ang2bohr
    nGrid = int(len(energies))

    #fill DVR matrix
    hamGrid = np.zeros((nGrid, nGrid))

    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            hamGrid[i, j] = T(i, j, dx, redMass) + d(i,j)*energies[i]

    #find eigenvalues, eigenvectors
    els, wfs = linalg.eigh(hamGrid)
    els *= hart2waveNumber
    wfs = wfs.T
    #return output
    return els, wfs

def dvr2D(dx_i, dx_j, energies, redMass_i, redMass_j):
    """
    Calculate vibrational quantum eigenvalues and eigenvectors
    given a bonding potential energy surface from sampling
    bond lengths of 2 bonds. Following Eqn 2.1 of
    https://doi.org/10.1063/1.462100 and Eqn 1 of
    https://doi.org/10.1063/1.468567.

     # are there cases where this equation shouldn't be used??

    Parameters
    ----------
     dx_i : float
        PES stepsize for first bond.
        In units of Angstroms.

     dx_j : float
        PES stepsize for second bond.
        In units of Angstroms.

    energies : numpy array
        Sampled energies of potential energy surface. Starting with all
        energies where the first bond is in its first bond length and going
        through each bond length of the second bond.
        Then moving to all energies where the first bond is in its second
        bond length and going through each bond length of the second bond,
        and continuing in this pattern.
        The points should be in Hartree.

    redMass_i : float
        Reduced mass of the first bond.
        In atomic mass units.

    redMass_j : float
        Reduced mass of the second bond.
        In atomic mass units.

    Returns
    -------
    eigenvalues : numpy array
        ordered energy levels of potential energy surface.

    eigenvectors : numpy array
        ordered matrix of wavefunctions
    """

    #conversion units
    hart2waveNumber = 2.195*10**5
    amu2eMass = 1822.888
    ang2bohr = 1.88973

    #conversions to hartree atomic units
    redMass_i = redMass_i*amu2eMass
    dx_i = dx_i*ang2bohr
    redMass_j = redMass_j*amu2eMass
    dx_j = dx_j*ang2bohr

    #4D to 2D transformation
    nGrid = int(len(energies)**(1/2))
    index = np.zeros((nGrid, nGrid))

    n = 0
    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            index[i, j] = n
            n += 1

    #fill DVR matrix
    hamGrid = np.zeros((nGrid**2, nGrid**2))
    energies = np.reshape(energies, (nGrid, nGrid))
    #print(energies)

    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            for ip in np.arange(nGrid):
                for jp in np.arange(nGrid):
                    ix = int(index[i, j])
                    iy = int(index[ip, jp])
                    hamGrid[ix, iy] = (d(i,ip)*T(j, jp, dx_j, redMass_j) +
                                       d(j,jp)*T(i, ip, dx_i, redMass_i) +
                                       d(i,ip)*d(j,jp)*energies[i, j]) #may need to be flipped
    #print(hamGrid)

    #find eigenvalues, eigenvectors
    els, wfs = linalg.eigh(hamGrid)
    els *= hart2waveNumber
    wfs = wfs.T
    #return output
    return els, wfs

def dvr2D_angle(dx_i, dx_j, theta, energies, redMass_i, redMass_j, m_0):
    """
    Calculate vibrational quantum eigenvalues and eigenvectors
    given a bonding potential energy surface from sampling
    bond lengths of 2 bonds. Following Eqn 2.1 of
    https://doi.org/10.1063/1.462100 and Eqn 1 of
    https://doi.org/10.1063/1.468567.

     # are there cases where this equation shouldn't be used??

    Parameters
    ----------
     dx_x : float
        PES stepsize for first bond.
        In units of Angstroms.

     dx_y : float
        PES stepsize for second bond.
        In units of Angstroms.

    theta :
        Bond Angle in degrees.

    energies : numpy array
        Sampled energies of potential energy surface. Starting with all
        energies where the first bond is in its first bond length and going
        through each bond length of the second bond.
        Then moving to all energies where the first bond is in its second
        bond length and going through each bond length of the second bond,
        and continuing in this pattern.
        The points should be in Hartree.

    redMass_x : float
        Reduced mass of the first bond.
        In atomic mass units.

    redMass_y : float
        Reduced mass of the second bond.
        In atomic mass units.

    m_0 : float
        Mass of the middle atom.
        In atomic mass units.

    Returns
    -------
    eigenvalues : numpy array
        ordered energy levels of potential energy surface.

    eigenvectors : numpy array
        ordered matrix of wavefunctions
    """

    #conversion units
    hart2waveNumber = 2.195*10**5
    amu2eMass = 1822.888
    ang2bohr = 1.88973

    #conversions to hartree atomic units
    redMass_i = redMass_i*amu2eMass
    dx_i = dx_i*ang2bohr
    redMass_j = redMass_j*amu2eMass
    dx_j = dx_j*ang2bohr
    m_0 = m_0*amu2eMass

    #4D to 2D transformation
    nGrid = int(len(energies)**(1/2))
    index = np.zeros((nGrid, nGrid))

    n = 0
    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            index[i, j] = n
            n += 1

    #fill DVR matrix
    hamGrid = np.zeros((nGrid**2, nGrid**2))
    energies = np.reshape(energies, (nGrid, nGrid))
    #print(energies)

    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            for ip in np.arange(nGrid):
                for jp in np.arange(nGrid):
                    ix = int(index[i, j])
                    iy = int(index[ip, jp])
                    hamGrid[ix, iy] = (d(i,ip)*T(j, jp, dx_j, redMass_j) +
                                       d(j,jp)*T(i, ip, dx_i, redMass_i) +
                                        -p(i, ip, dx_i)*p(j, jp, dx_j)*np.cos(theta*np.pi/180.0)/(m_0) + #there is a i*i in the denominator
                                       d(i,ip)*d(j,jp)*energies[i, j]) #may need to be flipped
    #print(hamGrid)

    #find eigenvalues, eigenvectors
    els, wfs = linalg.eigh(hamGrid)
    els *= hart2waveNumber
    wfs = wfs.T
    #return output
    return els, wfs

def dvr3D(dx_x, dx_y, dx_z, energies, redMass_x, redMass_y, redMass_z):
    """
    Calculate vibrational quantum eigenvalues and eigenvectors
    given a bonding potential energy surface from sampling
    bond lengths of 3 bonds. Following Eqn 2.1 of
    https://doi.org/10.1063/1.462100 and Eqn 1 of
    https://doi.org/10.1063/1.468567.

     # are there cases where this equation shouldn't be used??

    Parameters
    ----------
     dx_x : float
        PES stepsize for first bond.
        In units of Angstroms.

     dx_y : float
        PES stepsize for second bond.
        In units of Angstroms.

    dx_z : float
       PES stepsize for third bond.
       In units of Angstroms.

    energies : numpy array
        Sampled energies of potential energy surface. Starting with all
        energies where the first bond is in its first bond length and going
        through each bond length of the second bond.
        Then going through all energies where the third bond is in its first
        bond length.
        The points should be in Hartree.

    redMass_x : float
        Reduced mass of the first bond.
        In atomic mass units.

    redMass_y : float
        Reduced mass of the second bond.
        In atomic mass units.

    redMass_z : float
        Reduced mass of the third bond.
        In atomic mass units.


    Returns
    -------
    eigenvalues : numpy array
        ordered energy levels of potential energy surface.

    eigenvectors : numpy array
        ordered matrix of wavefunctions
    """
    #conversion units
    hart2waveNumber = 2.195*10**5
    amu2eMass = 1822.888
    ang2bohr = 1.88973

    #conversions to hartree atomic units
    redMass_x = redMass_x*amu2eMass
    dx_x = dx_x*ang2bohr
    redMass_y = redMass_y*amu2eMass
    dx_y = dx_y*ang2bohr
    redMass_z = redMass_z*amu2eMass
    dx_z = dx_z*ang2bohr

    #6D to 2D transformation
    nGrid = round((len(energies))**(1/3))
    index = np.zeros((nGrid, nGrid, nGrid))

    n = 0
    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            for k in np.arange(nGrid):
                index[i, j, k] = n
                n += 1

    #fill DVR matrix
    hamGrid = np.zeros((nGrid**3, nGrid**3))
    energies = np.reshape(energies, (nGrid, nGrid, nGrid))

    for i in np.arange(nGrid):
        for j in np.arange(nGrid):
            for k in np.arange(nGrid):
                for ip in np.arange(nGrid):
                    for jp in np.arange(nGrid):
                        for kp in np.arange(nGrid):
                            ix = int(index[i, j, k])
                            iy = int(index[ip, jp, kp])
                            hamGrid[ix, iy] = (d(i,ip)*d(j, jp)*T(k, kp, dx_z, redMass_z) +
                                       d(i,ip)*d(k, kp)*T(j, jp, dx_y, redMass_y)+ d(j,jp)*d(k, kp)*T(i, ip, dx_x, redMass_x) +
                                       d(i,ip)*d(j,jp)*d(k,kp)*energies[i, j, k])

    #find eigenvalues, eigenvectors
    els, wfs = linalg.eigh(hamGrid)
    els *= hart2waveNumber
    wfs = wfs.T
    #return output
    return els, wfs
