# Define inputs in the run script
from pyqchem import Structure, QchemInput
from pyqchem.qc_input import CustomSection
import numpy as np
import argparse

#Define necessary variables
FORCEFIELD_PATH = "TEA_updated.xml"
TRAJECTORY_PATH = "PAC_in_TEA_prod.h5"
functional = 'TPSS'
basis_set = '6-311++G**'
probe_resname = 'PAC'
solv_resname = 'TEA'
c_idx = [12, 12] #list of atom indices to center the traj, i.e. probe terminal hydrogen
ref_idx = [12] #probe reference point, i.e. probe terminal hydrogen
N = 2 #number of solvent molecules in qm region

parser = argparse.ArgumentParser()
parser.add_argument('frame', type=int)
parser.add_argument('point', type=int)
frame = parser.parse_args().frame
point = parser.parse_args().point

framePath = "frame" + str(frame) + "/"


def getQCInput(geometrySymbols, geom_pt, frame, functional, basis_set):

    molecule = Structure(coordinates=geom_pt,
                         symbols=geometrySymbols,
                         charge=0,
                         multiplicity=1)


    #qchem resource options
    qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange=functional,
                      basis=basis_set,
                      scf_convergence=10,
                      thresh=14,
                      max_scf_cycles=200,
                      mem_total=100000,
                      mem_static=2000,
                      extra_rem_keywords={"DFT_D":"D3",
                                          "XC_SMART_GRID":"True",
                                          "XC_GRID":"000099000590"})
    #Calculate energy at given point
    qc_input_text = qc_input.get_txt()

    return qc_input_text


def startingGeometry():
    import MDprune
    #extract qm and mm region files from trajctory at given snapshot index
    MDprune.mdprune(TRAJECTORY_PATH, framePath, frame, probe_resname, solv_resname, c_idx, ref_idx, N, mm_cutoff, FORCEFIELD_PATH)

    import PES_CCH
    steps = 20 # change to steps wanted
    qmin = -0.3
    qmax = +0.5
    R  = 4 #0-indexed atom numbers in qm .xyz file
    CR = 5
    CH = 6
    H  = 12

    g = open(framePath + "qmregion_frame"+str(frame)+".xyz", "r")
    geom = []
    for line in g:
        if len(line.split())!=4:
            pass
        else:
            lineString = '     ' .join([str(line.split()[0][0]), line.split()[1], line.split()[2], line.split()[3], '\n'])
            geom.append(lineString)
    
    geometrySymbols = []
    geometryCoordinates = []
    for line in geom:
        parts = line.split()
        geometrySymbols.append(parts[0])
        xyz = [float(parts[1]),float(parts[2]),float(parts[3])]
        geometryCoordinates.append(xyz)
    
    geometryCoordinates = np.array(geometryCoordinates)

    pes, qm_N = PES_CCH.PES_CCH(geometryCoordinates, steps, qmin, qmax, H, CH, CR, R) #refer to PES_CCH for variables
    #write frame to master file to keep track of progress
    frameData = str(frame) + "\n"
    frameFile = open("frameTracker.txt", "a")
    frameFile.write(frameData)
    frameFile.close()

    return geometrySymbols, pes, qm_N

try:
    points, energies = np.loadtxt(framePath + 'energyFile.txt', dtype=float, delimiter = ',', unpack = True)
    pointsList = list(np.uint32(points).tolist())
except:
    pointsList = []
if point not in pointsList:
    
    geometrySymbols, pes, qm_N = startingGeometry()

    geom_pt = pes[point*qm_N:point*qm_N+qm_N]

    #Call function to calculate single-point energy for a given point on the surface
    qc_input_text = getQCInput(geometrySymbols, geom_pt, frame, functional, basis_set)

    qcFile = open(framePath + "point{:.0f}QCInput.qcin".format(point), "w")
    qcFile.write(qc_input_text)
    qcFile.close()
