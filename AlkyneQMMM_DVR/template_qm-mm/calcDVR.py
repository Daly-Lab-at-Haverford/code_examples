import numpy as np
import dvr_v3 as dvr
import argparse

def calcDVR(framePath):

    steps=20 # change to steps wanted
    redmass=6.28664 #change to reduced mass wanted
    qmin=-0.3
    qmax=+0.5
    q = np.linspace(qmin, qmax, steps)

    #In what way were the SPEs stored?
    energyFile = open(framePath + 'energyFile.txt', "r")
    line0 = energyFile.read().split("\n")[0]
    energyFile.close()
    parts = line0.split(",")

    if len(parts) == 2:
        points, energies = np.loadtxt(framePath + 'energyFile.txt', dtype=float, delimiter = ',', unpack = True)

        points = points[-20:]
        energies = energies[-20:]
        orderedEnergies = np.zeros(20)
        for i in range(len(points)):
            orderedEnergies[int(points[i])] = energies[i]
        energies = orderedEnergies
    else:
        energies = np.loadtxt(framePath + 'energyFile.txt', dtype=float, unpack = True)

    for energy in energies:
        if energy == 0:
            return None

    # compute eigenvalues and eigenvectors for normal mode PES
    els, wfs = dvr.dvr1D(q[1]-q[0], energies, redmass)

    # save eigenvalues and eigenvectors for normal mode PES
    np.savetxt(framePath + 'pac_els_frame' + str(frame) + '.txt', els)
    np.savetxt(framePath + 'pac_wfs_frame' + str(frame) + '.txt', wfs)

    freq = els[1]-els[0]

    return (els, wfs, freq)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('frame', type=int)
    frame = parser.parse_args().frame

    framePath = "frame" + str(frame) + "/"

    dVRResults = calcDVR(framePath)
    if dVRResults != None:
        els, wfs, freq = dVRResults

        #save frameid and frequency to file
        freqData = frame, float(freq)
        value = ', '.join(map(str, freqData))
        freqStr = value + "\n"
        freqFile = open("snapshotDataPrecise.txt", "a")
        freqFile.write(freqStr)
        freqFile.close()
