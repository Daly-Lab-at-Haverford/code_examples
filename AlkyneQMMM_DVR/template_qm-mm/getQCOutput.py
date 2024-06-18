import argparse

parser = argparse.ArgumentParser()
parser.add_argument('frame', type=int)
parser.add_argument('point', type=int)
frame = parser.parse_args().frame
point = parser.parse_args().point

framePath = "frame" + str(frame) + "/"

try:
    points, energies = np.loadtxt(framePath + 'energyFile.txt', dtype=float, delimiter = ',', unpack = True)
    pointsList = list(np.uint32(points).tolist())
except:
    pointsList = []
if point not in pointsList:
    qcFile = open(framePath + "point{:.0f}QCOutput.qcout".format(point), "r")
    qc_output_text = qcFile.read()
    qcFile.close()

    enum = qc_output_text.find('Total energy in the final basis set')
    spe = float(qc_output_text[enum : enum+100].split()[8])

    #append energy to a txt file
    energyFile = open(framePath + "energyFile.txt", "a")
    energyLine = str(point) + "," + str(spe) + "\n"
    energyFile.write(energyLine)
    energyFile.close()
