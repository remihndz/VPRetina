# coding: utf-8
import os 
import sys

## TODO: Update this so that it moves the CRA and CRV to the combined root

try:
    arteryFile = sys.argv[1] 
    veinFile = sys.argv[2]
    outputFile = sys.argv[3]
except IndexError:
    print("You did not specify two files to marge. Use 'python MergingArteryVeins.py arteryFile.cco veinFile.cco output.cco")
    sys.exit(1)

# Stores the files' information, line by line.
with open(arteryFile, 'r') as f:
    artery = f.readlines()
with open(veinFile, 'r') as f:
    veins = f.readlines()

# Number of segments in each file
nArteries = int(artery[4])
nVeins = int(veins[4])

# Find location of the CRA/CRV and create an artificial root to link them together for further CCO uses
rootVein = tuple(float(x) for x in veins[1].split(' ')[:3])
rootArtery = tuple(float(x) for x in artery[1].split(' ')[:3])
rootCombined = tuple((a+v)/2. for a,v in zip(rootArtery[:-1], rootVein[:-1]))+(-0.2,)

# The combined tree's info (one line)
newTreeInfo = ("*Tree\n",)+tuple(str(x)+' ' for x in rootCombined) + tuple(str(float(a)+float(v))+' ' for a,v in zip(artery[1].split(' ')[3:-1], veins[1].split(' ')[3:-1]))+("\n\n",)

# The combined tree's vessels info (nArteries+nVeins+3 segments) 
newVesselsInfo = ("*Vessels\n", str(nArteries+nVeins+3)+'\n') + tuple(a for a in artery[5:5+nArteries]) + tuple(' '.join((str(int(v.split(' ')[0])+nArteries), *v.split(' ')[1:]))  for v in veins[5:5+nVeins])

# Add the vessels linking the artificial root to the CRA and CRV
root = ' '.join((str(x) for x in [nArteries+nVeins, *rootCombined, *rootCombined[:-1], rootCombined[-1]+0.05, 0.0, 0.0, 0.0, 0.0, 0, newTreeInfo[-2], newTreeInfo[3], 0.0, 0.0, 0.0, 1, 0.0, 0.0, -2])) + '\n'
rootToCRA = ' '.join((str(x) for x in [nArteries+nVeins+1, *rootCombined[:-1], rootCombined[-1]+0.05, *rootArtery[:3], *artery[5].split(' ')[7:]]))
rootToCRV = ' '.join((str(x) for x in [nArteries+nVeins+2, *rootCombined[:-1], rootCombined[-1]+0.05, *rootVein[:3], *veins[5].split(' ')[7:]]))
newVesselsInfo += (root, rootToCRA, rootToCRV, '\n')

# The combined tree's connectivity information
connectivityArtery = artery[-nArteries-1:]
CRAConnectivity = connectivityArtery[1].split(' ')
connectivityArtery[1] = f"{CRAConnectivity[0]} {nArteries+nVeins+1} {' '.join(CRAConnectivity[2:])}" # Update the CRA connectivity: it is no longer the tree's root

connectivityVein = list(' '.join(str(int(n) + nArteries) for n in v.split(' '))+'\n' for v in veins[-nVeins:])
CRVConnectivity = connectivityVein[0].split(' ')
connectivityVein[0] = f"{CRVConnectivity[0]} {nArteries+nVeins+2} {' '.join(CRVConnectivity[2:])}" # Update the CRV connectivity: it is no longer the tree's root

newConnectivityInfo = tuple(connectivityArtery) + ('\n',) + tuple(connectivityVein) + (f"{nArteries+nVeins} -1 {nArteries+nVeins+1} {nArteries+nVeins+2}\n", f"{nArteries+nVeins+1} {nArteries+nVeins} 0\n", f"{nArteries+nVeins+2} {nArteries+nVeins} {nArteries}\n")

with open(outputFile, 'w') as f:
    f.write(''.join(str(x) for x in newTreeInfo+newVesselsInfo+newConnectivityInfo))

print("Merged tree written in", outputFile)
