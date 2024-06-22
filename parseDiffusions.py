#PM env deepmd3

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

allDiffs = [] # size = samples/docs x 3 (bottom (pt/o surface),mid,top (pt surf))

fnames = sorted(glob.glob('out*'))
for fname in fnames:
    # fname = 'out1'

    with open(fname,'r') as f:
        rangeDiffs = [] # single file sample range limited component diffusion values
        for iline, line in enumerate(f):
            # bottom, mid, top order of file
            
            if line.startswith('diffusivity_c_range_components'):
                botDiff = line.split()
                # print(line.split())
                # print(float(line.split()[1][1:]))
                xTemp = float(line.split()[1][1:])
                yTemp = float(line.split()[2])
                zTemp = float(line.split()[3][:-1])
            if line.startswith('total'):
                totalTemp =float(line.split()[1])
                rangeDiffs.append([xTemp, yTemp, zTemp, totalTemp])
                # print([xTemp,yTemp,zTemp])

            # if iline==:
            #     botDiff = line.split()
            #     # print(line.split())
            #     # print(float(line.split()[1][1:]))
            #     xTemp = float(line.split()[1][1:])
            #     yTemp = float(line.split()[2])
            #     zTemp = float(line.split()[3][:-1])
            #     rangeDiffs.append([xTemp,yTemp,zTemp])
            #     # print([xTemp,yTemp,zTemp])

    # print(rangeDiffs)
    allDiffs.append(rangeDiffs)
    
# print(np.array(allDiffs).shape) # samples, c range, [x,y,z,total]
allDiffNP = np.array(allDiffs)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
print(colors)

labels = ['Pt/O Surface','Bulk','Pt Surface']

for i in range(3):
    plt.plot(allDiffNP[:,i,3],color=colors[i],label=labels[i])
    plt.hlines(np.mean(allDiffNP[:,i,3]),0,len(allDiffNP)-1,color=colors[i],linestyle='--')
plt.legend()
plt.savefig('totalDiff.pdf')
plt.close()


for i in range(3):
    plt.plot(allDiffNP[:,i,2])
plt.legend(['Pt/O Surface','Bulk','Pt Surface'])
plt.savefig('zDiff.pdf')
