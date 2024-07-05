#!/usr/bin/env python
# RDF only deepmd3
import numpy as np
import sys
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import gaussian_filter1d



allData = [] # compile each data
lenData = []

plotFile = "AavgAllRDF.pdf"
plt.figure(figsize=(4,3))
# for ifile, file in enumerate(sorted(glob.glob('*7.rdf.dat*'))):
for ifile, file in enumerate(sorted(glob.glob('*.rdf.dat*'))):
    data=np.genfromtxt(file, usecols=(0,1,2,3), names=True)
    labels = data.dtype.names[1:]
    
    allData.append(data['gSS'])
    lenData.append(len(data['gSS']))

    plotEach = False
    if plotEach:
    # for label in labels:
        for label in ['gSS']:
            g = gaussian_filter1d(data[label],3) # was 5
            plt.plot(data['r'], g, label=f'Iter-{ifile}')
        # label=f'{ifile*10}-{(ifile+1)*10}ps'
        # label=f'{ifile*10}-{(ifile+1)*40}ps'
        # plt.plot(data['r'], smoothg,label=label)

plotAll = True
if plotAll:
    # print(lenData)
    sampleNum = len(allData)
    minLen = min(lenData)
    allXeven = data['r'][:minLen]
    allYeven = np.array([allData[i][:minLen] for i in range(sampleNum)])
    allYmean = gaussian_filter1d(np.mean(allYeven,axis = 0),3)
    allYstd = np.std(allYeven,axis = 0)/sampleNum
    plt.plot(allXeven,allYmean)
    plt.fill_between(allXeven,allYmean-allYstd,allYmean+allYstd,alpha=0.5,color='red')

# plt.plot(rMid, rdf)
# plt.xlim(0, 10)
# plt.ylim(0, 20)
plt.xlabel('r [A]')
plt.ylabel('g(r)')
# plt.legend(['MgNa', 'NaCl', 'MgCl', 'ClCl'])
# plt.legend(['0-40ps', '40-80ps'])
# plt.legend(['0000', '0001'])
# plt.legend(labels)
plt.legend(['S-S'])
plt.savefig(plotFile, bbox_inches='tight')
# plt.close()

# print('DONE')
