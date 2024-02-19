#!/usr/bin/env python
# RDF only deepmd3
import numpy as np
import sys
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import gaussian_filter1d


plotFile = "rdf.pdf"
for ifile, file in enumerate(sorted(glob.glob('*.rdf.dat*'))):
    data=np.genfromtxt(file, usecols=(0,1,2,3), names=True)
    labels = data.dtype.names[1:]
    for label in labels:
        g = gaussian_filter1d(data[label],3) # was 5
        plt.plot(data['r'], g)
    # label=f'{ifile*10}-{(ifile+1)*10}ps'
    # label=f'{ifile*10}-{(ifile+1)*40}ps'
    # plt.plot(data['r'], smoothg,label=label)

# plt.plot(rMid, rdf)
plt.xlim(0, 10)
plt.ylim(0, 20)
plt.xlabel('r [A]')
plt.ylabel('g(r)')
# plt.legend(['MgNa', 'NaCl', 'MgCl', 'ClCl'])
# plt.legend(['0-40ps', '40-80ps'])
plt.legend(labels)
# plt.legend()
plt.savefig(plotFile, bbox_inches='tight')
# plt.close()

# print('DONE')
