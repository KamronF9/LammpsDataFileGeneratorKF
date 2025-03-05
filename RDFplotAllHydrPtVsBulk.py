#!/usr/bin/env python
# RDF only deepmd3
import numpy as np
import sys
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import sys

folderPrefixes = [
    '../allHydrRDFswPt/',
    '../allBulkRDFs/'
]

plotFile = "RDFallHydrTogetherQuad.pdf"

fig, axs = plt.subplots(2, 2, figsize=(8, 9), dpi=300, sharex=False, sharey=True)

for folderPrefix in folderPrefixes:

for i, ax in enumerate(axs.flat): 
    plt.sca(ax)

    # sys.exit(0)

    # if i in [0]

    # fnames = sorted(glob.glob(f'*SS*RDF.csv'))
    
    hydrLevels = [9, 12, 15]
    # hydrLevels = ['12Pt', '12Bulk', '12OFF']
    # hydrLevels = [9, 12, '12OFF', 15]
    # hydrLevels = [9, 12, '12OFF', '6JINN', '12JINN'] # 15
    # hydrLevels = [12, '12OFF', '12JINN'] # 15
    allHydrData = []

    for atPairs in ['SS', 'FF']:

        for hydrLevel in hydrLevels:
            fname = folderPrefix + f'hydr_{str(hydrLevel)}_{atPairs}_avg_RDF.csv'
            print(fname)

            df = pd.read_csv(fname)
            # allHydrData.append(df.copy())
            
            # g = gaussian_filter1d(df['g'],3) 
            g = df['g']
            plt.plot(df['r'], g, label=f'$\lambda$={str(hydrLevel)}')
            
        # plt.plot(rMid, rdf)
        if atPairs == 'SS':
            plt.xlim(0, 15)
            plt.ylim(0, 2)
        else:
            plt.xlim(0, 8)
            plt.ylim(0, 5)
        # plt.ylim(0, 20)
        plt.xlabel('r [A]')
        plt.ylabel(f'$g_{{{atPairs}}}(r)$')
        # plt.legend(['MgNa', 'NaCl', 'MgCl', 'ClCl'])
        # plt.legend(['0-40ps', '40-80ps'])
        # plt.legend(['0000', '0001'])
        # plt.legend(labels)
        plt.legend()
        plt.savefig(plotFile, bbox_inches='tight')
        plt.close()

