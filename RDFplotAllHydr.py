#!/usr/bin/env python
# RDF only deepmd3
import numpy as np
import sys
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import gaussian_filter1d
import pandas as pd

# for atPairs in ['CC','SS','OO', 'FF']:
for atPairs in ['SS', 'FF']:
# for atPairs in ['SS']:

# atPairs = 'OO' 
# atPairs = 'SS' 
# atPairs = sys.argv[1]
# hydrLevel = sys.argv[2]

    # fnames = sorted(glob.glob(f'*SS*RDF.csv'))
    plotFile = atPairs+"RDFallHydr.pdf"
    hydrLevels = [9, 12, 15]
    # hydrLevels = ['12Pt', '12Bulk', '12OFF']
    # hydrLevels = [9, 12, '12OFF', 15]
    # hydrLevels = [9, 12, '12OFF', '6JINN', '12JINN'] # 15
    # hydrLevels = [12, '12OFF', '12JINN'] # 15
    allHydrData = []

    plt.figure(figsize=(3,2))

    for hydrLevel in hydrLevels:
        fname = f'hydr_{str(hydrLevel)}_{atPairs}_avg_RDF.csv'
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

sys.exit(1)

atPairs = 'OO' 
# atPairs = 'SS' 
# atPairs = sys.argv[1]
# hydrLevel = sys.argv[2]

# fnames = sorted(glob.glob(f'*SS*RDF.csv'))
plotFile = atPairs+"StructureAllHydr.pdf"
hydrLevels = [9, 12, 15]
allHydrData = []

plt.figure(figsize=(4,3))

for hydrLevel in hydrLevels:
    
    fname = f'hydr_{str(hydrLevel)}_OO_avg_StructureFactor.csv'
    print(fname)
    df = pd.read_csv(fname)
    # allHydrData.append(df.copy())
    
    # Intensity = gaussian_filter1d(df['Intensity'],3) 
    Intensity = df['Intensity']
    plt.plot(df['q'], Intensity, label=f'$\lambda$={str(hydrLevel)}')
    
# plt.plot(rMid, rdf)
plt.xlim(0, 4)
plt.ylim(0, 2000)
plt.xlabel('q($\\AA^{-1}$)')
plt.ylabel('Intensity')
# plt.legend(['MgNa', 'NaCl', 'MgCl', 'ClCl'])
# plt.legend(['0-40ps', '40-80ps'])
# plt.legend(['0000', '0001'])
# plt.legend(labels)
plt.legend()
plt.savefig(plotFile, bbox_inches='tight')
plt.close()


# sys.exit(1)
'''
allData = [] # compile each data
lenData = []

plotFile = atPairs+"avgAllRDF.pdf"
plt.figure(figsize=(4,3))
# for ifile, file in enumerate(sorted(glob.glob('*7.rdf.dat*'))):
# for ifile, file in enumerate(sorted(glob.glob(f'{atPairs}_*.rdf.dat*'))):
#     data=np.genfromtxt(file, usecols=(0,1,2,3), names=True)
#     labels = data.dtype.names[1:]
    
    allData.append(data['g'+atPairs])
    lenData.append(len(data['g'+atPairs]))

    plotEach = False
    if plotEach:
    # for label in labels:
        for label in ['g'+atPairs]:
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
plt.legend([atPairs])
plt.savefig(plotFile, bbox_inches='tight')
plt.close()



# save data:
df = pd.DataFrame({'r':allXeven,'g':allYmean})
df.to_csv(f'hydr_{hydrLevel}_{atPairs}_avg_RDF.csv')


# print('DONE')


# Strucure factor:
from scipy.fft import fft

r = data['r']
g_r = allYmean
rho = g_r[-1]

# Method 1
k = np.linspace(0.01, max(r), 1000)  # Example k values
S_k = np.zeros_like(k)
for i, k_val in enumerate(k):
    integrand = r * (g_r - 1) * np.sin(k_val * r) / (k_val * r)
    S_k[i] = 1 + 4 * np.pi * rho * np.trapz(integrand, r)

# # Method 2
# # Calculate the structure factor
# k = np.fft.fftfreq(len(r), r[1] - r[0])
# S_k = 1 + rho * np.fft.fft(g_r * r**2) * (r[1] - r[0])

# # Normalize
S_k = np.real(S_k) / len(r)

plt.plot(k, S_k)
plt.xlim(0., 1.)
# plt.ylim(0, 20)
plt.xlabel('q($\\AA^{-1}$)')
plt.ylabel('Intensity')
plt.savefig(f'{atPairs}_StructureFactor.pdf', bbox_inches='tight')
plt.close()

# save data:
df = pd.DataFrame({'q':k,'Intensity':S_k})
df.to_csv(f'hydr_{hydrLevel}_{atPairs}_avg_StructureFactor.csv')
'''