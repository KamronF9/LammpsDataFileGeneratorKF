#!/usr/bin/env python
# RDF only deepmd3
import numpy as np
import sys
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import gaussian_filter1d
import pandas as pd

# for atPairs in ['SPt','SS','OO']:
# for atPairs in ['SS']:

# # atPairs = 'OO' 
# # atPairs = 'SS' 
# # atPairs = sys.argv[1]
# # hydrLevel = sys.argv[2]

# fnames = sorted(glob.glob(f'*SS*RDF.csv'))
plotFile = "zProfileAllHydr.pdf"

hydrLevels = [9]
# hydrLevels = [9, 12, 15]
allHydrData = []
zMax = [63.2, 68.6, 72.6]
ihydr = int(sys.argv[1]) # to choose which zmax to use
print('ihydr set zmax to', ihydr)

# plt.figure(figsize=(4,3))
plt,ax1 = plt.subplots(figsize=(4,3))


for hydrLevel in hydrLevels:
    
    # fname = f'hydr_{str(hydrLevel)}_{atPairs}_avg_RDF.csv'
    # fnames = ['densAll.txt', 'densSulfurAll.txt']
    
    # for fname in fnames:
    fname = f'zpolymer.txt'
    print(fname)
    df = pd.read_csv(fname, skiprows=4, index_col=False,  names=['index', 'fraction','number','density'], delim_whitespace=True)
    #,sep='\s+')
    
    print(df)
    print(df.columns)
    print(df.iloc[2])
    # allHydrData.append(df.copy())
    
    # g = gaussian_filter1d(df['g'],3) 
    # g = df['g']
    color = 'tab:red'
    # ax1.set_xlabel('Z fractional position')
    ax1.set_xlabel('Z (Ang)')
    ax1.set_ylabel('Mass Density ($g/cm^3$)')
    # ax1.plot(df['fraction']*zMax[ihydr], gaussian_filter1d(df['density'],1), label=f'Polymer')
    ax1.plot(df['fraction']*zMax[ihydr], df['density'], label=f'Polymer')
    ax1.set_ylim(0.05, 10)
    ax1.set_yscale('log')
    

    
    fname = f'zwaterHydron.txt'
    print(fname)
    color = 'tab:blue'
    
    df = pd.read_csv(fname, skiprows=4, index_col=False,  names=['index', 'fraction','number','density'], delim_whitespace=True)
    # ax2 = ax1.twinx() 
    # ax2.set_ylabel('Water Mass Density ($g/cm^3$)', color=color)
    ax1.plot(df['fraction']*zMax[ihydr], df['density'], label=f'Water')
    # ax2.set_ylim(0.05, 1)
    # ax2.set_ylim(0.01, 10)
    # ax1.set_yscale('log')
    ax1.legend()

# plt.plot(rMid, rdf)
# plt.xlim(0, 15)
# plt.ylim(0, 20)
# plt.xlabel('r [A]')
# plt.ylabel('g(r)')
# plt.legend(['MgNa', 'NaCl', 'MgCl', 'ClCl'])
# plt.legend(['0-40ps', '40-80ps'])
# plt.legend(['0000', '0001'])
# plt.legend(labels)
# plt.legend()
plt.savefig(plotFile, bbox_inches='tight')
plt.close()

# sys.exit(1)
'''
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