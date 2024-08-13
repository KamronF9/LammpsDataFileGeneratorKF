#!/usr/bin/env python
# RDF only
import numpy as np
import sys
import glob

from pymatgen.io.lammps import data, outputs
from pymatgen.core import SETTINGS, Element, Lattice, Structure
import numpy as np
# from pymatgen.analysis.diffusion import analyzer
# from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer
import glob
# >>> data.
# data.ATOMS_HEADERS      data.LammpsBox(         data.Molecule(          data.Structure(         data.clean_lines(       data.pd                 
# data.CLASS2_KEYWORDS    data.LammpsData(        data.Path(              data.SymmOp(            data.itertools          data.re                 
# data.CombinedData(      data.Lattice(           data.SECTION_HEADERS    data.Topology(          data.lattice_2_lmpbox(  data.warnings           
# data.Element(           data.MODULE_DIR         data.SECTION_KEYWORDS   data.YAML(              data.loadfn(            data.zopen(             
# data.ForceField(        data.MSONable(          data.StringIO(          data.annotations        data.np                 


lammpsData = data.LammpsData.from_file('07polyHydroniumWaterPtExtHydr9Compressed.data')
# >>> lammpsData.
# lammpsData.REDIRECT                 lammpsData.from_dict(               lammpsData.masses                   lammpsData.unsafe_hash(
# lammpsData.as_dict(                 lammpsData.from_ff_and_topologies(  lammpsData.save(                    lammpsData.validate_monty_v1(
# lammpsData.atom_style               lammpsData.from_file(               lammpsData.set_charge_atom(         lammpsData.validate_monty_v2(
# lammpsData.atoms                    lammpsData.from_structure(          lammpsData.set_charge_atom_type(    lammpsData.velocities
# lammpsData.box                      lammpsData.get_partial_json(        lammpsData.structure                lammpsData.write_file(
# lammpsData.disassemble(             lammpsData.get_string(              lammpsData.to_json(                 
# lammpsData.force_field              lammpsData.load(                    lammpsData.topology 
# print(np.array(lammpsData.as_dict()['topology']['Bonds'])) 
# [[   2    1  128]
#  [   2    1    2]
#  [   3    1 1704]
#  ...
#  [   1 4686 5532]
#  [   1 4687 5530]
#  [   1 4690 5523]]
# print(np.array(lammpsData.as_dict()['topology']['Bonds'])[:,1]) 
# [   1    1    1 ... 4686 4687 4690]

atIDbonded1 = np.array(lammpsData.as_dict()['topology']['Bonds'])[:,1]
atIDbonded2 = np.array(lammpsData.as_dict()['topology']['Bonds'])[:,2]
atIDbondedAll = np.concatenate((atIDbonded1,atIDbonded2))
# print(np.shape(atIDbondedAll)) #(8664,)

# atID = 5000 # 1
# if atID not in atIDbondedAll:
#     print(atID, 'not bonded') # works 5000 is pt not bonded, 1 is.

# sys.exit(1)

# need at least two timesteps
# if len(sys.argv)<2:
#     print('Usage: RDF.py <lammps dump>')
#     exit(1)

# inFile = sys.argv[1] #lammps dump file
# RFile = sys.argv[2] #lattice dump file

for ifile, file in enumerate(sorted(glob.glob('test100fs*'))):
    inFile = file
    intervalNum = str(ifile) #inFile[-6] # interval int value 
    outFile = 'RDF'+ intervalNum #basename for RDF files


    # os.chdir(r'/home/kamron/NaCl_MgCl2/integrate')
    # inFile = '050MgCl2nnp_last.dump'
    # RFile = 'xx'
    # # RFile = 'lammps_R.out'
    # outFile = 'test'
        
    print('Make sure dump output order is -> ITEM: ATOMS id type x y z')


    # legLabel = 'Mg-Cl'  # label for the RDF plot later

    refLine = 0
    nSteps = 0 #number of processed steps
    # nEvery = 10 #select this many frames
    preprocess = True

    latvecActive = False #Whether reading lattice vectors
    # tricLat = False     #if latvec is orthogonal = False, triclinic =True
    stepActive = False #Whether to process current data - tailor later
    atomsActive = False #Whether to read total atoms
    atposActive = False #Whether reading atomic positions
    rdfInited = False #Whether RDF params have been initiated 


    saveIntermedRDFs = False
    # splitTimestepToRecalcRDF = 80000/4 # total steps in 4 chunks 0.5fs perstep is 10ps each
    # eachStep = 200 # 100fs at 0.5dt
    resetRDFnSteps = 100 # 10ps at 100fs per step
    saveRDF = False # if RDF should be saved yet
    saveNum = 0 # append this to data file for each interval

    # TO DO
    # delete first TIMESTEP and moved to end of file

    # # For testing lines
    # iLine = 0
    # f = open(inFile)
    # #----
    # line = f.readline()
    # iLine += 1 
    # line

    # Preprocess file - to make sure last line has TIMESTEP keyword to know when to stop
    if preprocess:
        import subprocess
        def tail(f, n, offset=0):
            proc = subprocess.Popen(['tail', '-n', str(n), f], stdout=subprocess.PIPE)
            lines = proc.stdout.readlines()
            return lines #[:, -offset]
        if not 'TIMESTEP' in str(tail(inFile, 1, 0)):  
            with open(inFile, 'a') as f:  # append mode
                f.write('TIMESTEP')


    for iLine,line in enumerate(open(inFile)):
        


        #Lattice vectors:
        if latvecActive and iLine<refLine+3:
            iRow = iLine-refLine
            # if tricLat and iRow == 2:
                
            #     for iRline,Rline in enumerate(open(RFile)):
            #         if iRline == 2:
            #             Tric = [ float(tok) for tok in Rline.split()[1:] ]  # a b c alpha beta gamma
            #             for iname,name in enumerate(['Ta','Tb','Tc','Talpha','Tbeta','Tgamma']):
            #                 globals()[name] = Tric[iname]

            #             a = np.array([ Ta, 0., 0. ]) # ax,ay,az
            #             b = np.array([ Tb*np.cos(Tgamma * np.pi/180.), Tb*np.sin(Tgamma  * np.pi/180.), 0. ]) #bx = xy, by, bz
            #             cx = Tc*np.cos(Tbeta * np.pi/180.)
            #             cy = (Tb*Tc*np.cos(Talpha * np.pi/180.) - b[0]*cx)/b[1]
            #             cz = np.sqrt(Tc**2 - cx**2 - cy**2)
            #             c = np.array([cx, cy, cz])
            #             # cx also xz (Tc*np.cos(Tbeta)), cy = yz (dot(b,c)-bx*cx)/by
            #             # also see https://lammps.sandia.gov/doc/Howto_triclinic.html for formula

            #             R = np.vstack((a,b,c)).T 
            
            # if not tricLat:
            # read each line
            # ITEM: BOX BOUNDS xy xz yz
            # xlo_bound xhi_bound xy
            # ylo_bound yhi_bound xz
            # zlo_bound zhi_bound yz
            bounds[iRow] = [ float(tok) for tok in line.split() ]

            if iRow==2:
                latvecActive = False

                R = np.array([[bounds[0,1]-bounds[0,0], 0. , 0.],
                    [0., bounds[1,1]-bounds[1,0] , 0.],
                    [0., 0., bounds[2,1]-bounds[2,0] ]])

        if line.startswith('ITEM: BOX BOUNDS'):
            latvecActive = True
            # assuming ortho
            # if line.find('xy xz yz') > 0:
            #     tricLat = True
            # else:
            #     tricLat = False
            refLine = iLine+1
            R = np.zeros((3,3))
            Tric = np.zeros((6))
            bounds = np.zeros((3,3)) # use depending on the dump style of box
            # bounds = np.zeros((3,2)) 
        # Atomic positions
        if atposActive and iLine<refLine+nAtoms:
            iRow = iLine-refLine
            tokens = line.split()
            atID = int(tokens[0])
            
            # add atoms to RDF only if not bonded
            if atID not in atIDbondedAll:
                atNames.append(tokens[2])  # index for name str Pt S
                atpos.append([ float(tok) for tok in tokens[3:6] ])
            # if q in dump  then use 3:6 - 3,4,5 adjusted since now charge is included 
            if iRow+1==nAtoms:
                atposActive = False
                atNames = np.array(atNames)
                atpos = np.array(atpos)
        if line.startswith('ITEM: ATOMS'):
            atposActive = True
            refLine = iLine+1 # start of where to read in atom positions
            # atpos = np.zeros((nAtoms,3))
            atpos = []
            atNames = []
        # Number of atoms
        if atomsActive:
            nAtoms = int(line.split()[0])
            atomsActive = False
        if line.find('NUMBER OF ATOMS') > 0:
            atomsActive = True
        # Final processing
        if (line.find('TIMESTEP') > 0) and (iLine>5):  # once you get to the end/beginning of the next tally up the RDF
            # reset RDF etc
            if nSteps%resetRDFnSteps==0 and nSteps>0 and saveIntermedRDFs:
                saveNum += 1
                # rdfInited = False 
                saveRDF = True

            # RDF initialize
            if not rdfInited:
                rMax = 0.5 * np.mean(np.diag(R))
                dr = 0.01
                rBins = np.arange(0., rMax, dr)
                rBins[0] = 0.01*dr #ignore self
                rMid = 0.5*(rBins[:-1]+rBins[1:])
                binVol = (4*np.pi/3)*(rBins[1:]**3 - rBins[:-1]**3)
                # print(binVol)
                numRDFs = 4
                rdf = np.zeros((len(rMid),numRDFs))
                rdfInited = True

            x = np.dot(atpos, np.linalg.inv(R.T))   # normalize positions to lattice shape
            # xS = x[np.where(atNames=='S')[0]]
            # xF = x[np.where(atNames=='F')[0]]
            xO = x[np.where(atNames=='O')[0]]
            # xNa = x[np.where(atNames==3)[0]]
            def getRDF(x1, x2):
                dx = x1[None,:,:] - x2[:,None,:]  # None adds a dimension 
                dx -= np.floor(0.5+dx) #minimum image convention
                r = np.linalg.norm(np.dot(dx, R.T), axis=-1).flatten()  
                # maybe done to cast relative coords onto coord basis
                # norm -1 takes -> min(sum(abs(x), axis=0))
                return np.histogram(r, rBins)[0] * (np.linalg.det(R) / (binVol * len(x1) * len(x2))) # local / bulk density
            rdf[:,0] += getRDF(xO, xO)
            # rdf[:,0] += getRDF(xS, xS)
            # rdf[:,1] += getRDF(xF, xF)
            # rdf[:,2] += getRDF(xO, xO)
            # rdf[:,2] += getRDF(xMg, xCl)
            # rdf[:,3] += getRDF(xCl, xCl)

            if saveRDF and saveIntermedRDFs:
                rdfFile = outFile+".rdf.dat"+str(saveNum)
                rdf *= (1./resetRDFnSteps)
                np.savetxt(rdfFile, np.hstack((rMid[:,None], rdf)), header='r gOO', comments='') #  gFF gOO gSS
                rdfInited = False # reset rdf
                saveRDF = False # reset save flag

            nSteps += 1

    # save all
    if not saveIntermedRDFs:
        rdfFile = outFile+".rdf.datAll"
        rdf *= (1./nSteps)
        np.savetxt(rdfFile, np.hstack((rMid[:,None], rdf)), header='r gOO', comments='') #  gFF gOO
        # rdfInited = False # reset rdf
        # saveRDF = False # reset save flag        

    print('DONE')
