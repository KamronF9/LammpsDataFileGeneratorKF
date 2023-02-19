#!/usr/bin/env python

# 2-2-2023 updated to convert output to both jdftx and outcar for dpmd or snn use

# module load python/3.11 openmpi lammps
# swap for "openmpi-gpu" if you run on a single core or on gpu3001 - like with deepmd 
# module load python/3.11 openmpi-gpu lammps

import numpy as np
# import sys
import gzip
import os
import glob

# if len(sys.argv)<2:
#     print('Usage: aim2LammpsLaunch.py <jdftx-log.out as md.o######>')
#     exit(1)
# titles = []
# os.chdir(r'/home/kamron/DFT/LQ_lammps/')
pattern = r'*gz'
i=0

# Cycle through cwd for all aimd outcar files (delta in this case)
for pathAndFilename in sorted(glob.iglob(os.path.join(os.getcwd(), pattern))):
    i+=1
    inFile = os.path.basename(pathAndFilename)
    title, ext = os.path.splitext(inFile)

    # titles.append(title[:-4])
    # os.chdir(r'/home/kamron/DFT/LQ_lammps/')
    # inFile = 'md.o105016'
    # baseOutFile = 'o105016'
    # inFile = sys.argv[1] #JDFTx format filename
#     baseOutFile = inFile[3:]
    baseOutFile = title
    outFileVasp = baseOutFile + '.vasp' # vasp poscar
    outFileJdftx = baseOutFile + '.jdftxout' # jdftxoutput

    # Dictionary for atoms
    atomDict = {'Na':1, 'Cl':2}

    # Units:
    eV = 1./27.2114  # div to go from Hartree to eV    27eV = 1eH
    Angstrom = 1/0.5291772  # div to go from Bohr to Ang.  0.5 Ang =1 Bohr

    # E_ClMinus = -15.102402370381217  # reference energy of Cl- ion
    # E_NaPlus = -47.435762646473520  # reference energy of Na+ ion

    # Read JDFTx input and convert to VASP
    nSteps = 0  # number of processed steps
    nEvery = 10  # select this many frames
    stepActive = False  # Whether to process current data
    latvecActive = False  # Whether reading lattice vectors
    stressActive = False  # Whether reading stress tensor
    atposActive = False  # Whether reading atomic positions
    forcesActive = False  # Whether reading forces
    headerWrittenoutcar = False
    headerWrittenjdftx = False
    rdfInited = False


    def writeLatticeOutcar(R, fp):
        print("      direct lattice vectors                 reciprocal lattice vectors", file=fp)
        latVecs = np.hstack((R.T, (2*np.pi)*np.linalg.inv(R)))  # direct and recip vectors
        for i in range(3):
            print('    {:12.9f} {:12.9f} {:12.9f}    {:12.9f} {:12.9f} {:12.9f}'.format(*tuple(latVecs[i])), file=fp)

    def writeLatticeJdftx(R, fp):
        print('R = ', file=fp)
        for i in range(3):
            print('[{:13.4f}{:13.4f}{:13.4f} ]'.format(*tuple(R[i])), file=fp)
                  
    fpoutOutcar = open(outFileVasp, "w")
    fpoutJdftx = open(outFileJdftx, "w")
    fpIn = gzip.open(inFile, 'rt') if (inFile[-3:]==".gz") else open(inFile,'r')
    for iLine, line in enumerate(fpIn):
        if line.find('total atoms') > 0:
            nAtoms = int(line.split()[4])
        if line.startswith('IonicDynamics: Step:'):
            tokens = line.split()
            iStep = int(tokens[2])
            stepActive = (iStep % nEvery == 0)
            PE = (float(tokens[4]) - (nAtoms/2)*(E_NaPlus+E_ClMinus))/eV  # how does this compare to the PE in LAMMPS.  I assume that is not adjusted either.  Is this just part of the vasp definition?  number of at pairs/total E ref.
            # total PE - ref energy for pairs
        # Lattice vectors:
        if latvecActive and iLine<refLine+3:
            iRow = iLine-refLine
            R[iRow] = [ float(tok)/Angstrom for tok in line.split()[1:-1] ]
            if iRow==2:
                latvecActive = False
        if line.startswith('R ='):
            latvecActive = True
            refLine = iLine+1
            R = np.zeros((3,3))
        # Stress tensor:
        if stressActive and iLine<refLine+3:
            iRow = iLine-refLine
            stress[iRow] = [ float(tok)/(eV/Angstrom**3) for tok in line.split()[1:-1] ]
            if iRow==2:
                stressActive = False
        if stepActive and line.startswith('# Stress tensor in'):
            stressActive = True
            refLine = iLine+1
            stress = np.zeros((3,3))
        # Atomic positions:
        if atposActive and iLine<refLine+nAtoms:
            iRow = iLine-refLine
            tokens = line.split()
            atNames.append(tokens[1])
            atNamesInd.append(atomDict[tokens[1]])
            atpos[iRow] = [ float(tok) for tok in tokens[2:5] ]
            if iRow+1==nAtoms:
                atposActive = False
                if coordsType == "cartesian":
                    atpos *= 1./Angstrom
                else:
                    atpos = np.dot(atpos, R.T) #convert to Cartesian (Angstrom)
                atNames = np.array(atNames)
        if stepActive and line.startswith('# Ionic positions in '):
            atposActive = True
            refLine = iLine+1
            atpos = np.zeros((nAtoms,3))  # units are XX
            atNames = []
            atNamesInd = []
            coordsType = line.split()[4]
        # Forces:
        if forcesActive and iLine < refLine + nAtoms:
            iRow = iLine-refLine
            tokens = line.split()
            forces[iRow] = [float(tok) for tok in tokens[2:5]]  
            if iRow+1==nAtoms:
                forcesActive = False
                if coordsType == "Cartesian":
                    forces *= 1./(eV/Angstrom)  # convert to eV/Ang from Hartree/Bohr
                else:
                    forces = np.dot(forces, np.linalg.inv(R)/eV)  # convert to Cartesian (eV/Angstrom)
        if stepActive and line.startswith('# Forces in '):
            forcesActive = True
            refLine = iLine+1
            forces = np.zeros((nAtoms,3))
            coordsType = line.split()[3]

        # Below marks when we complete a section of data and are ready to output data to file and run lammps
        if stepActive and line.startswith("# Energy components:"):
            # Write Lammps data file here (for each time step) with positions 
            # Save sample points to out file
            # First check that the structure is orthorhombic XX

            activeName = baseOutFile + '_' + str(nSteps)
            activeNameFile = activeName + r'.data'
            with open(activeNameFile,"w") as f:
                f.write("ITEM: TIMESTEP\n")
                f.write("  %s atoms\n"% (str(nAtoms)))
                f.write("  2 atom types\n")
                f.write("0 %s xlo xhi\n"% (str(R[0,0])))
                f.write("0 %s ylo yhi\n"% (str(R[1,1])))
                f.write("0 %s zlo zhi\n"% (str(R[2,2])))
                # add later tilts from lammps triclinic a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo).
                f.write("\n")
                f.write("Atoms\n")
                f.write("\n")
                for i in range(nAtoms):
                    f.write("  %s %s %s %s %s %s\n"% (str(i+1),
                                                      str(atNamesInd[i]),
                                                      str(0),
                                                      str(atpos[i][0]),
                                                      str(atpos[i][1]),
                                                      str(atpos[i][2])))  # index, type, 0(charge), x, y, z
                f.write("\n")

            # Run lammps on the single timestep to get Lammps to generate the forces and PE
            os.system("lmp < lammps.in -v struc %s -v base %s"% (str(activeNameFile), str(activeName)))
            # os.system("'/home/shankar/Software/LAMMPS/lammps-stable/build/lmp' < lammps.in -v struc %s -v base %s"% (str(activeNameFile), str(activeName)))

            # So now we can extract the force on each file here and sort the atoms based on id then store and compute diff against DFT

            def LammpForceExtract(inFile, nAtoms):
                # Go through Lammps dump file and gather the forces and sort based on ID
                refLine = 0
                atposActive = False  # Whether reading atomic positions
                for iLine, line in enumerate(open(inFile)):
                    # Extract ID and forces
                    if atposActive and iLine<refLine+nAtoms:
                        iRow = iLine-refLine
                        tokens = line.split()
                        atID.append(int(tokens[0]))  
                        forcesLmp[iRow] = [ float(tok) for tok in tokens[5:8] ]  # 5,6,7   units of force = eV/Angstrom
                        if iRow+1 == nAtoms:
                            atposActive = False
                    if line.startswith('ITEM: ATOMS'):
                        atposActive = True
                        refLine = iLine+1  # start of where to read in atom info
                        atID = []
                        forcesLmp = np.zeros((nAtoms,3))

                # Sort the forces based on an ordered atom index
                atID = np.expand_dims(atID, axis=1)  # expand dimension from 64, to 64,1 so to allow concatenation in the 2nd dimension
                IDForcesLmp = np.concatenate((atID, forcesLmp), axis=1)  # shape = nAtoms x 4
                IDForcesLmp = IDForcesLmp[IDForcesLmp[:, 0].argsort()]  # sort the forces in order by ID
                return IDForcesLmp[:, 1:]  # return the forces only portion no including the ID since it is sorted now

            def LammpEnergyExtract():
                # Go through Lammps log file and gather the PE of the system
                PEActive = False  # Whether to read PE
                # Read data for PE line
                for iLine, line in enumerate(open('log.lammps')):
                    if PEActive:
                        tokens = line.split()
                        PELammps = (float(tokens[1]))  # units of PE = eV
                        PEActive = False
                    # if line.startswith('Step PotEng Volume Press Temp'):
                    if line.startswith('   Step         PotEng'):
                        PEActive = True
                return PELammps

            # Call functions with the new file names
            inFile = activeName + '.dump'
            forcesLmp = LammpForceExtract(inFile, nAtoms)
            # Test line:
            # forcesLmp = LammpForceExtract('o105016_0.dump', nAtoms)

            # Take the differences b/w the lammps forces output and DFT, they are in the same units now then output into the format for outcar
            deltaForces = forces - forcesLmp  # eV/ang

            # Determine delta with PE.  PE in eV for both
            PELammps = LammpEnergyExtract()
            deltaPE = PE - PELammps  # eV

            # Clean up the written data and dump files after extraction
    #         os.chdir(r'/home/kamron/DFT/LQ_lammps/')
    #         os.remove(os.path.join(os.getcwd(),'o105016_0'+'.dump'))
            os.remove(os.path.join(os.getcwd(),activeName + '.dump'))
            os.remove(os.path.join(os.getcwd(),activeName + '.data'))
            os.remove(os.path.join(os.getcwd(),'log.lammps'))

            
            
            # Frame complete: write to OUTCAR
            
            if(not headerWrittenoutcar):
                print('Writing header')
                for i in range(2):
                    atomCounts = ""
                    for atName in np.unique(atNames):
                        print(" POTCAR:    PAW_PBE ", atName, "31Dec1999", file=fpoutOutcar)
                        atomCounts += "  " + str(len(np.where(atNames == atName)[0]))
                print("   ions per type =            ", atomCounts, file=fpoutOutcar)
                writeLatticeOutcar(R, fpoutOutcar)
                headerWrittenoutcar = True
            print('Writing frame for time step:', iStep)
            print('--------------------------------------- Iteration{:7d}(   1)  ---------------------------------------'.format(iStep), file=fpoutOutcar)
            print('  FORCE on cell =-STRESS in cart. coord.  units (eV):', file=fpoutOutcar)
            print('  Direction    XX          YY          ZZ          XY          YZ          ZX', file=fpoutOutcar)
            print('  --------------------------------------------------------------------------------------', file=fpoutOutcar)
            print('  Alpha Z {:11.5f} {:11.5f} {:11.5f}'.format(0,0,0), file=fpoutOutcar)
            for component in ['Ewald  ', 'Hartree', 'E(xc)  ', 'Local  ', 'n-local', 'augment', 'Kinetic', 'Fock   ']:
                print('  {:s} {:11.5f} {:11.5f} {:11.5f} {:11.5f} {:11.5f} {:11.5f}'.format(component,0,0,0,0,0,0), file=fpoutOutcar)
            print('  -------------------------------------------------------------------------------------', file=fpoutOutcar)
            for component in ['Total  ', 'in kB  ']:
                print('  {:s} {:11.5f} {:11.5f} {:11.5f} {:11.5f} {:11.5f} {:11.5f}'.format(component,0,0,0,0,0,0), file=fpoutOutcar)
            print('  external pressure =        0.00 kB  Pullay stress =        0.00 kB', file=fpoutOutcar)
            writeLatticeOutcar(R, fpoutOutcar)
            print(' POSITION                                       TOTAL-FORCE (eV/Angst)', file=fpoutOutcar)
            print(' -----------------------------------------------------------------------------------', file=fpoutOutcar)
            for i in np.argsort(atNames):
                posForce = tuple(atpos[i]) + tuple(deltaForces[i])
                print(' {:12.5f} {:12.5f} {:12.5f}    {:13.6f} {:13.6f} {:13.6f}'.format(*posForce), file=fpoutOutcar)
            print(' -----------------------------------------------------------------------------------', file=fpoutOutcar)
            print('    total drift:                                0.000000      0.000000      0.000000\n\n', file=fpoutOutcar)
            print('--------------------------------------------------------------------------------------------------------', file=fpoutOutcar)
            print('  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)', file=fpoutOutcar)
            print('  ---------------------------------------------------', file=fpoutOutcar)
            print('  free  energy   TOTEN  =  {:17.8f} eV\n'.format(deltaPE), file=fpoutOutcar)
            print('  energy  without entropy= {:17.8f}  energy(sigma->0) = {:17.8f}'.format(deltaPE,deltaPE), file=fpoutOutcar)
            
            
# KEY jdftx output lines used in deepmd dpdata import
# Initialized 2 species with 64 total atoms.

# IonicDynamics: Step:   0  PE: -2007.310584  KE:   0.089181  T[K]:  298.000  P[Bar]:      nan  tMD[fs]:      0.00  t[s]:     79.67

# R = 
# [      17.1025            0            0  ]
# [            0      17.1025            0  ]
# [            0            0      17.1025  ]

# # Ionic positions in cartesian coordinates:
# ion Na  17.809059999999999  13.866210000000001  18.914269999999998 v  -0.000147115007877   0.000189629242139   0.000131818380166 1

# # Forces in Cartesian coordinates:
# force Na  -0.025408645767808  -0.004065536242905   0.016991626218423 1            
            
    
    
#  [ ]  need to undo conversions to stay in jdftx


            # Frame complete: write to JDFTXOUT

            # output was converted to cartesian 
            # convert to jdftx units
            # delta forces is ev/ang need Hartree/Bohr
            # atpos is Ang need Bohr
            # R *= Angstrom but not here
            atpos *= Angstrom 
            deltaForces *= (eV/Angstrom)  # mult to go away from unit and to eH/Bohr
    
            if(not headerWrittenjdftx):
                print('Writing jdftx header')
                print("Initialized ", len(np.unique(atNames)), " species with ", len(atNames)," total atoms.", file=fpoutJdftx)
                headerWrittenjdftx = True
            # print('PE: {:11.6}'.format(-2007.310584))
            writeLatticeJdftx(R*Angstrom, fpoutJdftx)
            print('Writing frame for time step:', iStep)
            print('IonicDynamics: Step:   {:d}  PE: {:11.5f}'.format(iStep, deltaPE*eV ), file=fpoutJdftx)  # need to convert from eV to Eh
            

            


            print('# Ionic positions in cartesian coordinates:', file=fpoutJdftx)
            for i in np.argsort(atNames):
                # atNamePosTuple = atNames[i] + tuple(atpos[i])
                # print('ion {:s}{:20.15f}{:20.15f}{:20.15f}'.format(*tuple(atpos[i])), file=fpoutJdftx)
                print('ion ', atNames[i],'{:20.15f}{:20.15f}{:20.15f}'.format(*tuple(atpos[i])), file=fpoutJdftx)

            print('# Forces in Cartesian coordinates:', file=fpoutJdftx)
            for i in np.argsort(atNames):
                print('force ', atNames[i],'{:20.15f}{:20.15f}{:20.15f}'.format(*tuple(deltaForces[i])), file=fpoutJdftx)
            print('     Etot =', file=fpoutJdftx)

            # Transition steps below to the next segment in DFT output file
            nSteps += 1
            stepActive = False

    fpoutOutcar.close()
    fpoutJdftx.close()

    #exit(1)
