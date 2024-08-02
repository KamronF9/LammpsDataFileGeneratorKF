# coding: utf-8
# Python 3.10.4
# Distributed under the terms of the MIT License.
# If you make an update that you feel others would need, feel free to make
# a merge request to the main branch with your update
# (https://github.com/cpueschel/Lammps-Data-File-Generator.git).

# dpdata1 local env
#numpy                     1.23.3                   pypi_0    pypi
#pymatgen                  2023.1.30                pypi_0    pypi

#PM env deepmd3

#TODO  add read in of structure from poscar instead of config

from __future__ import print_function
import numpy as np
from numpy import linalg as LA
import math
# from pymatgen.io.xyz import XYZ
# from pymatgen.core import structure
# from pymatgen.io.lammps import outputs
# from pymatgen.core.sites import PeriodicSite
# from pymatgen.core import Site
# from pymatgen.core import SETTINGS, Element, Lattice, Structure
# from pymatgen import Lattice, Structure, Molecule
import yaml
# from pymatgen.io.vasp import Poscar, sets
from string import digits
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import numpy as np

# [H_PtTotal,O_PtTotal,O_OTotal,O_HTotal] # append to df

df = pd.read_csv('data.csv')

print(df)

plt.figure(figsize=(7,5), dpi=300)

for i,col in enumerate(df.columns):
    # if i==1:continue
    data=gaussian_filter1d(df[col],2)
    # data = np.cumsum(df[col])
    plt.plot(data, label=col)

plt.xlabel('100fs interval')
plt.ylabel('Count')
plt.legend()
plt.grid(True)

plt.savefig('surfaceReactionCounts.png')
