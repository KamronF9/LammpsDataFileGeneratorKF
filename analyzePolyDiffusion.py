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
from pymatgen.io.xyz import XYZ
from pymatgen.core import structure
from pymatgen.io.lammps import outputs
from pymatgen.core.sites import PeriodicSite
from pymatgen.core import Site
from pymatgen.core import SETTINGS, Element, Lattice, Structure
# from pymatgen import Lattice, Structure, Molecule
import yaml
from pymatgen.io.vasp import Poscar, sets
from string import digits
import sys
import pandas as pd

from pymatgen.analysis.diffusion import analyzer

# Returns known atom types in a list

# os.chdir(r'agso3Ex')  # adjust later



def known_atom_types():
    types = []
    for each in config['types']:
        types.append(str(each))
    return types
# Returns known atom information in a list


def known_atom_types_general():
    types = known_atom_types()
    known_atom_types_general = []
    for every in types:
        charges = config['types'][every]['CHARGES']
        lj_parameters = (config['types'][every]['LJ_PARAMETERS'][
                         0], config['types'][every]['LJ_PARAMETERS'][1])
        mass = config['types'][every]['MASS']
        required_connectors = config['types'][every]['REQUIRED_CONNECTORS']
        known_atom_types_general.append(
            (every, charges, lj_parameters, mass, required_connectors))
    return known_atom_types_general

# Returns list_bond_types: list of known atom bonds(list of sorted tuples)


def known_bond_types():
    list_bond_types = []
    for every in range(0, len(config['bonds'])):
        i = config['bonds'][every]['i']
        j = config['bonds'][every]['j']
        energy = config['bonds'][every]['ENERGY']
        length = config['bonds'][every]['LENGTH']
        list_bond_types.append((i, j, energy, length))
    return list_bond_types

# Returns list_angle_types: list of known angles(list of sorted tuples)


def known_angle_types():
    list_angle_types = []
    for every in range(0, len(config['angles'])):
        i = config['angles'][every]['i']
        j = config['angles'][every]['j']
        k = config['angles'][every]['k']
        energy = config['angles'][every]['ENERGY']
        theta = config['angles'][every]['THETA']
        list_angle_types.append((i, j, k, energy, theta))
    return list_angle_types

# Returns known_torsion_types: list of known torsions(list of sorted tuples)


def known_torsion_types():
    known_torsion_types = []
    for every in range(0, len(config['propertorsions'])):
        i = config['propertorsions'][every]['i']
        j = config['propertorsions'][every]['j']
        k = config['propertorsions'][every]['k']
        l = config['propertorsions'][every]['l']
        energy = config['propertorsions'][every]['ENERGY']
        angle = config['propertorsions'][every]['ANGLE']
        multiplicity = config['propertorsions'][every]['MULTIPLICITY']
        known_torsion_types.append((i, j, k, l, energy, angle, multiplicity))
    return known_torsion_types

# Max possible bond length for the system with tolerance factor


def max_bond_length():
    max_bond_length = 0
    # print(known_bond_types())
    for each in known_bond_types():
        bond_length = each[3]
        if bond_length > max_bond_length:
            max_bond_length = bond_length
        # print(each, max_bond_length)
    return float(max_bond_length + max_bond_length *
                 config['bond_length_tolerance_factor'])

# Checks if smaller list is in another list


def is_slice_in_list(s, l):
    len_s = len(s)
    return any(s == l[i:len_s + i] for i in xrange(len(l) - len_s + 1))

# Does this list have digits


def hasNumbers(inputLIST):
    for each in inputLIST:
        if any(char.isdigit() for char in each):
            return any(char.isdigit() for char in each)
    return False

# Assigns the order for execution of atom types.


def type_assignment_execution_order():
    """
            Returns an array for the order of execution of the type assignment
            eg. C1 is assigned before C2
            Order:
            None types assigned first
            Those with just "X" are assigned 2nd, X is assumed to be equal to X1
            Those with "C2" assigned second, etc.
    """
    def position_order(types, order):
        position_order = []
        for each in order:
            for i in range(0, len(types)):
                if types[i] == each:
                    position_order.append(i)
        return position_order

    # intialize
    atom_types_general = known_atom_types_general()
    i, types, dependencies, type_execution_order = 0, [
        None] * len(atom_types_general), [None] * len(atom_types_general), []

    for each in atom_types_general:
        types[i] = each[0]
        dependencies[i] = each[4].values()
        i += 1

    # 1 Check for items without dependencies.
    switch = 0
    while len(type_execution_order) < len(types):
        if switch == 0:
            for i in range(0, len(types)):
                if 'None' in dependencies[i] and types[
                        i] not in type_execution_order:
                    type_execution_order.append(types[i])

        if switch == 1:
            for i in range(0, len(types)):
                if 'None' in dependencies[i]:
                    continue
                if is_slice_in_list(dependencies[i], type_execution_order) and types[
                        i] not in type_execution_order:
                    type_execution_order.append(types[i])
        if switch == 2:
            for i in range(0, len(types)):
                if 'None' in dependencies[i]:
                    continue
                if is_slice_in_list(dependencies[i], type_execution_order) and types[
                        i] not in type_execution_order:
                    continue
                if not hasNumbers(dependencies[i]) and types[
                        i] not in type_execution_order:
                    type_execution_order.append(types[i])

        if switch == 3:
            for i in range(0, len(types)):
                if types[i] not in type_execution_order:
                    type_execution_order.append(types[i])

        switch = switch + 1
        if switch > 10:
            print("Arguments for input not properly ordered")
            break

    # Allow user to verify order
    print_types = []
    for each in (position_order(types, type_execution_order)):
        print_types.append(types[each])
    print("Order of execution for type assignment:", print_types)
    return position_order(types, type_execution_order)

# Given a atom type (eg. C3) returns number


def atom_type_NUMBER(type_of_atom):
    i = 0
    for each in known_atom_types_general():
        if type_of_atom == each[0]:
            return i
        i += 1
    return None

# Determines if two sites are bonded together


def site_bonded(site1, site2, bond_length):
    """
    Requirements for a bond:
            Must have a bond length = required bond length +/- varriance
            Must be the correct atomic species_string.
            Then it is identified as a "None" type of bond. Since specific atomic types need to be determined based on bonds.
    """
    site1type = atom_sites[site1].type
    site2type = atom_sites[site2].type
    # print(site1type,site2type)
    # upper_range = bond_length + bond_length * \
    #     float(config['bond_length_tolerance_factor'])
    # lower_range = bond_length - bond_length * \
    #     float(config['bond_length_tolerance_factor'])

    bond_types_general = known_bond_types()
    # print(bond_types_general)
    # print('--',site1type, site2type, upper_range, lower_range, bond_types_general, bond_length)
    for each in bond_types_general:
        # print(each) # ('C', 'S', 680.0, 1.48)
        # print(each[0], each[0].translate(digits), site1type)      

        if each[0] == site1type and each[1] == site2type:
            if bond_length <= each[3] * (1+float(config['bond_length_tolerance_factor'])):
                return True
        # if each[0].translate(digits) == site2type and each[1].translate(digits) == site1type:
        #     if bond_length <= upper_range:
                # return True
    return False

# Checks for digits within an input string


def contains_digits(s):
    return any(char.isdigit() for char in s)

# Assigns bond types


def bond_type_assignment(siteval):
    atom_number = 0
    for bonded_atom in atom_sites[siteval].bonds:
        bond_type_counter = 0
        for bond_types in known_bond_types():
            if (bond_types[0] == atom_sites[bonded_atom].type and bond_types[1] == atom_sites[siteval].type) or (
                    bond_types[1] == atom_sites[bonded_atom].type and bond_types[0] == atom_sites[siteval].type):
                atom_sites[siteval].bond_types[atom_number] = bond_type_counter
                break
            bond_type_counter += 1
        atom_number += 1

# Removes bonds that are left untyped, we do not have information for this type of bonding behavior.
# This can happen if the bond tolerances are raised too high.


def bond_untype_removal(i, errors_list):
    indexes = []

    for ii in range(0, len(atom_sites[i].bonds)):
        if atom_sites[i].bond_types[ii] == None:
            # Remove This Bond
            indexes.append(ii)

            if not sorted([atom_sites[i].type, atom_sites[
                          atom_sites[i].bonds[ii]].type]) in errors_list:
                errors_list.append(
                    sorted([atom_sites[i].type, atom_sites[atom_sites[i].bonds[ii]].type]))

    if not indexes == []:
        for index in sorted(indexes, reverse=True):
            del atom_sites[i].bond_types[index]
            del atom_sites[i].bonds[index]
            del atom_sites[i].bond_lengths[index]

    return errors_list
# Assigns atom types


def type_assignment(siteval, level, order_assign_number):
    atom_types_general, bond_types_general = known_atom_types_general(), known_bond_types()
    i, types, dependencies = 0, [
        None] * len(atom_types_general), [None] * len(atom_types_general)
    for each in atom_types_general:
        types[i] = each[0]
        dependencies[i] = each[4].values()
        i += 1

    # Check if atomic species is correct for the order asign number.
    if types[order_assign_number].translate(digits) == sites[
            siteval].species_string:
        bonded_atoms, nn_types = atom_sites[siteval].bonds, []

        # Level 0: (None) Assign higher level types based on basic species,
        # requires no bonds.
        if level == 0:
            atom_sites[siteval].type = types[order_assign_number]

        for each in bonded_atoms:
            # Level 1: Remove type from atom to get species
            if level == 1:
                # Level 1: Assign based on atomic species.
                nn_types.append(atom_sites[each].type.translate(digits))
            # Level 2: Assigned based on type attached.
            if level == 2:
                nn_types.append(atom_sites[each].type)

        # Case 1: Lists are Identical, assign type
        if sorted(nn_types) == sorted(dependencies[order_assign_number]):
            # Assign new type
            atom_sites[siteval].type = types[order_assign_number]

        # Case 2: Lists are possible, determine if larger than prior
        # diescription
        if sorted(dependencies[order_assign_number]) < sorted(nn_types):
            if contains_digits(atom_sites[siteval].type):
                # Assign highest order assign number
                if len(
                    dependencies[order_assign_number]) >= len(
                    dependencies[
                        atom_type_NUMBER(
                            atom_sites[siteval].type)]):
                    atom_sites[siteval].type = types[order_assign_number]
            else:
                atom_sites[siteval].type = types[order_assign_number]
# Counts number of bonds


def count_BONDS():
    bonds, site_number, bonds_added = 0, 0, []
    for each in atom_sites:
        for every in each.bonds:
            if sorted([site_number, every]) not in bonds_added:
                bonds_added.append(sorted([site_number, every]))
                bonds += 1
        site_number += 1
    return bonds

# Counts number of angles


def count_ANGLES():
    angles = 0
    for each in atom_sites:
        angles += len(each.angles)
    return angles

# Counts number of dihedrals


def count_dihedrals():
    dihedrals = 0
    for each in atom_sites:
        dihedrals += len(each.dihedrals)
    return dihedrals

# Adds a specified number of blank lines


def save_BLANK_LINES(number_blank_lines, f):
    for i in range(0, number_blank_lines):
        print('', file=f)

# Generates the data file


def generate_DATA_FILE():
    """
    See: http://lammps.sandia.gov/doc/2001/data_format.html
    for an example of data file.
    """

    fil = filename.split(".", 1)
    f = open(str(fil[0] + ".data"), 'w')

    # LAMMPS Description           (1st line of file)
    print('LAMMPS DATA FILE FOR ' + fil[0] + ' in metal units', file=f)

    save_BLANK_LINES(1, f)

    print(str(len(atom_sites)) + " atoms", file=f)
    print(str(count_BONDS()) + " bonds", file=f)
    print(str(count_ANGLES()) + " angles", file=f)
    num_dihedrals = count_dihedrals()
    if num_dihedrals > 0:
        print(str(num_dihedrals) + " dihedrals", file=f)
    #print(len(atom_sites)+ " impropers", file=f)

    save_BLANK_LINES(1, f)

    print(str(len(known_atom_types())) + " atom types", file=f)
    print(str(len(known_bond_types())) + " bond types", file=f)
    if not onlyFindBonds:
        print(str(len(known_angle_types())) + " angle types", file=f)
        print(str(len(known_torsion_types())) + " dihedral types", file=f)
    print("0 improper types", file=f)

    save_BLANK_LINES(1, f)
    # print("ITEM: BOX BOUNDS xy xz yz pp pp pp", file=f)
    # save_BLANK_LINES(1, f)

    # X = config["lattice_parameters"]["vectors"]['i']
    # Y = config["lattice_parameters"]["vectors"]['j']
    # Z = config["lattice_parameters"]["vectors"]['k']
    R = structure.lattice.matrix
    # print(R[0])
    X = R[0]
    Y = R[1]
    Z = R[2]


    # alpha = math.radians(config["lattice_parameters"]["angles"]['alpha'])
    # beta = math.radians(config["lattice_parameters"]["angles"]['beta'])
    # gamma = math.radians(config["lattice_parameters"]["angles"]['gamma'])
    xlo = 0
    ylo = 0
    zlo = 0
    xhi = X[0]
    yhi = Y[1]
    zhi = Z[2]
    xy = Y[0]
    xz = Z[0]
    yz = Z[1]
    # xhi = LA.norm(X)
    # xy = LA.norm(Y) * math.cos(gamma)
    # yhi = LA.norm(Y) * math.sin(gamma)
    # xz = LA.norm(Z) * math.cos(beta)
    # Was yz changed to xz, veryify this was a typo. -Pueschel
    # yz = (np.dot(Y, Z) - xy * xz) / yhi
    # zhi = math.sqrt(math.pow(LA.norm(Z), 2) -
                    # math.pow(xz, 2) - math.pow(yz, 2))
    # -0.5 0.5 xlo xhi       #(for periodic systems this is box size,
    print(f'{str(xlo)} {str(xhi)} xlo xhi', file=f)
    # -0.5 0.5 ylo yhi       # for non-periodic it is min/max extent of atoms)
    print(f'{str(ylo)} {str(yhi)} ylo yhi', file=f)
    # -0.5 0.5 zlo zhi       #(do not include this line for 2-d simulations)
    print(f'{str(zlo)} {str(zhi)} zlo zhi', file=f)
    print(f'{str(xy)} {str(xz)} {str(yz)} xy xz yz', file=f)



    # ATOMS
    save_BLANK_LINES(1, f)
    print("Atoms # full", file=f)
    save_BLANK_LINES(1, f)
    i = 1
    katg = known_atom_types_general()
    molecule_number = 1
    sites = structure.sites
    for each in atom_sites:
        type_of_atom = atom_type_NUMBER(each.type)
        coords = sites[i - 1].coords
        katg_each = katg[type_of_atom]
        print(str(i) +
              " " +
              str(molecule_number) +
              " " +
              str(type_of_atom +
                  1) +
              " " +
              str(katg_each[1]) +
              " " +
              str(coords[0]) +
              " " +
              str(coords[1]) +
              " " +
              str(coords[2]), file=f)
            #   " 0 0 0", file=f)
        i += 1

    # Velocities, ADD IF NEEDED
    # save_BLANK_LINES(1, f)
    # print("Velocities", file=f)
    # save_BLANK_LINES(1, f)
    # i = 1
    # for each in atom_sites:
    #     print(str(i) + " 0 0 0", file=f)
    #     i += 1


    save_BLANK_LINES(1, f)

    print("Masses", file=f)
    save_BLANK_LINES(1, f)
    gen, i = known_atom_types_general(), 1
    for each in gen:
        print(str(i) + " " + str(each[3]), file=f)
        i += 1

    # save_BLANK_LINES(1, f)
    # print("Nonbond Coeffs",f)
    # Add this feature if you require it.

    # Bond Coeffs, ENERGY [eV/A^2], LENGTH [A]
    save_BLANK_LINES(1, f)
    print("Bond Coeffs", file=f)
    save_BLANK_LINES(1, f)
    bonds, i = known_bond_types(), 1
    for each in bonds:
        # bond K r0
        print(str(i) + " " + str(each[2]/2/kcalPermolToeV) + " " + str(each[3]), file=f)  # ***Note divide by two for putting in terms of lammps
        i += 1
    if not onlyFindBonds:
            
        #ANGLES (BENDING), ENERGY [eV/rad^2], THETA [deg]
        save_BLANK_LINES(1, f)
        print("Angle Coeffs", file=f)
        save_BLANK_LINES(1, f)
        angles, i = known_angle_types(), 1
        for each in angles:
            # angle K theta0
            print(str(i) + " " + str(each[3]/2/kcalPermolToeV) + " " + str(each[4]), file=f) # ***Note divide by two for putting in terms of lammps
            i += 1

        # PROPER TORSIONS, ENERGY [eV], ANGLE [deg]
        # ***Note divide by two for putting in terms of lammps
        # AND the 3 for the n harmonic dihedrals for dreiding nafion
        save_BLANK_LINES(1, f)
        print("Dihedral Coeffs", file=f)
        save_BLANK_LINES(1, f)
        torsions, i = known_torsion_types(), 1
        for each in torsions:
            # dihedral K d n
            print(str(i) + " " + str(each[4]/2/kcalPermolToeV) + " " +
                str(each[6]) + " 3", file=f) 
            i += 1

    # Bonds
    save_BLANK_LINES(1, f)
    print("Bonds", file=f)
    save_BLANK_LINES(1, f)
    i = 1
    site_number = 1
    bonds_added = []
    for each in atom_sites:
        for every in each.bonds:
            if sorted([site_number, every + 1]) not in bonds_added:
                print(str(i) + " " + str(each.bond_type(every) + 1) +
                      " " + str(site_number) + " " + str(every + 1), file=f)
                bonds_added.append(sorted([site_number, every + 1]))
                i += 1
        site_number = site_number + 1
    if not onlyFindBonds:
        # Angles
        save_BLANK_LINES(1, f)
        print("Angles", file=f)
        save_BLANK_LINES(1, f)
        i = 1
        site_number = 1
        for each in atom_sites:
            ii = 0
            for every in each.angles:
                print(f'{str(i)} {str(each.angle_types[ii] + 1)} {str(site_number)} '
                f'{str(every[1] + 1)} {str(every[2] + 1)} # '
                f'{atom_sites[site_number-1].type} {atom_sites[every[1]].type} '
                f'{atom_sites[every[2]].type}', file=f)  # index for atom sites is ahead by 1
                ii = ii + 1
                i += 1
            site_number = site_number + 1

        # Dihedrals
        save_BLANK_LINES(1, f)
        print("Dihedrals", file=f)
        save_BLANK_LINES(1, f)
        i = 1
        site_number = 1
        for each in atom_sites:
            ii = 0
            for every in each.dihedrals:
                print(f'{str(i)} {str(each.dihedral_types[ii] + 1)} {str(site_number)} '
                f'{str(every[1] + 1)} {str(every[2] + 1)} {str(every[3] + 1)} # '
                f'{atom_sites[site_number-1].type} {atom_sites[every[1]].type} '
                f'{atom_sites[every[2]].type} {atom_sites[every[3]].type}', file=f)  # index for atom sites is ahead by 1
                ii = ii + 1
                i += 1
            site_number = site_number + 1

    # Impropers
    # save_BLANK_LINES(1, f)
    # print("Impropers", file=f)
    # If you need this functionality, feel free to add this feature and then
    # request a merge to the main branch.


def generate_VDW_DATA_FILE():
    # ** note conversion from D0 to epsilon and R0 to sigma
    # mixing by geometric mean
    # R0=2^1/6*sigma
    # D0=4*epsilon XX seems like it is D0 is epsilon
    # don't need to mix now that lammps can just do it BUT not with hybrid
    atom_types_general = known_atom_types_general()
    f = open(str("VDW_LAMMPS"), 'w')
    for i in range(len(atom_types_general)):
        lj_i = atom_types_general[i][2]
        # f.write('pair_coeff ' + str(i + 1) + ' ' + str(i + 1) + \
        # ' ' + str(float(lj_i[1])/4/kcalPermolToeV) + ' ' + str(lj_i[0]/1.123) + '\n')
        for j in range(i,len(atom_types_general)):
            lj_j = atom_types_general[j][2]
            mixed_epsilon = np.sqrt(float(lj_i[1]) * float(lj_j[1]))
            mixed_sigma = np.sqrt(float(lj_i[0]) * float(lj_j[0]))
            # mixed_sigma = (float(lj_i[0]) + float(lj_j[0])) / float(2)
            f.write('pair_coeff ' + str(i + 1) + ' ' + str(j + 1) + \
                    ' lj/cut/coul/cut ' + str(mixed_epsilon/kcalPermolToeV) + ' ' + str(mixed_sigma/1.123) + '\n')  
    f.close()


class StructureSite:
    # bonds is a list of atoms that the atom is bonded to

    def __init__(self, type):
        self.type = type
        self.bonds = []  # list(set(bonds))
        self.bond_types = []
        self.bond_lengths = []
        self.angles = []
        self.angle_types = []
        self.dihedrals = []
        self.dihedral_types = []

    def add_bond(self, bond, type, bond_length):
        if bond not in self.bonds:
            self.bonds.append(bond)
            self.bond_types.append(type)
            self.bonds = self.bonds
            self.bond_types = self.bond_types
            self.bond_lengths.append(bond_length)

    def bond_type(self, atom2):
        i = 0
        bond_types = self.bond_types
        for each in self.bonds:
            if atom2 == each:
                return bond_types[i]
            i += 1
        return None

    def bonded(self, atom2):
        for each in self.bonds:
            if atom2 == each:
                return True

    def add_angle(self, atom1, atom2, atom3, angle_type):
        if [atom1, atom2, atom3] not in self.angles:
            if atom_sites[atom1].bonded(atom2):# or atom_sites[atom2].bonded(atom1):
                if atom_sites[atom2].bonded(atom3):# or atom_sites[atom3].bonded(atom2):
                    if angle_type is not None:
                        self.angle_types.append(angle_type)
                        self.angles.append([atom1, atom2, atom3])

    def check_if_angle(self, atom1, atom2, atom3):
        atoms_list = (atom_sites[atom1].type, atom_sites[
                      atom2].type, atom_sites[atom3].type)
        # print(atom1,atom2,atom_sites[atom1].bonded(atom2))
        if atom_sites[atom1].bonded(atom2):# or atom_sites[atom2].bonded(atom1):
            if atom_sites[atom2].bonded(atom3):# or atom_sites[atom3].bonded(atom2):
                # Check if angle type is defined in our list of angle types.
                known_angle_type = known_angle_types()
                angle_type_list = []
                for each in known_angle_type:
                    angle_type_list.append((each[0], each[1], each[2]))
                    # add permutation of arrangement
                    # angle_type_list.append((each[0], each[1], each[2]))
                # print(angle_type_list)
                # print(atom1, atom2, atom3)
                # print(atoms_list)
                count = 0
                for each in angle_type_list:
                    # print(atoms_list, each) #XX
                    if atoms_list == each:
                        return (True, count)
                    count += 1
        return (False, None)

    def add_dihedral(self, atom1, atom2, atom3, atom4, dihedral_type):
        # Assuming Only Proper Torsions for now
        if [atom1, atom2, atom3, atom4] not in self.dihedrals:
            if atom_sites[atom1].bonded(atom2):
                if atom_sites[atom2].bonded(atom3):
                    if atom_sites[atom3].bonded(atom4):
                        self.dihedral_types.append(dihedral_type)
                        self.dihedrals.append([atom1, atom2, atom3, atom4])

    def check_if_dihedral(self, atom1, atom2, atom3, atom4):
        # print(atom1, atom2, atom3, atom4)
        # Assuming Only Proper Torsions for now
        atoms_sorted_list = (
            atom_sites[atom1].type,
            atom_sites[atom2].type,
            atom_sites[atom3].type,
            atom_sites[atom4].type)
        if atom_sites[atom1].bonded(atom2):
            if atom_sites[atom2].bonded(atom3):
                if atom_sites[atom3].bonded(atom4):
                    # Check if dihedral type is defined in our list of dihedral types.
                    # i,j,k,l for proper torsions
                    known_torsion_type = known_torsion_types()
                    dihedral_type_list = []
                    for each in known_torsion_type:
                        dihedral_type_list.append(
                            (each[0], each[1], each[2], each[3]))
                    count = 0
                    for each in dihedral_type_list:
                        # print(atoms_sorted_list, each)
                        if atoms_sorted_list == each:
                            return (True, count)
                        count += 1
        return (False, None)

#-------------- Import Config File -----------
# if len(sys.argv) < 2:
# 	print('Usage: LDFG.py <dump File name>')
# 	exit(1)

# filename = sys.argv[1]  #read dump filename

# global config
# with open("config", 'r') as ymlfile:
#     # print(ymlfile)
#     config = yaml.safe_load(ymlfile)

global kcalPermolToeV
kcalPermolToeV = 23.06 #div by to go from kcal/mol to eV.  1 eV to 23.06 kcal/mol

#-------------- Load in Structure -------------

#TESTING -----
from pymatgen.io.lammps import outputs
from pymatgen.core import SETTINGS, Element, Lattice, Structure
import numpy as np
# from pymatgen.analysis.diffusion import analyzer
from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer
import glob
# frames=outputs.parse_lammps_dumps('/global/homes/k/kamron/Scratch/NNmd/Poly/MC11b_lammps6merWPtandO2inNVTwsmallModel/temp3largeandfloat32/pt2/test100fsInt.dump')
# frames=outputs.parse_lammps_dumps('test100fsInt.dump')
# frames=outputs.parse_lammps_dumps('pt2mc11bsorted100fs.dump')
# frames=outputs.parse_lammps_dumps('test100fsInt7.dump')

# frames=outputs.parse_lammps_dumps('../test5fsInt.dump')
#TESTING ^^^^

# frames=outputs.parse_lammps_dumps(filename)


# df = pd.DataFrame(columns=('Pt_H','Pt_O','O_O','O_H'))

structures = []
fnames = sorted(glob.glob('*.dump')) # [:2]

for filename in fnames:
    print(filename)
    frames=outputs.parse_lammps_dumps(filename)
    # loop over every frame
    for iframe, frame in enumerate(frames):
        pass
        # break
    iframe += 1 # to set to 1 basis
    endFrame = iframe-iframe%10
    print('totalframes ',iframe, ' end by ', endFrame)  
    
    frames=outputs.parse_lammps_dumps(filename) 
    #^^ need to reload to enable generator to reset - must be a better way to do this

    # test100fsInt1.dump
    # totalframes  250  end by  250
    # test100fsInt2.dump
    # totalframes  250  end by  250
    # test100fsInt3.dump
    # totalframes  226  end by  220
    # test100fsInt4.dump
    # totalframes  250  end by  250
    # test100fsInt5.dump
    # totalframes  225  end by  220
    # test100fsInt6.dump
    # totalframes  250  end by  250
    # test100fsInt7.dump
    # totalframes  224  end by  220
    # test100fsInt8.dump
    # totalframes  231  end by  230


    for iframe, frame in enumerate(frames):
        # if iframe==60: break
        if iframe+1 == endFrame: break
        print('current frames/total =', iframe+1, '/' , endFrame) # 1 basis
        # sys.exit(1)
        coords=np.stack((frame.data.x.to_numpy(),frame.data.y.to_numpy(),frame.data.z.to_numpy())).T
        # print(coords[6579])
        # atomic_symbols = frame.data.element
        type_dict={1:'C',2:'F',3:'H',4:'O',5:'Pt',6:'S'}
        # print('frame.data.type', frame.data.type)
        # account for if an atom disappears
        atomic_symbols = [type_dict[i] for i in frame.data.type]
        lattice=frame.box.to_lattice()
        global structure
        structure = Structure(
            lattice,
            atomic_symbols,
            coords,
            to_unit_cell=False,
            validate_proximity=False,
            coords_are_cartesian=True,
        )
        structures.append(structure)

    # print(coords) # matches dump
    # print(lattice) #

# print(structures[0])

# ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 1,100,initial_disp=None,initial_structure=structures[0],c_ranges=[(0.,.25)])



# cmax=116. # old
cmax = 99.

# cs = np.array([[30.,90.]])/cmax
# cs = np.array([[8.,30.],[30.,90.],[90.,112.]])/cmax
cs = np.array([[10.,35.],[35.,65.],[65.,90.]])/cmax # for new structure
cs = cs.tolist()
# allMSDs = []

for i,c in enumerate(cs):
    print('Region:', c)
    # from_structures(structures,specie,temperature,time_step,step_skip,initial_disp=None,initial_structure=None,**kwargs)
    # ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 1,100,initial_disp=None,initial_structure=structures[0],c_ranges=[(0.25,0.75)])
    # ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 1,100,initial_disp=None,initial_structure=structures[0],c_ranges=[(25.,75.)])
    # ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 1,100,initial_disp=None,initial_structure=structures[0],c_ranges=[(0.,25.),(25.,75.)])
    # ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 1,100,initial_disp=None,initial_structure=structures[0],c_ranges=[(0.,.25)])
    ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 100,1,initial_disp=None,initial_structure=structures[0],c_ranges=[c]) # ,smoothed=False
    # ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 5,10,initial_disp=None,initial_structure=structures[0],c_ranges=[c]) 
    print('diffusivity_c_range_components ',ob.diffusivity_c_range_components)
    print('std ',ob.diffusivity_c_range_components_std_dev)
    print('total ', ob.diffusivity_c_range)
    # print(ob.diffusivity_c_range_components)

    # ob=DiffusionAnalyzer.from_structures(structures, 'H', 300, 1,100,None,structures[0])
    # ob.msd_components

    # print(ob.get_summary_dict())
    plt, df = ob.get_msd_plot(mode='ranges')
    # allMSDs.append(df)
    df.to_csv(f'{i}MSD.csv')
    plt.savefig(f'{i}MSD.png')
    plt.close()


# results
# [0.06896551724137931, 0.25862068965517243]
# [0.00000000e+00 0.00000000e+00 1.44035715e-06]
# [0.00000000e+00 0.00000000e+00 2.00831611e-09]
# 0.00021766855659840408
# [0.25862068965517243, 0.7758620689655172]
# [0.00000000e+00 0.00000000e+00 1.85207845e-06]
# [0.00000000e+00 0.00000000e+00 2.47591558e-09]
# 0.00027753267740046383
# [0.7758620689655172, 0.9655172413793104]
# [0.00000000e+00 0.00000000e+00 1.16148073e-06]
# [0.00000000e+00 0.00000000e+00 1.99130287e-09]
# 0.0003285887615171238

# ob.REDIRECT                         ob.conductivity_components_std_dev  ob.dt                               ob.get_msd_plot(                    ob.msd                              ob.time_step
# ob.as_dict()                        ob.conductivity_std_dev             ob.export_msdt(                     ob.get_summary_dict(                ob.msd_components                   ob.to_json()
# ob.avg_nsteps                       ob.corrected_displacements          ob.framework_indices                ob.haven_ratio                      ob.plot_msd(                        ob.unsafe_hash()
# ob.chg_conductivity                 ob.diffusivity                      ob.from_dict(                       ob.indices                          ob.smoothed                         ob.validate_monty(
# ob.chg_conductivity_std_dev         ob.diffusivity_components           ob.from_files(                      ob.lattices                         ob.specie                           
# ob.chg_diffusivity                  ob.diffusivity_components_std_dev   ob.from_structures(                 ob.max_framework_displacement       ob.sq_disp_ions                     
# ob.chg_diffusivity_std_dev          ob.diffusivity_std_dev              ob.from_vaspruns(                   ob.max_ion_displacements            ob.step_skip                        
# ob.conductivity                     ob.disp                             ob.get_drift_corrected_structures(  ob.min_obs                          ob.structure                        
# ob.conductivity_components          ob.drift                            ob.get_framework_rms_plot(          ob.mscd                             ob.temperature       
sys.exit(1)


'''
    # structure_pos = Poscar.from_file(filename)
    # structure = structure_pos.structure
    # print(structure.lattice) XX
    sites = structure.sites
    position_order_assignment = type_assignment_execution_order()
    global atom_sites
    atom_sites = [None] * len(sites)
    for i in range(len(sites)):
        # This generates the class from above for site_in_Structure
        atom_sites[i] = StructureSite(str(sites[i].species_string))




    #-------------- Intialize Lists -------------
    atom_types_general = known_atom_types_general()
    i, types, dependencies = 0, [
        None] * len(atom_types_general), [None] * len(atom_types_general)

    for each in atom_types_general:
        types[i] = each[0]
        dependencies[i] = each[4].values()
        # print(each[4]) 0:none
        i += 1

    #--------- Find Nearest Neighbors -------------
    # nn_sites = structure.get_all_neighbors(max_bond_length(), include_index=True)
    nn_sites = structure.get_all_neighbors(r=max_bond_length())
    # nn_sites[65]  S 
    # nn_sites[36]  OOOC

    # print('s1 element',structure.species[27496-1])
    # print('nn_sites',nn_sites[27496-1])
    # # get all neigh returns [(site, dist) ...]
    # print('s2',nn_sites[27496-1][0][2])
    # print('bond dist',nn_sites[27496-1][0][1])
    # print('[coord] element',nn_sites[27496-1][0][0])
    # outputs
    # s1 element F
    # nn_sites [PeriodicSite: C (41.2143, 47.1704, 30.9792) [0.6546, 0.6827, 0.2667], PeriodicSite: C (42.4738, 46.4372, 30.3612) [0.6746, 0.6721, 0.2613], PeriodicSite: F (42.0048, 45.2571, 29.8931) [0.6672, 0.6550, 0.2573], PeriodicSite: F (43.3704, 46.1663, 31.3109) [0.6889, 0.6682, 0.2695]]
    # s2 27491
    # bond dist 2.4623963247211025
    # [coord] element [41.2143 47.1704 30.9792] C


    #------------- Bonded Atoms Assignment ----------------
    # Iterates through each position assignment for Bond Assignment.
    excludeFromBonds = config['excludeFromBonds']
    try:
        onlyFindBonds = config['onlyFindBonds']
    except:
        onlyFindBonds = False
    # print(config)
    print('onlyFindBonds-',onlyFindBonds)

    # len(np.argwhere(np.array([str(i) for i in structure.species]) == 'Pt'))
    # len(np.argwhere(np.array([str(i) for i in structure.species]) == 'Pt' )
    zVals = np.array(structure.cart_coords)[:,2]
    sitesNearSurfandSurf = np.argwhere((zVals<8.)&(6.<zVals))
    # print(sitesNearSurfandSurf[0][0])
    # sys.exit(1)

    # PtTopLayerSites = []
    # for siteval in range(len(sites)):
    #     if structure.species[siteval] == 'Element Pt' and 6.<structure.cart_coords[siteval][2]<7.:
    #         PtTopLayerSites.append(siteval)
    # print(PtTopLayerSites)

    H_PtTotal = 0
    O_PtTotal = 0
    O_OTotal = 0
    O_HTotal = 0
    OsiteOnSurf = []

    # siteval = 0
    for siteval in sitesNearSurfandSurf:
        siteval = siteval[0]

        if str(structure.species[siteval])=='H':
            each_site = nn_sites[siteval]      
            for each in each_site:
                site2 = each[2]
                atDistance = each[1]
                if str(structure.species[site2]) == 'Pt':
                    # print('H Pt dist ', atDistance)
                    if site_bonded(siteval, site2, atDistance):  # s1#, s2#,  length bw atoms
                        # print('H Pt bond')
                        H_PtTotal += 1

        if str(structure.species[siteval])=='O':
            each_site = nn_sites[siteval]      
            for each in each_site:
                site2 = each[2]
                atDistance = each[1]
                if str(structure.species[site2]) == 'Pt':
                    # print('H Pt dist ', atDistance)
                    if site_bonded(siteval, site2, atDistance):  # s1#, s2#,  length bw atoms
                        # print('O Pt bond')
                        O_PtTotal += 1
                        OsiteOnSurf.append(siteval)

    # go through O on Pt surface 
    for siteval in OsiteOnSurf:  
        # double check still O
        if str(structure.species[siteval])=='O':
            each_site = nn_sites[siteval]      
            for each in each_site:
                site2 = each[2]
                atDistance = each[1]

                if str(structure.species[site2]) == 'O':
                    # print('H Pt dist ', atDistance)
                    if site_bonded(siteval, site2, atDistance):  # s1#, s2#,  length bw atoms
                        # print('O Pt bond')
                        O_OTotal += 1
                
                if str(structure.species[site2]) == 'H':
                    # print('H Pt dist ', atDistance)
                    if site_bonded(siteval, site2, atDistance):  # s1#, s2#,  length bw atoms
                        # print('O Pt bond')
                        O_HTotal += 1

                        
                        


    df.loc[len(df.index)] = [H_PtTotal,O_PtTotal,O_OTotal,O_HTotal] # append to df

    print(df)
    df.to_csv('data.csv',index=False)
    # print(f"H_PtTotal {H_PtTotal}")
    # print(f"O_PtTotal {O_PtTotal}")
    # print(f"O_OTotal {O_OTotal}")
    # print(f"O_HTotal {O_HTotal}")

'''
