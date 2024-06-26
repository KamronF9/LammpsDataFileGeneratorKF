#LAMMPS Data File Generator
Generate data file for lammps using force field and a POSCAR (VASP) like input file.

##1. What does it do?
Quickly generate a LAMMPS data file to define a force field for a periodic system and Van Der Waals interactional parameters data file.
This includes:
* Bonds
* Angles
* Proper Torsions (dihedrals)
* Improper Torsions [Future]

The tolerances can be adjusted in the config file, since it is likely that unoptimized structures will have errors in bond lengths. 

##2. Required Dependencies
Make sure you have installed pymatgen. This can be done:
    `easy_install pymatgen`
or 
    `pip install pymatgen==2023.3.23`
Alternatively, go to their [website](http://pymatgen.org/) for installation instructions.

##3. Quick-start Guide
Obtain all of the necessary input files; 
>1. Config containing proper atom, bond, angle, dihedral, and filename.
>2. Molecule/System defined in a POSCAR (VASP). Additional support for just xyz is being implemented.
>3. Ensure you have pymatgen and python 2.7+
>4. The example folder includes the files to generate the data file for a CHA Zeolite.
To run simply edit the input files and run with them in the same folder;
`python LDFG.py`

