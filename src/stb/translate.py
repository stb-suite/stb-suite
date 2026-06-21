#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
try:
    from importlib.metadata import version as _pkg_version
    VERSION = _pkg_version("stb_suite")
except Exception:
    VERSION = "1.9.5"    
from stb.cli import COLORS, color_text, show_intro

import os
import sys
import warnings
import subprocess
import argparse
import textwrap
from typing import List, Dict
import argparse
from ase.io import read as ase_read
import numpy as np

def periodic_table():
    element_atomicnumber = {
        "H": 1,   
        "He": 2,  
        "Li": 3,  
        "Be": 4,  
        "B": 5,   
        "C": 6,   
        "N": 7,   
        "O": 8,   
        "F": 9,   
        "Ne": 10, 
        "Na": 11, 
        "Mg": 12, 
        "Al": 13, 
        "Si": 14, 
        "P": 15,  
        "S": 16,  
        "Cl": 17, 
        "Ar": 18, 
        "K": 19,  
        "Ca": 20, 
        "Sc": 21, 
        "Ti": 22, 
        "V": 23,  
        "Cr": 24, 
        "Mn": 25, 
        "Fe": 26, 
        "Co": 27, 
        "Ni": 28, 
        "Cu": 29, 
        "Zn": 30, 
        "Ga": 31, 
        "Ge": 32, 
        "As": 33, 
        "Se": 34, 
        "Br": 35, 
        "Kr": 36, 
        "Rb": 37, 
        "Sr": 38, 
        "Y": 39,  
        "Zr": 40, 
        "Nb": 41, 
        "Mo": 42, 
        "Tc": 43, 
        "Ru": 44, 
        "Rh": 45, 
        "Pd": 46, 
        "Ag": 47, 
        "Cd": 48, 
        "In": 49, 
        "Sn": 50, 
        "Sb": 51, 
        "Te": 52, 
        "I": 53,  
        "Xe": 54, 
        "Cs": 55, 
        "Ba": 56, 
        "La": 57, 
        "Ce": 58, 
        "Pr": 59, 
        "Nd": 60, 
        "Pm": 61, 
        "Sm": 62, 
        "Eu": 63, 
        "Gd": 64, 
        "Tb": 65, 
        "Dy": 66, 
        "Ho": 67, 
        "Er": 68, 
        "Tm": 69, 
        "Yb": 70, 
        "Lu": 71, 
        "Hf": 72, 
        "Ta": 73, 
        "W": 74,  
        "Re": 75, 
        "Os": 76, 
        "Ir": 77, 
        "Pt": 78, 
        "Au": 79, 
        "Hg": 80, 
        "Tl": 81, 
        "Pb": 82, 
        "Bi": 83, 
        "Po": 84, 
        "At": 85, 
        "Rn": 86, 
        "Fr": 87, 
        "Ra": 88, 
        "Ac": 89, 
        "Th": 90, 
        "Pa": 91, 
        "U": 92,  
        "Np": 93, 
        "Pu": 94, 
        "Am": 95, 
        "Cm": 96, 
        "Bk": 97, 
        "Cf": 98, 
        "Es": 99, 
        "Fm": 100,
        "Md": 101,
        "No": 102,
        "Lr": 103,
        "Rf": 104,
        "Db": 105,
        "Sg": 106,
        "Bh": 107,
        "Hs": 108,
        "Mt": 109,
        "Ds": 110,
        "Rg": 111,
        "Cn": 112,
        "Nh": 113,
        "Fl": 114,
        "Mc": 115,
        "Lv": 116,
        "Ts": 117,
        "Og": 118 
    }

    # Dictionary with all atomic numbers table (invert Dictionary)
    atomicnumber_element = {str(v): k for k, v in element_atomicnumber.items()}
    return element_atomicnumber, atomicnumber_element

def readfile(filedata):
    with open(filedata, 'r') as fil:
        data = [line.split() for line in fil if line.strip()
                ]
    return data

def dic_atoms_position(atomsposition):
    dic_atomspos = {}
    for elt in atomsposition:
        if elt[0] not in dic_atomspos:
            dic_atomspos[elt[0]] = []
        dic_atomspos[elt[0]].append([elt[1], elt[2], elt[3]])
    return dic_atomspos

############# Extract Data Functions ######################

# This Function define the number of atoms types for xyz file
def getatomsandvectors_xyz(dataxyz, latticedata):
    element, atomicnumber = periodic_table()
    xyz = readfile(dataxyz)[2:]
    lattice = readfile(latticedata)
  
    latticeparameter = lattice[0][0]
    typevectors = 'Cartesian'
    vectors = [lattice[1], lattice[2], lattice[3]]
  
    atoms = []
    j = 1
    getatoms = []
    for i in range(len(xyz)):
        if xyz[i][0] not in atoms:
            atoms.append(xyz[i][0])
            getatoms.append([j, element[xyz[i][0]], xyz[i][0]])
            j = j+1
    atomsposition = xyz
    for elem in getatoms:
        cont = 0
        for elxyz in xyz:
            if elem[2] == elxyz[0]:
                cont = cont+1
        elem.append(str(cont))
    atomic_position = dic_atoms_position(atomsposition)
    return typevectors, latticeparameter, vectors, getatoms, atomic_position

# This Function define the number of atoms types for vasp file
def getatomsandvectors_vasp(poscar):
    element, atomicnumber = periodic_table()
    datavasp = readfile(poscar)
    latticeparameter = datavasp[1][0]
    vectors = datavasp[2:5]
    typevectors = datavasp[7][0]
    getatoms = []
    for i in range(len(datavasp[5])):
        getatoms.append([i+1, element[datavasp[5][i]],datavasp[5][i], datavasp[6][i]])
    atomsposition = []
    cont = 8 
    for el in getatoms:
        for i in range(int(el[3])):
            atomsposition.append([el[2], datavasp[cont][0], datavasp[cont][1], datavasp[cont][2]])
            cont = cont+1
    atomic_position = dic_atoms_position(atomsposition)
    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_fhi(input_fhi):
    element, atomicnumber = periodic_table()
    datafhi = readfile(input_fhi)
    vectors = []
    atomdata = {}
    latticeparameter = "1.00"
    for lines in datafhi:
        if lines[0] == "atom_frac":
            atomdata.setdefault(lines[4], [])
            atomdata[lines[4]].append([
                f"{float(lines[1]):.8f}",
                f"{float(lines[2]):.8f}",
                f"{float(lines[3]):.8f}"])
            typevectors = 'Direct'
        elif lines[0] == "atom":
            atomdata.setdefault(lines[4], [])
            atomdata[lines[4]].append([
                f"{float(lines[1]):.8f}",
                f"{float(lines[2]):.8f}",
                f"{float(lines[3]):.8f}"])
            typevectors = 'Cartesian'
        elif lines[0] == 'lattice_vector':
            vectors.append([lines[1], lines[2], lines[3]])
    getatoms = []
    icont = 1
    for el in atomdata:
        getatoms.append([str(icont), str(element[el]),
                        el, str(len(atomdata[el]))])
        icont = icont+1
    atomic_position = atomdata

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_cif(input_cif):
    """
    Extract structure data from a .cif file using the ASE library.
    Returns data in the standard format expected by the script.
    """
    
    # 1. Get local element dictionary
    element, atomicnumber = periodic_table()
    
    try:
        structure = ase_read(input_cif)
    except NameError:
        print(color_text("[ERROR] Function 'ase_read' not found.", 'red'))
        print(color_text("       Make sure 'from ase.io import read as ase_read' is at the top of the script.", 'yellow'))
        sys.exit(1)
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read CIF file with ASE: {e}", 'red'))
        print(color_text("       Make sure 'ase' is installed (pip install ase)", 'yellow'))
        sys.exit(1)

    # 3. Extract lattice vectors — ASE .get_cell() returns full unscaled vectors in Angstroms
    vectors = []
    for vec in structure.get_cell():
        vectors.append([f"{vec[0]:.8f}", f"{vec[1]:.8f}", f"{vec[2]:.8f}"])

    # 4. Extract atomic positions — Cartesian coordinates in Angstroms
    
    symbols = structure.get_chemical_symbols()
    positions = structure.get_positions()    
    
    dicatoms = {} 
    for sym, pos in zip(symbols, positions):
        if sym not in dicatoms:
            dicatoms[sym] = []
      
        dicatoms[sym].append([f"{pos[0]:.8f}", f"{pos[1]:.8f}", f"{pos[2]:.8f}"])
    
    atomic_position = dicatoms

    # 5. Extract species information
    
    getatoms = []
    icont = 1
    
  
    for sym in dicatoms: 
        getatoms.append([
            f"{icont}",            # Indice (ex: '1')
            f"{element[sym]}",
            f"{sym}",
            f"{len(dicatoms[sym])}"
        ])
        icont += 1

    # 6. Set vector type and lattice parameter (ASE returns Cartesian, unscaled)
    typevectors = 'Cartesian'
    latticeparameter = '1.00'

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_siesta(input_siesta):
    element, atomicnumber = periodic_table()
    datasiesta = readfile(input_siesta)
    latticeparameter = "1.00"
    typevectors = "Direct"
    vectors = datasiesta[:3]
    dicatoms = {}
    for pos in datasiesta[4:]:
        if pos[0] not in dicatoms:
            dicatoms[pos[0]] = []
        dicatoms[pos[0]].append([pos[0], pos[1], atomicnumber[pos[1]]])
    getatoms = []
    for atoms in dicatoms:
        getatoms.append([dicatoms[atoms[0]][0][0],
                        dicatoms[atoms[0]][0][1],
                        dicatoms[atoms[0]][0][2],
                        str(len(dicatoms[atoms[0]]))])
    atomic_position = {}
    for position in datasiesta[4:]:
        if atomicnumber[position[1]] not in atomic_position:
            atomic_position[atomicnumber[position[1]]] = []
        atomic_position[atomicnumber[position[1]]].append(
            [position[2], position[3], position[4]])
    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_dftb(input_dftb):
    element, atomicnumber = periodic_table()
    datadftb = readfile(input_dftb)
    vectors = []
    getatoms = []
    atomic_position = {}
    dic_data = {}
    latticeparameter = "1.00"
    if datadftb[0][1] == 'S':
        typevectors = 'Cartesian'
    elif datadftb[0][1] == 'F':
        typevectors = 'Direct'
    else:
        print("[FAIL] Type of coordinate not define: Only S or F are accepted")
        exit()
    for i in range(1, 4):
        vectors.append([f"{float(datadftb[-i][0]):.8f}",
                        f"{float(datadftb[-i][1]):.8f}",
                        f"{float(datadftb[-i][2]):.8f}"])
    icont = 0
    dic_data = {f"{i + 1}": elem for i, elem in enumerate(datadftb[1])}
    for line in datadftb[2:-4]:
        atomic_position.setdefault(dic_data[line[1]], [])
        atomic_position[dic_data[line[1]]].append(
            [f"{float(line[2]):.8f}", f"{float(line[3]):.8f}", f"{float(line[4]):.8f}"])
    icont = 1
    for el in atomic_position:
        getatoms.append([f"{icont}",
                         f"{element[el]}",
                         f"{el}",
                         f"{len(atomic_position[el])}"])
        icont = icont+1

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_xsf(input_xsf):
    print("[WARNING] Only for PRIMVEC format")
    element, atomicnumber = periodic_table()
    dataxsf = readfile(input_xsf)
    vectors = []
    getatoms = []
    atomic_position = {}
    dic_data = {}
    latticeparameter = "1.00"
    typevectors = 'Cartesian'
    dataxsf = [v for v in dataxsf if not str(v[0]).startswith("#")]
    for j in range(len(dataxsf)):
        if dataxsf[j][0] == 'PRIMVEC':
            for i in range(1, 4):
                vectors.append([f"{float(dataxsf[j+i][0]):.8f}",
                                f"{float(dataxsf[j+i][1]):.8f}",
                                f"{float(dataxsf[j+i][2]):.8f}"])
    for j in range(len(dataxsf)):
        if dataxsf[j][0] == 'PRIMCOORD':
            na = dataxsf[j+1][0]
            for i in range(int(na)):
                atomic_position.setdefault(atomicnumber[dataxsf[j+i+2][0]], [])
                atomic_position[atomicnumber[dataxsf[j+i+2][0]]].append(
                    [f"{float(dataxsf[j+i+2][1]):.8f}",
                     f"{float(dataxsf[j+i+2][2]):.8f}",
                     f"{float(dataxsf[j+i+2][3]):.8f}"])
    icont = 1
    for el in atomic_position:
        getatoms.append([f"{icont}",
                         f"{element[el]}",
                         f"{el}",
                         f"{len(atomic_position[el])}"])
        icont = icont+1

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

def getatomsandvectors_fdf(input_fdf):
    """
    Extract structure data from a SIESTA .fdf file.
    """
    element, atomicnumber = periodic_table()
    datafdf = readfile(input_fdf)

    latticeparameter = "1.00"
    typevectors = None
    vectors = []
    species_map = {}  # maps index (str) -> symbol (str), e.g. '1' -> 'C'
    atomic_position_raw = []
    
    atomic_position = {} 
    getatoms = []

  
    in_vectors_block = False
    in_species_block = False
    in_coords_block = False

    for line in datafdf:
        if not line:
            continue
        
      
        key = line[0].lower()
        
      
        if key == '%block':
            block_name = line[1].lower() 
            if block_name == 'latticevectors':
                in_vectors_block = True
                continue
            elif block_name == 'chemicalspecieslabel':
                in_species_block = True
                continue
            elif block_name == 'atomiccoordinatesandatomicspecies':
                in_coords_block = True
                continue
        
      
        if key == '%endblock':
            in_vectors_block = False
            in_species_block = False
            in_coords_block = False
            continue

      
        if in_vectors_block:
            vectors.append([line[0], line[1], line[2]])
        elif in_species_block:
            species_map[line[0]] = line[2]
        elif in_coords_block:
            atomic_position_raw.append([line[0], line[1], line[2], line[3]])
        
      
        else:
            if key == 'latticeconstant':
                latticeparameter = line[1]
            elif key == 'atomiccoordinatesformat':
                if line[1].lower() == 'fractional':
                    typevectors = 'Direct'
                elif line[1].lower() == 'ang':
                    typevectors = 'Cartesian'
    
    # Default to Direct (Fractional) if format not specified
    if typevectors is None:
        print("[WARNING] AtomicCoordinatesFormat not found. Assuming 'Direct' (Fractional).")
        typevectors = 'Direct'

    for pos in atomic_position_raw:
        x, y, z, species_index = pos
        symbol = species_map[species_index] 
        
        if symbol not in atomic_position:
            atomic_position[symbol] = []
        atomic_position[symbol].append([x, y, z])

    # Build getatoms in the order defined by ChemicalSpeciesLabel
    sorted_species_indices = sorted(species_map.keys(), key=int)

    for index_str in sorted_species_indices:
        symbol = species_map[index_str]
        atomic_num = element[symbol]
        count = len(atomic_position.get(symbol, []))
        
      
        if count > 0:
             getatoms.append([index_str, str(atomic_num), symbol, str(count)])

    return typevectors, latticeparameter, vectors, getatoms, atomic_position

###################### Write functions ################

def writefilefdf(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )

    numberofatoms = 0
    for lin in getatoms:
        numberofatoms = numberofatoms+int(lin[3])
    outfile = []
    outfile.append(
        '# automatic create  using stb-translate (https://github.com/bastoscmo/stb-suite)\n\n')
    outfile.append(f"NumberOfSpecies    {len(getatoms)}")
    outfile.append(f"NumberofAtoms      {numberofatoms}\n\n")
    outfile.append("%block ChemicalSpeciesLabel")
    for atoms in getatoms:
        outfile.append(f" {atoms[0]}   {atoms[1]}   {atoms[2]}")
    outfile.append("%endblock ChemicalSpeciesLabel \n")
    outfile.append(f"LatticeConstant {latticeparameter} Ang \n")
    
    if final_type == 'Direct':
        outfile.append("AtomicCoordinatesFormat  Fractional \n\n")
    if final_type == 'Cartesian':
        outfile.append("AtomicCoordinatesFormat  Ang\n\n")

    outfile.append("%block LatticeVectors")
    for lin in vectors:
        outfile.append(f" {lin[0]}   {lin[1]}   {lin[2]} ")
    outfile.append("%endblock LatticeVectors\n\n")
    outfile.append("%block AtomicCoordinatesAndAtomicSpecies")
    
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"  {position[0]}   {position[1]}   {position[2]}   {elem[0]}  ")

    outfile.append("%endblock AtomicCoordinatesAndAtomicSpecies")
    np.savetxt(outfilename, outfile, fmt='%s')
    return

def writefileposcar(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )

    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    outfile.append(f"{latticeparameter}")
    for lin in vectors:
        outfile.append(f"{lin[0]}   {lin[1]}   {lin[2]} ")
    lineatoms = ""
    linenatoms = ""
    for atoms in getatoms:
        lineatoms = lineatoms + f"{atoms[2]}   "
        linenatoms = linenatoms+f"{atoms[3]}   "
    outfile.append(f"{lineatoms}")
    outfile.append(f"{linenatoms}")
    
    if final_type == 'Direct':
        outfile.append("Direct")
    if final_type == 'Cartesian':
        outfile.append("Cartesian")

    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"{position[0]}   {position[1]}   {position[2]}")

    np.savetxt(outfilename, outfile, fmt='%s')
    return

  
def angle(u, v):
        cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
        

def writefilecif(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    """
    Writes a CIF file in the standard format.
    Note: CIF *always* uses fractional (Direct) coordinates.
    """
    
    if coord_format and coord_format.lower() == 'cartesian':
        print("[WARNING] CIF format requires fractional (Direct) coordinates.")
        print("[INFO]    Ignoring '--coord-format cartesian' and converting to Direct.")
    
    # Force conversion to Direct, regardless of user input
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'direct' # Force 'direct'
    )

    vectors_np = np.array(vectors, dtype=float) * float(latticeparameter)

    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    outfile.append("data_generated")
    outfile.append("_symmetry_space_group_name_H-M   'P 1'")
    outfile.append("_symmetry_Int_Tables_number      1")
    outfile.append("_cell_length_a    {:.8f}".format(np.linalg.norm(vectors_np[0])))
    outfile.append("_cell_length_b    {:.8f}".format(np.linalg.norm(vectors_np[1])))
    outfile.append("_cell_length_c    {:.8f}".format(np.linalg.norm(vectors_np[2])))

    alpha = angle(vectors_np[1], vectors_np[2])
    beta = angle(vectors_np[0], vectors_np[2])
    gamma = angle(vectors_np[0], vectors_np[1])

    outfile.append("_cell_angle_alpha  {:.8f}".format(alpha))
    outfile.append("_cell_angle_beta   {:.8f}".format(beta))
    outfile.append("_cell_angle_gamma  {:.8f}".format(gamma))
    outfile.append(" ")

    outfile.append("loop_")
    outfile.append("_symmetry_equiv_pos_as_xyz")
    outfile.append("  'x, y, z'")
    outfile.append(" ")

    outfile.append("loop_")
    outfile.append("_atom_site_label")
    outfile.append("_atom_site_type_symbol")
    outfile.append("_atom_site_fract_x")
    outfile.append("_atom_site_fract_y")
    outfile.append("_atom_site_fract_z")

    # Since we know final_positions is 'Direct', the logic is simple
    for elem in getatoms:
        for pos_str_list in final_positions[elem[2]]:
            outfile.append(
                f"{elem[2]}   {elem[2]}   {pos_str_list[0]}   {pos_str_list[1]}   {pos_str_list[2]}"
            )

    np.savetxt(outfilename, outfile, fmt='%s')
    return

def writefilexyz(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    if coord_format and coord_format.lower() == 'direct':
        print("[WARNING] XYZ format requires Cartesian (Angstrom) coordinates.")
        print("[INFO]    Ignoring '--coord-format direct' and converting to Cartesian.")

    # Force conversion to Cartesian
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'cartesian' # Force 'cartesian'
    )
    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    outfile.append(f"{sum}")
    outfile.append(comment)
    
    # Since we know final_positions is 'Cartesian', the logic is simple
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"{elem[2]}   {position[0]}   {position[1]}   {position[2]}")

    np.savetxt(outfilename, outfile, fmt='%s')
    return

def writefiledftb(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )
    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    
    if final_type == 'Cartesian':
        outfile.append(f"{sum}   S") # S = Cartesian
    elif final_type == 'Direct':
        outfile.append(f"{sum}   F") # F = Fractional

    atoms = ""
    for i in range(len(getatoms)):
        atoms = atoms + f"{getatoms[i][2]}   "
    outfile.append(atoms)

    icont = 1
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"    {icont}  {elem[0]}  {position[0]}   {position[1]}   {position[2]}")
            icont = icont+1

    outfile.append(f"    0.00000000  0.00000000 0.00000000")
    vectors_np = np.array(vectors, dtype="float")*float(latticeparameter)
    for i in range(3):
        outfile.append(f"    {vectors_np[i][0]:.8f}   {vectors_np[i][1]:.8f}   {vectors_np[i][2]:.8f}")
    np.savetxt(outfilename, outfile, fmt='%s')
    return

def writefilexsf(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    if coord_format and coord_format.lower() == 'direct':
        print("[WARNING] XSF format (PRIMCOORD) requires Cartesian (Angstrom) coordinates.")
        print("[INFO]    Ignoring '--coord-format direct' and converting to Cartesian.")

    # Force conversion to Cartesian
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, 'cartesian' # Force 'cartesian'
    )
    
    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    vectors_np = np.array(vectors, dtype="float")*float(latticeparameter)
    outfile.append('# create by stb-translate ')
    outfile.append(f'# {comment}\n')
    outfile.append('CRYSTAL')
    outfile.append(f"PRIMVEC")
    for i in range(3):
        outfile.append(f"    {vectors_np[i][0]:.8f}   {vectors_np[i][1]:.8f}   {vectors_np[i][2]:.8f}")
    
    # Since we know final_positions is 'Cartesian':
    outfile.append(f"PRIMCOORD")
    outfile.append(f"{sum}  1")
    for elem in getatoms:
        for position in final_positions[elem[2]]:
            outfile.append(
                f"    {elem[1]}   {position[0]}   {position[1]}   {position[2]}")
    
    np.savetxt(outfilename, outfile, fmt='%s')
    return

def writefilefhi(typevectors, latticeparameter, vectors, getatoms, atomsposition, outfilename, coord_format=None):
    
    final_type, final_positions = convert_coordinates(
        typevectors, latticeparameter, vectors, atomsposition, coord_format
    )

    outfile = []
    outfile.append('# automatic create using stb-translate (https://github.com/bastoscmo/stb-suite)')
    sum = 0
    comment = ""
    for el in getatoms:
        comment = comment+f"{el[2]}{el[3]} "
        sum = sum+int(el[3])
    vectors_np = np.array(vectors, dtype="float")*float(latticeparameter)
    for vector in vectors_np:
        outfile.append(f"lattice_vector   {vector[0]:.8f}   {vector[1]:.8f}   {vector[2]:.8f}")
    outfile.append(" ")
    
    if final_type == 'Cartesian':
        for elem in getatoms:
            for position in final_positions[elem[2]]:
                outfile.append(
                    f"atom   {position[0]}   {position[1]}   {position[2]}   {elem[2]}")
    elif final_type == 'Direct':
        for elem in getatoms:
            for position in final_positions[elem[2]]:
                outfile.append(
                    f"atom_frac   {position[0]}   {position[1]}   {position[2]}   {elem[2]}")

    np.savetxt(outfilename, outfile, fmt='%s')
    return

##### NEW HELPER FUNCTION #####
def convert_coordinates(typevectors_in: str, 
                        latticeparameter: str, 
                        vectors: List[List[str]], 
                        atomsposition_in: Dict[str, List[List[str]]], 
                        coord_format_out: str) -> (str, Dict[str, List[List[str]]]):
    """
    Converts coordinates between Direct (Fractional) and Cartesian (Angstrom)
    based on the desired output format.

    Args:
        typevectors_in (str): The input format ('Direct' or 'Cartesian').
        latticeparameter (str): The lattice parameter (scaling factor).
        vectors (List[List[str]]): The lattice vectors.
        atomsposition_in (Dict): The input positions dictionary.
        coord_format_out (str): The desired output format ('direct' or 'cartesian', or None).

    Returns:
        Tuple[str, Dict]: (final_type, final_positions)
                         The final format and the converted positions dictionary.
    """
    
    # 1. Prepare the lattice matrix, scaled by the lattice parameter
    lattice_matrix = np.array(vectors, dtype=float) * float(latticeparameter)
    
    # 2. Determine input and output types
    # If the user did not specify a format, use the input format
    if coord_format_out is None:
        return typevectors_in, atomsposition_in

    # Normalize names to 'Direct' and 'Cartesian' for comparison
    type_in_norm = 'Direct' if typevectors_in.lower() == 'direct' else 'Cartesian'
    type_out_norm = 'Direct' if coord_format_out.lower() == 'direct' else 'Cartesian'

    # 3. Check if conversion is needed
    if type_in_norm == type_out_norm:
        # No conversion needed. Return original data, but with normalized type.
        return type_out_norm, atomsposition_in

    # 4. Perform conversion
    atomsposition_out = {}
    
    if type_in_norm == 'Direct' and type_out_norm == 'Cartesian':
        # --- Convert from Direct -> Cartesian ---
        print("[INFO] Converting coordinates from Direct -> Cartesian...")
        for symbol, positions in atomsposition_in.items():
            atomsposition_out[symbol] = []
            for pos in positions:
                pos_np = np.array(pos, dtype=float)
                # v_cart = M . v_direct
                cart_pos = np.dot(lattice_matrix, pos_np)
              
                atomsposition_out[symbol].append(
                    [f"{cart_pos[0]:.8f}", f"{cart_pos[1]:.8f}", f"{cart_pos[2]:.8f}"]
                )
        return 'Cartesian', atomsposition_out

    elif type_in_norm == 'Cartesian' and type_out_norm == 'Direct':
        # --- Convert from Cartesian -> Direct ---
        print("[INFO] Converting coordinates from Cartesian -> Direct...")
        # v_direct = M_inv . v_cart
        inv_lattice = np.linalg.inv(lattice_matrix)
        for symbol, positions in atomsposition_in.items():
            atomsposition_out[symbol] = []
            for pos in positions:
                pos_np = np.array(pos, dtype=float)
                direct_pos = np.dot(inv_lattice, pos_np)
              
                atomsposition_out[symbol].append(
                    [f"{direct_pos[0]:.8f}", f"{direct_pos[1]:.8f}", f"{direct_pos[2]:.8f}"]
                )
        return 'Direct', atomsposition_out
    
    # Fallback (should not happen)
    return typevectors_in, atomsposition_in
##### END NEW HELPER FUNCTION #####

def main():

    parser = argparse.ArgumentParser(
        description="File format converter using stb-translate."
    )

    # Atualizado o texto de ajuda para incluir 'fdf'
    parser.add_argument("-if", "--in-format", required=True, choices=INPUT_FORMATS,
                        help="Input file format (options: poscar, cif, siesta, xyz, fhi, dftb, xsf, fdf)")
    parser.add_argument("-i", "--in-file", required=True,
                        help="Path to the input file")
    parser.add_argument("-of", "--out-format", required=True, choices=OUTPUT_FORMATS,
                        help="Output file format (options: cif , xyz, poscar, fdf, dftb, xsf, fhi)")
    parser.add_argument("-o", "--out-file", required=True,
                        help="Path to the output file")
    
    ##### NEW ARGUMENT #####
    parser.add_argument(
        "-cf", "--coord-format", 
        choices=['direct', 'cartesian'], 
        default=None,
        help="Specify the output coordinate format (direct/fractional or cartesian). "
             "If not specified, uses the input format or the output format's default."
    )
    ##### END NEW ARGUMENT #####

    parser.add_argument(
        "--lattice", help="Lattice vectors file, required only for XYZ output")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-translate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro == True:
        show_intro()

    print("\n" + color_text("TRANSLATE:", 'bold'))
    print("-"*60)

  
    if args.in_format == "xyz" and not args.lattice:
        parser.error(
            "The --lattice argument is required when input format is XYZ.")

    print(f"\n[INFO] Converting {args.in_file} ({args.in_format}) to {args.out_file} ({args.out_format})...")

  
    if args.coord_format:
        print(f"[INFO] Requested output coordinate format: {args.coord_format}")

    if args.out_format == "xyz":
        print(f"[INFO] Lattice vector file: {args.lattice}")

    # Adicionado o 'case' para o novo formato 'fdf'
    match (args.in_format):
        case "poscar":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_vasp(
                args.in_file)
        case "cif":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_cif(
                args.in_file)
        case "siesta":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_siesta(
                args.in_file)
        case "xyz":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xyz(
                args.in_file, args.lattice)
        case "fhi":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fhi(
                args.in_file)
        case "dftb":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_dftb(
                args.in_file)
        case "xsf":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_xsf(
                args.in_file)
        case "fdf":
            typevectors, latticeparameter, vectors, getatoms, atomsposition = getatomsandvectors_fdf(
                args.in_file)
        

    print(f"[OK] Read the file {args.in_file} ({args.in_format})")

    # All write functions now receive 'args.coord_format'
    match (args.out_format):
        case "xyz":
            writefilexyz(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)
        case "poscar":
            writefileposcar(typevectors, latticeparameter, vectors,
                            getatoms, atomsposition, args.out_file, args.coord_format)
        case "fdf":
            writefilefdf(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)
        case "dftb":
            writefiledftb(typevectors, latticeparameter, vectors,
                          getatoms, atomsposition, args.out_file, args.coord_format)
        case "xsf":
            writefilexsf(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)
        case "fhi":
            writefilefhi(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)

        case "cif":
            writefilecif(typevectors, latticeparameter, vectors,
                         getatoms, atomsposition, args.out_file, args.coord_format)

    print(f"[OK] Writing the file {args.out_file} ({args.out_format})")
    
    print("[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Converting input files is 10% coding, 90% crying.\n\n", 'bold'))

if __name__ == "__main__":
    main()
