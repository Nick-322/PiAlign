#!/usr/bin/env python3

# The master python script of PiAlign.
# Author: Yusaku Nitta
# E-mail: <ynitta6@gatech.edu>

import time
import argparse
import os
import subprocess
import sys
import numpy as np
import re
import itertools
import platform

sys.dont_write_bytecode = True
import filter_identical_interfaces
import find_large_interfaces

from collections import defaultdict
from prody import *
from prody.proteins import compare, localpdb
from filter_identical_interfaces import *
from find_large_interfaces import *


confProDy(verbosity='none')
confProDy(auto_secondary=True)

parser = argparse.ArgumentParser(
    description="PiAlign extracts interfaces from protein complexes and performs structural alignment of interface pairs.",
    epilog="""Example usage:
  PiAlign.py [options] -p1 <file1.pdb> -c1a <recptor_chainID(s)> -c1b <ligand_chainID(s)> -p2 <file2.pdb> -c2a <receptor_chainID(s)> -c2b <ligand_chainID(s)>
Practical case:
  PiAlign.py -transPDB -writeVMD -p1 7e5r.pdb -c1a HL -c1b A -p2 7zfa.pdb -c2a DE -c2b AB
"""
)

parser.add_argument('-f1', '--pdb1', dest='pdb_file1', help="First PDB file")
parser.add_argument('-f2', '--pdb2', dest='pdb_file2', help="Second PDB file")
parser.add_argument('-c1a', '--chain1a', dest='pchains1a', help="Protein chain(s) from the receptor side of an interface chain pair of PDB1")
parser.add_argument('-c1b', '--chain1b', dest='pchains1b', help="Protein chain(s) from the ligand side of an interface chain pair of PDB1")
parser.add_argument('-c2a', '--chain2a', dest='pchains2a', help="Protein chain(s) from the receptor side of an interface chain pair of PDB2")
parser.add_argument('-c2b', '--chain2b', dest='pchains2b', help="Protein chain(s) from the ligand side of an interface chain pair of PDB2")
parser.add_argument('-o', '--outdir', dest='output_directory', type=str, default='.', help="Output directory (default: current directory)")

# Advanced arguments group
advanced = parser.add_argument_group('Advanced options')
advanced.add_argument('-p', '--parsedPDB', dest='write_parsedPDB', action='store_true', help="Write a parsed PDB file")
advanced.add_argument('-t', '--transPDB', dest='trans_parsedPDB', action='store_true', help="Write PDB file 1 transformed to give the best interface alignment to PDB file 2")
advanced.add_argument('-w', '-writeVMD', dest='write_VMD', action='store_true', help="Write a VMD file of two interfaces that are superimposed")
advanced.add_argument('-s1','-searchIntCh1', dest='search_interface_chain1', action='store_true', help="For PDB1, perform a rigorous, strict search of interface chain combinations that capture a unique/large interface")
advanced.add_argument('-s2', '--searchIntCh2', dest='search_interface_chain2', action='store_true', help="For PDB2, perform a rigorous, strict search of interface chain combinations that capture a unique/large interface")

# Parse the arguments
args = parser.parse_args()

# Access the parsed arguments
pdb_file1 = args.pdb_file1
pdb_file2 = args.pdb_file2
pchains1a = args.pchains1a
pchains1b = args.pchains1b
pchains2a = args.pchains2a
pchains2b = args.pchains2b
write_parsedPDB = args.write_parsedPDB
trans_PDB = args.trans_parsedPDB
write_vmd = args.write_VMD
search_interface_chain1 = args.search_interface_chain1
search_interface_chain2 = args.search_interface_chain2
output_directory = args.output_directory

start = time.time()
if len(sys.argv) == 1:
    print('Error: invalid argument.')
    parser.print_help()
    sys.exit(1)

input_stem1 = str(os.path.basename(pdb_file1)).replace(".pdb", "").split('_')[0]
input_stem2 = str(os.path.basename(pdb_file2)).replace(".pdb", "").split('_')[0]

# 20 standard amino acids
standard_aa_names = ['ALA', 'CYS', 'ASP', 'GLU',
                     'PHE', 'GLY', 'HIS', 'ILE', 
                     'LYS', 'LEU', 'MET', 'ASN', 
                     'PRO', 'GLN', 'ARG', 'SER', 
                     'THR', 'VAL', 'TRP', 'TYR']

min_in_res = 20                             # Minimum number of protein residues of a interface
interface_distance_threshold = 4.5          # Distance threshold for interfaces
inter_interface_distance_threshold = 5.5    # Distance threshold for larger interfaces
seq_id_threshold = 95                       # 95% sequence identity as a threshold for identical protein chains


######################################################################
######################################################################
#########           Stage 1: parse raw PDB files          ############
######################################################################
######################################################################

def ExtractChainAtoms(chain_atoms, pchainsa, pchainsb):
    ''' Extract all atoms that belong to specified chains and returns their list. '''

    atom_list = []
    chain_ids = pchainsa + pchainsb 
    for chain_id in chain_ids:
        for atom in chain_atoms: 
            if atom.getChid() != chain_id: 
                continue
            if atom.getResname() not in standard_aa_names:
                continue
            if atom.getName() == 'OXT':
                continue
            if atom.getAltloc() not in ('A', ' '):  # There should be better way to filter (altloc can be other characters)
                continue
            atom_num = atom.getSerial()
            atom_name = atom.getName()
            altloc = atom.getAltloc()
            res_name = atom.getResname()
            chain_ID = atom.getChid()
            res_num = atom.getResnum()
            insertion_code = atom.getIcode()
            x, y, z = atom.getCoords()
            occupancy = atom.getOccupancy() 
            b_factor = atom.getBeta()
            entry = (
                f"{'ATOM  ':>4}"
                f"{atom_num:>5}"
                f"{' ':>2}"
                f"{atom_name:<3}"
                f"{altloc:>1}"
                f"{res_name:>3}"
                f"{' ':>1}"
                f"{chain_ID:>1}"
                f"{res_num:>4}"
                f"{insertion_code:>1}"
                f"{' ':>3}"
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                f"{occupancy:>6.2f}"
                f"{b_factor:>6.2f}\n"
            )
            atom_list.append(entry)
        atom_list.append("TER\n")
    return atom_list


def ParsePDB(pdb_file, pchainsa, pchainsb):
    ''' Parse PDB files and return atom records of target chains as AtomGroup class. '''

    start = time.time()
    structure = parsePDB(pdb_file)
    chain_atoms = structure.select(f"not hetatm and "
                                   f"(altloc A or altloc _) and "
                                   f"chain {' '.join(pchainsa+pchainsb)} and "
                                   f"not resname nonstdaa and "
                                   f"not name OXT")
    if write_parsedPDB:
        atom_list = ExtractChainAtoms(chain_atoms, pchainsa, pchainsb)   
        parsed_pdb = "".join(atom_list)
        input_stem = str(pdb_file).replace(".pdb", "").split("_")[0]
        output_file_name = f"{input_stem}{pchainsa}_{pchainsb}_parsed.pdb"
        output_path = os.path.join(output_directory, output_file_name)
        with open(output_path, 'w') as output_file:
            output_file.write(parsed_pdb)
        end = time.time()
        print(f"Time for creating a parsed pdb file: {end - start}s")             

    return chain_atoms


######################################################################
######################################################################
#########        Stage 2: extract interface residues        ##########
######################################################################
######################################################################

def FindIntAtoms(parsed_structure, pchainsa, pchainsb):
    ''' Extract atoms from the interface based on the distance threshold of 4.5Å. Returns a list of interface atom pairs. '''

    # Extracting atoms that belong to specified chains
    atoms_pchainsa = parsed_structure.select(f"chain {' '.join(pchainsa)}")
    atoms_pchainsb = parsed_structure.select(f"chain {' '.join(pchainsb)}")

    # Set the distance threshold for the interface (in Angstroms)
    distance_threshold = 4.5
    interface_atom_pairs = measure.contacts.findNeighbors(atoms_pchainsa, distance_threshold, atoms_pchainsb)
    if len(interface_atom_pairs) == 0:
        print(f'Error: No interface found between {pchainsa} and {pchainsb}.')
        sys.exit(1)

    return interface_atom_pairs


def ExtractIntAtoms(interface_atoms_set, seen_entry_set, additional_seen_entry_set, parsed_structure, Interface_list_set, chID):
    ''' Extract CA atom records from interface_atom_set and append it to a list. '''

    # DSSP secondary stucture assignment
    # 1->coil, 2->helix, 3->turn, 4->strand  
    secstr_mapping = {
    "C": 1,
    "S": 1,  
    "H": 2,  
    "G": 2,  
    "I": 2,  
    "T": 3,  
    "B": 4,  
    "E": 4
    }

    residue_index = interface_atoms_set.getResnum()
    insertion_code = interface_atoms_set.getIcode()
    chainID = interface_atoms_set.getChid()
    altloc = interface_atoms_set.getAltloc()

    # Filter1: for duplicated atoms, which have multiple atoms from the same residues
    # The other side of interface chain(s) within the specified distance
    if (residue_index, chainID, insertion_code, altloc) not in seen_entry_set:

        if altloc == " ":
            altloc = '_'
        if insertion_code == "":
            insertion_code = '_'

        CA_atoms = parsed_structure.select(f"name CA and "
                                            f"altloc {altloc} and "
                                            f"chain {chainID} and "
                                            f"resnum {residue_index}{insertion_code}")

        if CA_atoms is not None:
            for atom in CA_atoms:

                # Filter2: for selecting a record from only one atom from each interface residue
                if atom.getSerial() in additional_seen_entry_set:
                    continue
                atom_num = atom.getSerial()
                atom_name = atom.getName()
                altloc = atom.getAltloc()
                res_name = atom.getResname()
                if chID == "_":
                    chain_ID = atom.getChid()
                else:
                    chain_ID = chID
                res_num = atom.getResnum()
                insertion_code = atom.getIcode()
                x, y, z = atom.getCoords()
                secstr = secstr_mapping.get(atom.getSecstr(), 1)  # If Secstr assignment is unknown, 1 (blank space) is automatically assigned
                b_factor = atom.getBeta()
                entry = (
                    f"{'ATOM  ':>4}"
                    f"{atom_num:>5}"
                    f"{' ':>2}"
                    f"{atom_name:<3}"
                    f"{altloc:>1}"
                    f"{res_name:>3}"
                    f"{' ':>1}"
                    f"{chain_ID:>1}"
                    f"{res_num:>4}"
                    f"{insertion_code:>1}"
                    f"{' ':>3}"
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                    f"{secstr:>6.2f}"
                    f"{b_factor:>6.2f}\n"
                )
                seen_entry_set.add((residue_index, chainID, insertion_code, altloc))
                additional_seen_entry_set.add(atom.getSerial())
                Interface_list_set.append(entry)

            if chID != "_":
                return atom.getChid(), chain_ID, f"{res_num}{insertion_code}", res_name


def ExtractIntCA(interface_atom_pairs, parsed_structure, reachid1=None, reachid2=None):
    ''' Creates and returns a list of atom records of all interface CA atoms in a PDB format with no "TER" line(s). Counts and returns the total number of interface residues. '''

    Interface_list_set1 = []
    Interface_list_set2 = []
    interface_chain_pair = {}  # For VMD script generation
    chain_reassign_dict = {}
    chain_reassign_list_raw = []
    seen_entry_set1 = set()
    seen_entry_set2 = set()
    additional_seen_entry_set1 = set()
    additional_seen_entry_set2 = set()

    for pair in interface_atom_pairs:
        interface_atoms_set1 = pair[0] 
        interface_atoms_set2 = pair[1]
        interface_chain_pair[interface_atoms_set1.getChid()] = interface_atoms_set2.getChid()

        if reachid1 != '_':
            chid_dict = ExtractIntAtoms(interface_atoms_set1, seen_entry_set1, additional_seen_entry_set1, parsed_structure, Interface_list_set1, reachid1)
            chain_reassign_list_raw.append((chid_dict))
        else:
            ExtractIntAtoms(interface_atoms_set1, seen_entry_set1, additional_seen_entry_set1, parsed_structure, Interface_list_set1, reachid1)

        if reachid2 != '_':
            chid_dict2 = ExtractIntAtoms(interface_atoms_set2, seen_entry_set2, additional_seen_entry_set2, parsed_structure, Interface_list_set2, reachid2)
            chain_reassign_list_raw.append((chid_dict2))
        else:
            ExtractIntAtoms(interface_atoms_set2, seen_entry_set2, additional_seen_entry_set2, parsed_structure, Interface_list_set2, reachid2)

    Interface_list = Interface_list_set1 + Interface_list_set2
    Interface_list.sort()  # Sort the res entries by atom number
    num_cont_res = len(Interface_list)
    if num_cont_res < min_in_res:
        print("Error: the number of interface residues is too small. Please check your input files.")
        sys.exit(1)

    if len(chain_reassign_list_raw) != 0:
        chain_reassign_list = [element for element in chain_reassign_list_raw if element is not None]

        for item in chain_reassign_list:
            key = item[1:]
            value = item[0]
            if key not in chain_reassign_dict:
                chain_reassign_dict[key] = value 
            elif isinstance(chain_reassign_dict[key], list):
                chain_reassign_dict[key].append(value)  # If key already has a list, append the value
            else:
                chain_reassign_dict[key] = [chain_reassign_dict[key], value]  # If key has a single value, convert it to a list and append the new value

    if chain_reassign_dict:
        return Interface_list, num_cont_res, chain_reassign_dict, interface_chain_pair
    else:
        return Interface_list, num_cont_res, interface_chain_pair


def EditIntCA(Interface_list):
    ''' Add "TER" lines to the interface atom records list. '''

    Interface_list_edited = []
    current_chain = None

    for i in Interface_list:
        if (current_chain is not None and i[21] != current_chain):
            Interface_list_edited.append("TER\n")

        Interface_list_edited.append(i)
        current_chain = i[21]

    Interface_list_edited.append("TER\n")
    Int_CA_list = "".join(Interface_list_edited)

    return Int_CA_list


def WriteIntPDB(pdb_file, Int_CA_list, pchainsa, pchainsb):
    ''' Write a PDB file of interface atoms. '''

    full_path_pdb_file = os.path.abspath(pdb_file)
    input_stem = str(os.path.basename(full_path_pdb_file)).replace(".pdb", "").split('_')[0]
    output_file_name = f"{input_stem}{pchainsa}_{pchainsb}_int.pdb"
    output_path = os.path.join(output_directory, output_file_name)
    with open(output_path, 'w') as output_file:
        output_file.write(Int_CA_list)
    return output_file_name


######################################################################
######################################################################
#########      Stage 3: creating a contact residue list     ##########
######################################################################
######################################################################

def CountContact(interface_atom_pairs):
    ''' Counts and returns the total number of atomic contacts and residues contacts. Returns a list of residue-residue contact pairs with a new index assignment. '''

    residue_pair_1 = []
    residue_pair_2 = []

    for pair in interface_atom_pairs: 
        # Reassign indices of interface residues that each atom belongs to
        residue_pair_1.append(pair[0].getResindex())
        residue_pair_2.append(pair[1].getResindex())
    res_res_contact = list(zip(residue_pair_1, residue_pair_2))
    res_res_contact_list = set(res_res_contact)

    num_atomic_cont = len(interface_atom_pairs)
    num_res_cont = len(res_res_contact_list)

    return num_atomic_cont, num_res_cont, res_res_contact_list


def ReassignContRes(res_res_contact_list):
    ''' Create and returns residue-to-residue contact pairs with a new residue index assignment (from 1 in an ascending order) for the WriteConLst function. '''

    res_res_contact_sort = sorted(res_res_contact_list)
    # Assign new residue indices (from 1 in an ascending order) to each element of residue-to-residue contact pairs
    unique_values = sorted(set(value for pair in res_res_contact_sort for value in pair))

    # Create a mapping from each unique value to its assigned index
    value_to_index_mapping = {value: index + 1 for index, value in enumerate(unique_values)}

    # Update both elements of each pair using the mapping
    res_res_contact_updated = [(value_to_index_mapping[pair[0]], value_to_index_mapping[pair[1]]) for pair in res_res_contact_sort]

    # Consider residue-to-residue contacts from the other side of the pair
    flipped_res_res_contact = [(pair[1], pair[0]) for pair in res_res_contact_updated]
    flipped_res_res_contact_sort = sorted(flipped_res_res_contact)
    full_res_res_contact = res_res_contact_updated + flipped_res_res_contact_sort
    return full_res_res_contact


def WriteContLst(pdb_file, pchainsa, pchainsb, num_atomic_cont, num_cont_res, num_res_cont, res_res_contact_list):
    ''' Write and output a cont.lst file. '''

    receptor_chain = pchainsa
    ligand_chain = pchainsb

    full_path_pdb_file = os.path.abspath(pdb_file)
    input_stem = str(os.path.basename(full_path_pdb_file)).replace(".pdb", "").split('_')[0]
    output_name = f"{input_stem}{receptor_chain}_{ligand_chain}_con.lst"
    contact_distance_cutoff = 4.50
    secondary_structure_flag = 1
    full_res_res_contact = ReassignContRes(res_res_contact_list)

    output_path = os.path.join(output_directory, output_name)
    # Write the information to the lst file
    with open(output_path, 'w') as lst_file:
        lst_file.write("Atomic and residue contact lists\n\n")
        lst_file.write(f"Input PDB file                :: {pdb_file}\n")
        lst_file.write("Receptor Chain(s)             :: {}\n".format(receptor_chain))
        lst_file.write("Ligand Chain(s)               :: {}\n".format(ligand_chain))
        lst_file.write(f"Contact output                :: {output_name}\n")
        lst_file.write("Contact distance cutoff       :: {:<6.2f}\n".format(contact_distance_cutoff))
        lst_file.write("Secondary structure flag      :: {:<6d}\n".format(secondary_structure_flag))
        lst_file.write(f"Total atomic contacts found   :: {num_atomic_cont}\n")
        lst_file.write("Total interacting residues    :: {:<d}\n".format(num_cont_res))

        # Residue-residue contacts
        lst_file.write("\nTotal residue-residue contacts found  :: {:<d}\n".format(num_res_cont))
        lst_file.write(" ResIndex  NumCont Contact_Residues (starts 1 from the first interface residue)\n")

        # Create a dictionary to store pairs grouped by the first element
        grouped_pairs = {}

        # Iterate through the pairs and group them
        for pair in full_res_res_contact:
            key = pair[0]
            if key not in grouped_pairs:
                grouped_pairs[key] = []
            grouped_pairs[key].append(pair[1])

        # Print the formatted output
        for key in sorted(grouped_pairs.keys()):
            values = grouped_pairs[key]
            lst_file.write("RES  {:>4d}{:>4d}".format(key, len(values)))
            for value in values:
                lst_file.write("{:>5d}".format(value))
            lst_file.write("\n")
    return output_name


######################################################################
######################################################################
#########          Stage 4: Structural Alignment            ##########
######################################################################
######################################################################

## @@@@  use subroutines  @@@@@
# Call IS-align Fortran code to run structural alignment of interfaces 

def Int_Align(int_pdb_name1, int_pdb_name2, con_lst_name1, con_lst_name2):
    ''' Call and execute IS-align subroutines via subprocess. '''

    # Determine the path to the IS-align executable based on the OS
    system = platform.system()
    if system == 'Darwin':  # macOS
        is_align_executable = os.path.join(os.path.dirname(__file__), 'IS-align_mac')
    elif system == 'Linux':  # Linux
        is_align_executable = os.path.join(os.path.dirname(__file__), 'IS-align_linux')
    else:
        raise OSError(f"Unsupported operating system: {system}")

    # Define paths to the input files
    int_pdb_path1 = os.path.join(output_directory, int_pdb_name1)
    int_pdb_path2 = os.path.join(output_directory, int_pdb_name2)
    con_lst_path1 = os.path.join(output_directory, con_lst_name1)
    con_lst_path2 = os.path.join(output_directory, con_lst_name2)

    # Run the IS-align executable with the appropriate arguments
    result = subprocess.run([is_align_executable, '-v', '2', int_pdb_path1, int_pdb_path2, con_lst_path1, con_lst_path2], capture_output=True)
    stdout = result.stdout

    return stdout


######################################################################
######################################################################
#########               Parse IS-align Stdout               ##########
######################################################################
######################################################################
def ExtractVecMat(alnout):
    ''' Extract and return tansformation matrix and vector values from IS-align stdout. '''
    output = alnout.split('\n')
    # Remove unimportant information for end user
    while output and not output[0].startswith('Structure 1'):
        output.pop(0)

    trans_vec = []
    rot_mat = []
    fields = []
    flag = False
    count = 0
    for line in output:
        if 'Transformation matrix' in line:
            flag = True
            continue
        if flag and line.strip().split()[0].isdigit():
            fields.append(line.split())
            count += 1 
        if count == 3:
            break

    for i in range(3):
        trans_vec.append(float(fields[i][1]))
        rot_mat.append([float(i) for i in fields[i][2:5]])

    out = '\n'.join(output)

    return trans_vec, rot_mat, out


######################################################################
######################################################################
#########         Advanced option: Write a VMD file         ##########
######################################################################
######################################################################

def WriteVMD(pdb_file1, pchains1a, pchains1b, pdb_file2,  pchains2a, pchains2b, parsed_out, trans_vec, rot_mat, interface_chain_pair, interface_chain_pair2):
    ''' Generates a VMD file that visualizes an superimposed image of two interfaces. '''

    lines = parsed_out.split('\n')
    output_name = f"{input_stem1}{pchains1a}{pchains1b}_{input_stem2}{pchains2a}{pchains2b}.vmd"
    output_path = os.path.join(output_directory, output_name)

    full_path_pdb_file1 = os.path.abspath(pdb_file1)
    full_path_pdb_file2 = os.path.abspath(pdb_file2)

    vmd_rep = "NewCartoon 0.3 6.0 4.5 0"

    # Initialize interface residue lists
    int1 = [{'chain': '-'}]
    int2 = [{'chain': '-'}]
    flag = False

    # Find aligned interfacial residues
    for line in lines:
        if line.startswith(' Index Ch1 Resid1'):
            flag = True
            continue

        if flag and line.strip().split()[0].isdigit():
            ch1, res1, _, ch2, res2 = line.split()[1:6]

            if ch1 != int1[-1]['chain']:
                int1.append({'chain': ch1, 'lst': []})
            if ch2 != int2[-1]['chain']:
                int2.append({'chain': ch2, 'lst': []})

            int1[-1]['lst'].append(res1)
            int2[-1]['lst'].append(res2)

        if line.startswith(' ":" ') or line.startswith(' "*" '):
            break

    # Extract chain and residue information
    int1 = int1[1:] 
    int2 = int2[1:] 
    chain1_list = {}
    aln1_list = {}
    for i in range(1, len(int1)):
        chain1_list['ch1_'+str(i)] = int1[i]['chain']
        aln1_list['aln1_'+str(i)] = ' '.join(int1[i]['lst'])

    chain2_list = {}
    aln2_list = {}
    for i in range(1, len(int2)):
        chain2_list['ch2_'+str(i)] = int2[i]['chain']
        aln2_list['aln2_'+str(i)] = ' '.join(int2[i]['lst'])

    interface_chain_pair_swapped = {}
    for key, value in interface_chain_pair.items():
        if value in interface_chain_pair_swapped:
            interface_chain_pair_swapped[value].append(key)
        else:
            interface_chain_pair_swapped[value] = [key]

    interface_chain_pair_swapped2 = {}
    for key, value in interface_chain_pair2.items():
        if value in interface_chain_pair_swapped2:
            interface_chain_pair_swapped2[value].append(key)
        else:
            interface_chain_pair_swapped2[value] = [key]

    with open(output_path, 'w') as lst_file:
        lst_file.write("#!/usr/local/bin/vmd\n")
        lst_file.write("### VMD script for visualizing the alignment generated by PiAlign\n\n")

    # Script for loading the first molecule
        lst_file.write(f"mol new {full_path_pdb_file1} type pdb\n")
        lst_file.write("mol delrep 0 top\n")

        for i in int1:
            if i['chain'] in list(pchains1a):
                chain_id1b = interface_chain_pair[i['chain']]
                lst_file.write(f"mol representation {vmd_rep}\n")
                lst_file.write("mol color ColorID 10\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and same residue as protein within 4.5 of chain {chain_id1b}}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")

                lst_file.write(f"mol representation {vmd_rep}\n")
                lst_file.write("mol color ColorID 10\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and not same residue as protein within 4.5 of chain {chain_id1b}}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")
            elif i['chain'] in list(pchains1b):
                chain_id1a = interface_chain_pair_swapped[i['chain']]
                for chaina in chain_id1a:
                    lst_file.write(f"mol representation {vmd_rep}\n")
                    lst_file.write("mol color ColorID 0\n")
                    lst_file.write(f"mol selection {{chain {i['chain']} and same residue as protein within 4.5 of chain {chaina}}}\n")
                    lst_file.write("mol material BrushedMetal\n")
                    lst_file.write("mol addrep top\n")

                    lst_file.write(f"mol representation {vmd_rep}\n")
                    lst_file.write("mol color ColorID 0\n")
                    lst_file.write(f"mol selection {{chain {i['chain']} and not same residue as protein within 4.5 of chain {chaina}}}\n")
                    lst_file.write("mol material BrushedMetal\n")
                    lst_file.write("mol addrep top\n")
        for i in int1:
            if i['chain'] in list(pchains1a):
                #### show aligned Ca in VDW rep
                lst_file.write("mol representation VDW 0.6 8.0\n")
                lst_file.write("mol color ColorID 10\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and resid {' '.join(i['lst'])} and name CA}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")
            elif i['chain'] in list(pchains1b):
                lst_file.write("mol representation VDW 0.6 8.0\n")
                lst_file.write("mol color ColorID 0\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and resid {' '.join(i['lst'])} and name CA}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")

        ##### script for transformation
        lst_file.write("\nset transMat {}\n")
        lst_file.write(f"lappend transMat {{{' '.join([str(rot_mat[0][0]), str(rot_mat[0][1]), str(rot_mat[0][2]), str(trans_vec[0])])}}}\n")
        lst_file.write(f"lappend transMat {{{' '.join([str(rot_mat[1][0]), str(rot_mat[1][1]), str(rot_mat[1][2]), str(trans_vec[1])])}}}\n")
        lst_file.write(f"lappend transMat {{{' '.join([str(rot_mat[2][0]), str(rot_mat[2][1]), str(rot_mat[2][2]), str(trans_vec[2])])}}}\n")
        lst_file.write("lappend transMat {0.0 0.0 0.0 1.0}\n")
        lst_file.write("set myMol [atomselect top all]\n")
        lst_file.write("$myMol move $transMat\n\n")

        ##### script for loading the second molecule
        lst_file.write(f"mol new {full_path_pdb_file2} type pdb\n")
        lst_file.write("mol delrep 0 top\n")

        for i in int2:
            if i['chain'] in list(pchains2a):
                chain_id2b = interface_chain_pair2[i['chain']]
                lst_file.write(f"mol representation {vmd_rep}\n")
                lst_file.write("mol color ColorID 3\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and same residue as protein within 4.5 of chain {chain_id2b}}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")

                lst_file.write(f"mol representation {vmd_rep}\n")
                lst_file.write("mol color ColorID 3\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and not same residue as protein within 4.5 of chain {chain_id2b}}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")

            elif i['chain'] in list(pchains2b):
                chain_id2a = interface_chain_pair_swapped2[i['chain']]
                for chaina in chain_id2a:
                    lst_file.write(f"mol representation {vmd_rep}\n")
                    lst_file.write("mol color ColorID 1\n")
                    lst_file.write(f"mol selection {{chain {i['chain']} and same residue as protein within 4.5 of chain {chaina}}}\n")
                    lst_file.write("mol material BrushedMetal\n")
                    lst_file.write("mol addrep top\n")

                    lst_file.write(f"mol representation {vmd_rep}\n")
                    lst_file.write("mol color ColorID 1\n")
                    lst_file.write(f"mol selection {{chain {i['chain']} and not same residue as protein within 4.5 of chain {chaina}}}\n")
                    lst_file.write("mol material BrushedMetal\n")
                    lst_file.write("mol addrep top\n")

        for i in int2:
            if i['chain'] in list(pchains2a):
                #### show aligned Ca in VDW rep
                lst_file.write("mol representation VDW 0.6 8.0\n")
                lst_file.write("mol color ColorID 3\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and resid {' '.join(i['lst'])} and name CA}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")
            elif i['chain'] in list(pchains2b):
                lst_file.write("mol representation VDW 0.6 8.0\n")
                lst_file.write("mol color ColorID 1\n")
                lst_file.write(f"mol selection {{chain {i['chain']} and resid {' '.join(i['lst'])} and name CA}}\n")
                lst_file.write("mol material BrushedMetal\n")
                lst_file.write("mol addrep top\n")

    return output_name


######################################################################
######################################################################
#########       Advanced option: Transform  PDB file        ##########
######################################################################
######################################################################

def transform_pdb(pdb_file1, pchains1a, pchains1b, parsed_structure, rot_mat, trans_vec):
    ''' Generates a PDB file of the 1st complex with modified cartesian coordinates based on transform matrix and vector values. '''

    input_stem = str(pdb_file1).replace(".pdb", "").split('_')[0]
    output_name = f"{input_stem}_trans.pdb"
    output_path = os.path.join(output_directory, output_name)

    # Write the information to the lst file
    atom_list = []
    chain_ids = pchains1a + pchains1b 
    for chain_id in chain_ids: 
        for atom in parsed_structure: 
            if atom.getChid() != chain_id: 
                continue
            atom_num = atom.getSerial()
            atom_name = atom.getName()
            altloc = atom.getAltloc()
            res_name = atom.getResname()
            chain_ID = atom.getChid()
            res_num = atom.getResnum()
            insertion_code = atom.getIcode()
            x, y, z = atom.getCoords()
            tx = trans_vec[0] + x * rot_mat[0][0] + y * rot_mat[0][1] + z * rot_mat[0][2]
            ty = trans_vec[1] + x * rot_mat[1][0] + y * rot_mat[1][1] + z * rot_mat[1][2]
            tz = trans_vec[2] + x * rot_mat[2][0] + y * rot_mat[2][1] + z * rot_mat[2][2]

            occupancy = atom.getOccupancy() 
            b_factor = atom.getBeta()
            entry = (
                f"{'ATOM  ':>4}"
                f"{atom_num:>5}"
                f"{' ':>2}"
                f"{atom_name:<3}"
                f"{altloc:>1}"
                f"{res_name:>3}"
                f"{' ':>1}"
                f"{chain_ID:>1}"
                f"{res_num:>4}"
                f"{insertion_code:>1}"
                f"{' ':>3}"
                f"{tx:>8.3f}{ty:>8.3f}{tz:>8.3f}"
                f"{occupancy:>6.2f}"
                f"{b_factor:>6.2f}\n"
            )
            atom_list.append(entry)
        atom_list.append("TER\n")
    trans_pdb_list = "".join(atom_list)
    with open(output_path, 'w') as output_file:
        output_file.write(trans_pdb_list)
    return output_name


#####################################################################################
#####################################################################################
###  Advanced option: filtering identical interface and identify super interface  ###
#####################################################################################
#####################################################################################

def filter_iden_int_and_find_large_int(input_stem_name, pdb_file, pchainsa, pchainsb):
    print(f" Filtering identical interface ChainIDs and arranging chainID assignment for larger interfaces from {input_stem_name}.pdb... ")
    unique_int_chain_dict = {}  # Before final_int_chain_dict dictionary. This dictionary does not consider multi-chain interface
    chain_atoms = create_chain_atom_obj(pdb_file, pchainsa, pchainsb)

    chain_selection_obj_dict = {}
    for chain in pchainsa+pchainsb:
        chain_selection_obj_dict[chain] = chain_atoms.select(f"chain {chain}")


    ### 1. All-against-all pair-wise neighbor search betweeen receptors chain ###
    '''This block of code performs all-against-all interface pair search returns a list of interface chainID pairs'''
    ### Why do we do this?: to correctly assign less significantly interfacing chains (ex. light chains from antibodies)
    ## a. Find interacting receptor pairs.
    int_receptor_pair = all_against_all_int_search(chain_selection_obj_dict, pchainsa)


    ## b. Check if there are any un-paired receptor (ex. nanobody)
    orphan_receptor_list = find_orphan_receptor(int_receptor_pair, pchainsa)


    ### 2. Identifying interacting ligand-receptors sets
    ### Create a dictionary of chains that are interacting with each of antigen/ligand chains (key: antigen/ligand chainIDs, value: antibody/receptor chainIDs).
    ## a. Initiate a dictionary that contains ligand chainID as key and its interacting receptor chainIDs as value. There could be duplicate (miss assignment of receptor chains).
    int_chain_dict = create_init_ligand_receptor_dict(chain_selection_obj_dict, int_receptor_pair, pchainsb, min_in_res, interface_distance_threshold)

    if len(orphan_receptor_list) > 0:
        assign_orphan_receptor(orphan_receptor_list, chain_selection_obj_dict, pchainsb, interface_distance_threshold, int_chain_dict)


    ## b. Identify duplicates (miss assignment of receptor chains).
    duplicates = find_duplicates(int_chain_dict)

    if len(pchainsa) == 1:
        duplicates = None

    ## c. Filter out duplicate by calculating the ratio of much the atomic interaction happens over all the possible atomic interaction between duplicate chain and all the rest chain assigned in the same groups.
    if duplicates is not None and len(duplicates) != 0:
        filter_duplicate_chain(duplicates, int_chain_dict, chain_selection_obj_dict, int_receptor_pair)


    ### 3. perform the sequence alignments of chains to identify/filter out identical ligand and receptor chains.
    ### Perform the sequence alignments of chains to identify/filter out identical ligand and receptor chains.
    ### What does step return?: unique_int_chain_dict... unique interface ligand-receptor chain sets
    ## a. Check the sequence identity of ligand chains.
    # i. Create a dictionary that stores ligand chainIDs as keys and selection objects (used for alignment) as values.
    ligand_chain_res_dict = {}  # KEEP!!
    identical_ligand_list = []  # KEEP!!

    ligand_chain_res_dict = create_ligand_chain_res_dict(int_chain_dict, chain_selection_obj_dict)  # keys: ligand chainID, values: selection objects (used for alignment)


    # ii. Perform sequence identity check of ligand chains if there's more than one ligand chains in PDB.
    if len(ligand_chain_res_dict) > 1:
        identical_ligand_list = seq_id_check_ligand(chain_selection_obj_dict, int_chain_dict, seq_id_threshold)


    ## b. Check the sequence identity of receptor chains.
    # Find identical receptor chains among different interface receptor-ligand sets. Create a disctionary of identical receptor sets (key: one chain from identical receptor chain group, value: all the rest of chains from identical receptor chain group.
    if len(int_chain_dict) > 1:
        identical_receptor_dict = seq_id_check_receptor(chain_selection_obj_dict, int_chain_dict, seq_id_threshold)

    else:
        identical_receptor_dict = {}


    ## c. Identify unique interface chain sets.
    # i. Generate a nested list: each sublist contains a set of identical receptor chains (ex. [[receptor1, receptor2], [...]...]
    identical_receptor_comb_list = [[receptor_prime] + receptors_secondary for receptor_prime, receptors_secondary in identical_receptor_dict.items()]

    if len(identical_receptor_comb_list) != 0:
        identical_interaction_set_dict = map_identical_interactions(int_chain_dict, identical_receptor_comb_list)

        if len(identical_interaction_set_dict) > 1:
            # Handle the edge case for having differnet number of identical interaction identical interaction set
            identical_interaction_set_dict = reorder_identical_interaction_set_dict(identical_interaction_set_dict)
            refine_identical_interaction_set_dict(identical_interaction_set_dict)


        # ii. For each set of identical interfaces, find a single interface that has the most number of interfacing chains'
        best_unique_int_dict = find_best_interfaces(identical_interaction_set_dict, chain_selection_obj_dict)


        # iii. Keep track of identical receptor chains with a "non_unique_receptors" list.
        non_unique_receptors = track_non_unique_receptors(identical_receptor_dict, best_unique_int_dict)


        # iv. Use int_chain_dict to keep one of the identical chains and filter out the rest of identical chains to create a dictionary of unique ligand-receptor pairs
        # unique_int_chain_dict will have a ligand as key and its interacting receptors + the receptors' partner(s) as values, but not all partners necesarily have significant interaction with ligand
        # Handle edge case where orphan receptor-ligand is involved in the identical interactions
        filter_identical_interface(int_chain_dict, identical_ligand_list, unique_int_chain_dict, non_unique_receptors, identical_interaction_set_dict)


    # When int_chain_dict already capture all the unique interfaces
    if len(unique_int_chain_dict) == 0:
        unique_int_chain_dict = int_chain_dict


    ### 4. Super interface detection
    ## Calculate the pairwise distance of residues in two interfaces to check if there are any super interface (involve more than 1) (threshold: 5.5)
    super_int_chain_dict = find_large_interface(unique_int_chain_dict, chain_atoms, chain_selection_obj_dict)

    if len(super_int_chain_dict) > 1:
        super_int_chain_dict = find_large_interface_across_ligand(super_int_chain_dict, chain_atoms)


    ### 5. Final check: confirm the number of interfacing residues is greater than or equal to 20. Finalize the interface ligand-receptor sets
    valid_super_int_chain_dict = find_valid_super_int(chain_selection_obj_dict, super_int_chain_dict)

    if len(valid_super_int_chain_dict) == 0:
        print(f'Error: No valid interface found.')
        sys.exit(1)

    print(f" Done! \n {{ligand chain: [receptor chains]}}... {valid_super_int_chain_dict}")

    return valid_super_int_chain_dict


###################################################################################
###################################################################################
#########   executing functions for interface extraction and alignment   ##########
###################################################################################
###################################################################################

def interface_extraction_and_alignment(pdb_file1, pchains1a, pchains1b, pdb_file2, pchains2a, pchains2b):
    print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("+++++              Step 1: Processing PDB files ...              +++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

    parsed_structure = ParsePDB(pdb_file1, pchains1a, pchains1b)
    parsed_structure2 = ParsePDB(pdb_file2, pchains2a, pchains2b)

    print(f"{pdb_file1}:  found {len(pchains1a+pchains1b)} protein chains {pchains1a} {pchains1b}")
    print(f"{pdb_file2}:  found {len(pchains2a+pchains2b)} protein chains {pchains2a} {pchains2b}\n")



    print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("+++++      Step 2: Extracting protein-protein interfaces ...     +++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    print("Using only Calpha atoms for defining protein-protein interfaces.\n")
    print(f"Minimum number of interface residue required: 20 AAs")

    # 1st PDB file
    interface_atom_pairs = FindIntAtoms(parsed_structure, pchains1a, pchains1b)
    total_atomic_cont, total_res_cont, res_res_contact_list = CountContact(interface_atom_pairs)

    if len(pchains1a) > 1:
        if len(pchains1b) > 1:
            Interface_list, total_cont_res, chain_reassign_dict, interface_chain_pair = ExtractIntCA(interface_atom_pairs, parsed_structure, '@', '#') 
        elif len(pchains1b) == 1:
            Interface_list, total_cont_res, chain_reassign_dict, interface_chain_pair = ExtractIntCA(interface_atom_pairs, parsed_structure, '@', '_') 
    elif len(pchains1a) == 1: 
        if len(pchains1b) > 1:
            Interface_list, total_cont_res, chain_reassign_dict, interface_chain_pair = ExtractIntCA(interface_atom_pairs, parsed_structure, '_', '#') 
        elif len(pchains1b) == 1:
            Interface_list, total_cont_res, interface_chain_pair = ExtractIntCA(interface_atom_pairs, parsed_structure, '_', '_')

    Int_CA_list = EditIntCA(Interface_list)
    int_pdb_name1 = WriteIntPDB(pdb_file1, Int_CA_list, pchains1a, pchains1b)
    con_lst_name1 = WriteContLst(pdb_file1, pchains1a, pchains1b, total_atomic_cont, total_cont_res, total_res_cont, res_res_contact_list)
    print(f"{pdb_file1}: found 1 valid PPI(s) {''.join(pchains1a+pchains1b)} {total_cont_res}")


    # 2nd PDB file
    interface_atom_pairs2 = FindIntAtoms(parsed_structure2, pchains2a, pchains2b)
    total_atomic_cont2, total_res_cont2, res_res_contact_list2 = CountContact(interface_atom_pairs2)

    if len(pchains2a) > 1:
        if len(pchains2b) > 1:
            Interface_list2, total_cont_res2, chain_reassign_dict2, interface_chain_pair2 = ExtractIntCA(interface_atom_pairs2, parsed_structure2, '$', '%') 
        elif len(pchains2b) == 1:
            Interface_list2, total_cont_res2, chain_reassign_dict2, interface_chain_pair2 = ExtractIntCA(interface_atom_pairs2, parsed_structure2, '$', '_') 
    elif len(pchains2a) == 1: 
        if len(pchains2b) > 1:
            Interface_list2, total_cont_res2, chain_reassign_dict2, interface_chain_pair2 = ExtractIntCA(interface_atom_pairs2, parsed_structure2, '_', '%') 
        elif len(pchains2b) == 1:
            Interface_list2, total_cont_res2, interface_chain_pair2 = ExtractIntCA(interface_atom_pairs2, parsed_structure2, '_', '_')
    Int_face_PDB2 = EditIntCA(Interface_list2)
    int_pdb_name2 = WriteIntPDB(pdb_file2, Int_face_PDB2, pchains2a, pchains2b)
    con_lst_name2 = WriteContLst(pdb_file2, pchains2a, pchains2b, total_atomic_cont2, total_cont_res2, total_res_cont2, res_res_contact_list2)
    print(f"{pdb_file2}: found 1 valid PPI(s) {''.join(pchains2a+pchains2b)} {total_cont_res2}\n")



    print("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("+++++ Step 3: Performing protein-protein interface alignment ... +++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    print(f">>>{input_stem1}{pchains1a}_{pchains1b} vs {input_stem2}{pchains2a}_{pchains2b}\n")

    out = Int_Align(int_pdb_name1, int_pdb_name2, con_lst_name1, con_lst_name2)
    alnout = out.decode('utf-8')
    alnout_lines = alnout.splitlines(keepends=True)
    num_lines = len(alnout_lines)

    # Remove unimportant information for end user
    for i in range(num_lines):
        if alnout_lines[0].startswith(' Error'):
            print("Error: the number of contacts is too small. Please check your input files.")
            sys.exit(1)
        elif alnout_lines[0].startswith('Structure 1'):
            break
        else:
            alnout_lines.pop(0)

    num_lines_update = len(alnout_lines)


    start_modification = False

    # Perform the chainID reassignment when any side of interface has more than 1 chain 
    if any(len(chain1) > 1 for chain1 in [pchains1a, pchains1b]) or any(len(chain2) > 1 for chain2 in [pchains2a, pchains2b]):
        for i in range(num_lines_update):
            if alnout_lines[i].startswith('Structure 1:'):
                if len(pchains1a) > 1 and len(pchains1b) > 1:
                    alnout_lines[i] = alnout_lines[i].replace("@ #", f"{pchains1a} {pchains1b}")
                    # In case where the masked chain order was swapped 
                    alnout_lines[i] = alnout_lines[i].replace("# @", f"{pchains1b} {pchains1a}")

                elif len(pchains1a) > 1 and len(pchains1b) == 1:
                    alnout_lines[i] = alnout_lines[i].replace("@", f"{pchains1a}")
                elif len(pchains1a) == 1 and len(pchains1b) > 1:
                    alnout_lines[i] = alnout_lines[i].replace("#", f"{pchains1b}")

            elif alnout_lines[i].startswith('Structure 2:'):
                if len(pchains2a) > 1 and len(pchains2b) > 1:
                    alnout_lines[i] = alnout_lines[i].replace("$ %", f"{pchains2a} {pchains2b}") 
                    alnout_lines[i] = alnout_lines[i].replace("% $", f"{pchains2b} {pchains2a}") 
                elif len(pchains2a) > 1 and len(pchains2b) == 1:
                    alnout_lines[i] = alnout_lines[i].replace("$", f"{pchains2a}") 
                elif len(pchains2a) == 1 and len(pchains2b) > 1:
                    alnout_lines[i] = alnout_lines[i].replace("%", f"{pchains2b}") 


            if alnout_lines[i].startswith(" Index"):
                start_modification = True
                continue 

            if alnout_lines[i].startswith(' ":"'):
                start_modification = False
                break

            elif start_modification:
                chain1, Resid1, AA1, chain2, Resid2, AA2 = alnout_lines[i].strip().split()[1:7]
                if chain1 == '@' or chain1 == '#': 
                    current_chain1 = chain_reassign_dict[(chain1, Resid1, AA1)]
                    current_resid1 = Resid1

                    # Check if multiple original chainIDs are called from the same set of temporary ChiD, residue number, and amino acid type calls 
                    if len(chain_reassign_dict[(chain1, Resid1, AA1)]) == 1:  # If there's no duplicate assinments for int residues 
                        if isinstance(chain_reassign_dict[(chain1, Resid1, AA1)], list):  # This is required for case where original chainIDs that idenetical set of tempID. resnum, AA but the chainIDs were removing from the list to one. 
                            modified_line = alnout_lines[i].replace(chain1, chain_reassign_dict[(chain1, Resid1, AA1)][0])
                            current_chain1 = chain_reassign_dict[(chain1, Resid1, AA1)][0]

                        else:  # If there's no duplicate assinments for int residues 
                            modified_line = alnout_lines[i].replace(chain1, chain_reassign_dict[(chain1, Resid1, AA1)])
                            current_chain1 = chain_reassign_dict[(chain1, Resid1, AA1)]
                            current_resid1 = Resid1

                    elif len(chain_reassign_dict[(chain1, Resid1, AA1)]) > 1:  # Avoid the mis assignment of ChID due to duplicate original Chain IDs
                        chain1pv, Resid1pv = alnout_lines[i-1].strip().split()[1:3]
                        current_chain1 = chain1pv
                        current_resid1 = Resid1pv
                        if current_chain1 == 'Ch1' and current_resid1 == 'Resid1':
                            current_chain1 = chain_reassign_dict[(chain1, Resid1, AA1)]
                            current_resid1 = Resid1
                        for chain_1 in chain_reassign_dict[(chain1, Resid1, AA1)]:
                            Resid1_int = int(re.match(r'\d+', Resid1).group())
                            current_resid1_int = int(re.match(r'\d+', current_resid1).group())
                            if chain_1 == current_chain1 and Resid1_int > current_resid1_int:  # If there's a preceeding residue entry that has the same chainid (not superchain) as the current residue
                                chain1_var = alnout_lines[i].strip().split()[1]
                                modified_line = alnout_lines[i].replace(chain1_var, chain_1)  
                                current_chain1 = chain_1
                                current_resid1 = Resid1
                                chain_reassign_dict[(chain1, Resid1, AA1)].remove(chain_1)
                                break


                            else:  # Back tracking if there's no previous residue entry that has the same chainid (not superchain) as the current residue to refer 
                                if alnout_lines[i+1].startswith(' ":"'):  # If the residue is at the end of the int list
                                    modified_line = alnout_lines[i].replace(chain1, chain_1)

                                else:  # If the residue is not at the end of the int list --> res is in the middle, but previous res has a diff ChID
                                    chain1nx, Resid1nx, AA1nx = alnout_lines[i+1].strip().split()[1:4]
                                    try:
                                        chain_reassign_dict[(chain1nx, Resid1nx, AA1nx)]
                                        chainnx1 = chain_reassign_dict[(chain1nx, Resid1nx, AA1nx)]
                                        Resid1nx_int = int(re.match(r'\d+', Resid1nx).group())

                                        if chain_1 == chainnx1 and Resid1nx_int > Resid1_int: 
                                            chain1_var = alnout_lines[i].strip().split()[1]
                                            modified_line = alnout_lines[i].replace(chain1_var, chain_1)
                                            current_chain1 = chain_1
                                            current_resid1 = Resid1nx
                                            chain_reassign_dict[(chain1, Resid1, AA1)].remove(chain_1)
                                            break

                                        else:  # Ex. if chainid, AA, and resid are identical within a superchain 
                                            modified_line = alnout_lines[i].replace(chain1, chain_1)
                                    except KeyError: 
                                        if chain_1 == current_chain1 and Resid1_int <= current_resid1_int:
                                            continue
                                        modified_line = alnout_lines[i].replace(chain1, chain_1)        

                else: 
                    modified_line = alnout_lines[i]

                if chain2 == '$' or chain2 == '%':
                    current_chain2 = chain_reassign_dict2[(chain2, Resid2, AA2)]
                    current_resid2 = Resid2

                    # Check if multiple original chainIDs are called from the same set of temporary ChiD, residue number, and amino acid type calls 
                    if len(chain_reassign_dict2[(chain2, Resid2, AA2)]) == 1:  # If there's no duplicate assinments for int residues
                        if isinstance(chain_reassign_dict2[(chain2, Resid2, AA2)], list):
                            modified_line = modified_line.replace(chain2, chain_reassign_dict2[(chain2, Resid2, AA2)][0])
                            current_chain2 = chain_reassign_dict2[(chain2, Resid2, AA2)][0]

                        elif not isinstance(chain_reassign_dict2[(chain2, Resid2, AA2)], list):
                            modified_line = modified_line.replace(chain2, chain_reassign_dict2[(chain2, Resid2, AA2)])
                            current_chain2 = chain_reassign_dict2[(chain2, Resid2, AA2)]
                            current_resid2 = Resid2

                    elif len(chain_reassign_dict2[(chain2, Resid2, AA2)]) > 1:  # Avoid the mis assignment of ChID due to duplicate original Chain IDs
                        chain2pv, Resid2pv = alnout_lines[i-1].strip().split()[4:6]
                        current_chain2 = chain2pv
                        current_resid2 = Resid2pv
                        if current_chain2 == 'Ch2' and current_resid2 == 'Resid2':
                            current_chain2 = chain_reassign_dict2[(chain2, Resid2, AA2)]
                            current_resid2 = Resid2
                        for chain_2 in chain_reassign_dict2[(chain2, Resid2, AA2)]:
                            Resid2_int = int(re.match(r'\d+', Resid2).group())
                            current_resid2_int = int(re.match(r'\d+', current_resid2).group())
                            if chain_2 == current_chain2 and Resid2_int > current_resid2_int:  # If there's a preceeding residue entry that has the same chainid (not superchain) as the current residue
                                chain2_var = modified_line.strip().split()[4]
                                modified_line = modified_line.replace(chain2_var, chain_2)
                                current_chain2 = chain_2
                                current_resid2 = Resid2
                                chain_reassign_dict2[(chain2, Resid2, AA2)].remove(chain_2)
                                break


                            else:  # Back tracking if there's no previous residue entry that has the same chainid (not superchain) as the current residue to refer 
                                if alnout_lines[i+1].startswith(' ":"'):  # If the residue is at the end of the int list
                                    modified_line = modified_line.replace(chain2, chain_2)

                                else:  # If the residue is not at the end of the int list --> res is in the middle, but previous res has a diff ChID
                                    chain2nx, Resid2nx, AA2nx = alnout_lines[i+1].strip().split()[4:7]
                                    try:
                                        chain_reassign_dict2[(chain2nx, Resid2nx, AA2nx)]
                                        chainnx2 = chain_reassign_dict2[(chain2nx, Resid2nx, AA2nx)]
                                        Resid2nx_int = int(re.match(r'\d+', Resid2nx).group())

                                        if chain_2 == chainnx2 and Resid2nx_int > Resid2_int:  # If a residue after has the same ChID
                                            chain2_var = modified_line.strip().split()[4]
                                            modified_line = modified_line.replace(chain2_var, chain_2)
                                            current_chain2 = chain_2
                                            current_resid2 = Resid2
                                            chain_reassign_dict2[(chain2, Resid2, AA2)].remove(chain_2)
                                            break

                                        else:  # Ex. if chainid, AA, and resid are identical within a superchain 
                                            modified_line = modified_line.replace(chain2, chain_2)

                                    except KeyError:  # If chain2nx, Resid2nx, AA2nx is beyond the scope of chain_reassign_dict2
                                        if chain_2 == current_chain2 and Resid2_int <= current_resid2_int:
                                            continue
                                        modified_line = modified_line.replace(chain2, chain_2)

                alnout_lines[i] = modified_line  # *****DON'T CHANGE THE POSITION*******  

    print(''.join(alnout_lines))  # *****DON'T CHANGE THE POSITION******* 

    if trans_PDB:
        parsed_out = ''.join(alnout_lines)
        trans_vec, rot_mat, out = ExtractVecMat(parsed_out)
        output_name = transform_pdb(pdb_file1, pchains1a, pchains1b, parsed_structure, rot_mat, trans_vec)
        print(f"Transformed info was saved to {output_name}\n")

    if write_vmd:
        parsed_out = ''.join(alnout_lines)
        trans_vec, rot_mat, out = ExtractVecMat(parsed_out)
        out_name = WriteVMD(pdb_file1, pchains1a, pchains1b, pdb_file2,  pchains2a, pchains2b, out, trans_vec, rot_mat, interface_chain_pair, interface_chain_pair2)
        print(f"VMD script was saved to {out_name}\n")


if pdb_file1 and pdb_file2:
    # Fetch a PDB file if the file is not in a current working directory
    if not os.path.exists(pdb_file1):
        pdb_id = str(pdb_file1).replace(".pdb", "")
        localpdb.fetchPDB(pdb_id, compressed=False)
        print(f"––Since {pdb_file1} was not found in your current working directory, the PDB file is fetched now.")

    # Fetch a PDB file if the file is not in a current working directory
    if not os.path.exists(pdb_file2):
        pdb_id2 = str(pdb_file2).replace(".pdb", "")
        localpdb.fetchPDB(pdb_id2, compressed=False)
        print(f"––Since {pdb_file2} was not found in your current working directory, the PDB file is fetched now.")

    if search_interface_chain1 and not search_interface_chain2:
        receptor_ligand_pair = []
        valid_super_int_chain_dict1 = filter_iden_int_and_find_large_int(input_stem1, pdb_file1, pchains1a, pchains1b)
        
        for ligands, receptors in valid_super_int_chain_dict1.items():
            receptor_ligand_pair_per_key = list(itertools.product(receptors, [ligands]))
            receptor_ligand_pair.extend(receptor_ligand_pair_per_key)

        for receptor, ligand in receptor_ligand_pair:
            interface_extraction_and_alignment(pdb_file1, receptor, ligand, pdb_file2, pchains2a, pchains2b)
            
    elif search_interface_chain2 and not search_interface_chain1:
        receptor_ligand_pair = []
        valid_super_int_chain_dict2 = filter_iden_int_and_find_large_int(input_stem2, pdb_file2, pchains2a, pchains2b)

        for ligands, receptors in valid_super_int_chain_dict2.items():
            receptor_ligand_pair_per_key = list(itertools.product(receptors, [ligands]))
            receptor_ligand_pair.extend(receptor_ligand_pair_per_key)

        for receptor, ligand in receptor_ligand_pair:
            interface_extraction_and_alignment(pdb_file1, pchains1a, pchains1b, pdb_file2, receptor, ligand)

    elif search_interface_chain1 and search_interface_chain2:
        receptor_ligand_pair1 = []
        receptor_ligand_pair2 = []
        valid_super_int_chain_dict1 = filter_iden_int_and_find_large_int(input_stem1, pdb_file1, pchains1a, pchains1b)
        valid_super_int_chain_dict2 = filter_iden_int_and_find_large_int(input_stem2, pdb_file2, pchains2a, pchains2b)

        for ligands, receptors in valid_super_int_chain_dict1.items():
            receptor_ligand_pair_per_key = list(itertools.product(receptors, [ligands]))
            receptor_ligand_pair1.extend(receptor_ligand_pair_per_key)

        for ligands, receptors in valid_super_int_chain_dict2.items():
            receptor_ligand_pair_per_key = list(itertools.product(receptors, [ligands]))
            receptor_ligand_pair2.extend(receptor_ligand_pair_per_key)

        for receptor1, ligand1 in receptor_ligand_pair1:
            for receptor2, ligand2 in receptor_ligand_pair2:
                interface_extraction_and_alignment(pdb_file1, receptor1, ligand1, pdb_file2, receptor2, ligand2)

    else:
        interface_extraction_and_alignment(pdb_file1, pchains1a, pchains1b, pdb_file2, pchains2a, pchains2b)


end = time.time()
print(f"Total Running Time: {end - start}s")
