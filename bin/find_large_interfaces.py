#!/usr/bin/env python3

import filter_identical_interfaces

from prody import *
from filter_identical_interfaces import *


min_in_res = 20                             # Minimum number of protein residues of a interface
interface_distance_threshold = 4.5          # Distance threshold for interfaces
inter_interface_distance_threshold = 5.5    # Distance threshold for larger interfaces
seq_id_threshold = 95                       # 95% sequence identity as a threshold for identical protein chains



####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####
### 1. Calculate the pairwise distance of residues in two interfaces to check if there are any super interface   ###
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####

# A super interface is formed when there are more than interacting chains, and more than 2 res from each chain) (threshold: 5.5)

def find_close_pairs(value, chain_selection_obj_dict, distance_threshold=5.5, interaction_threshold=20):
    '''
    Finds and returns pairs of receptors that have a number of interactions above a specified threshold.
    Returns a list of receptor pairs that have interactions above the specified threshold.
    '''
    rec_rec_dist_check_pairs = list(itertools.combinations(value, 2))
    rec_rec_close_pairs = []
    
    for pair in rec_rec_dist_check_pairs:
        interface_atoms_chain_ligand1 = chain_selection_obj_dict[pair[0]]
        interface_atoms_chain_ligand2 = chain_selection_obj_dict[pair[1]]

        interacting_ligand_atom_pairs = measure.contacts.findNeighbors(interface_atoms_chain_ligand1, distance_threshold, interface_atoms_chain_ligand2)
        num = count_unique_residue_findNeighbor(interacting_ligand_atom_pairs)
        
        if num >= interaction_threshold:
            rec_rec_close_pairs.append(pair)
    
    return rec_rec_close_pairs


def count_chain_residue_composition_findNeighbor(neighbor_pair_list):
    '''
    Calculates the residue composition for each chain involved in p-p interfaces.

    Parameters:
    neighbor_pair_list (list): A list of neighbor pairs generated by findNeighbors() function

    Returns a dictionary where keys are chain identifiers and values are lists of unique residue numbers involved in interactions for each chain.
    '''
    composition = {}
    for neighbor_pair in neighbor_pair_list:
        chain_flag1 = neighbor_pair[0].getChid()
        chain_flag2 = neighbor_pair[1].getChid()
        if chain_flag1 not in composition:
            composition[chain_flag1] = []

        if chain_flag2 not in composition:
            composition[chain_flag2] = []
       
        composition[chain_flag1].append(neighbor_pair[0].getResnum())
        composition[chain_flag2].append(neighbor_pair[1].getResnum())

    for key, values in composition.items():
        composition[key] = list(set(values))

    return composition


def find_large_interface(unique_int_chain_dict, chain_atoms, chain_selection_obj_dict):
    """
    Identifies and constructs larger interface receptor chains within the same ligand based on unique ligand-receptor interactions.
    
    Parameters:
    unique_int_chain_dict (dict): Dictionary with ligands as keys and their corresponding unique receptor chains as values.
    chain_atoms: Atom selection object containing chain atoms.
    chain_selection_obj_dict (dict): Dictionary with chain identifiers as keys and their corresponding atom selection objects as values.
    
    Returns a dictionary with ligands as keys and their corresponding super interface chains as values.
    """
    super_int_chain_dict = {}
    # Run a for loop for each unique ligand(key)-receptor(value)
    for ligand, receptors in unique_int_chain_dict.items():
        temp_receptors_check_list = []

        # a. receptor-receptor distance check --> pick up a combination of ligands for large interface calculation
        rec_rec_close_pairs = find_close_pairs(receptors, chain_selection_obj_dict)

        count = 0
        while receptors != temp_receptors_check_list:
            temp_chain_check_list = list(receptors.copy())

            # Create a list of interface chainID pairs
            if count == 0:
                chainIDs_interface_alignment_pairs = list(itertools.product(ligand, rec_rec_close_pairs))

            if count > 0:  # Even larger interface
                chainIDs_interface_alignment_pairs = list(itertools.combinations(receptors, 2))

            # Create a list that track the previous receptors list. Used as a flag to stop the while loop when there's no additional super interface chain combos
            temp_receptors_check_list = list(receptors)

            # b. Check if there are pairs of interfaces that are interacting with each other (finding a larger interface).
            for pair in chainIDs_interface_alignment_pairs:
                interacting_interface_atom_pairs = None
                composition_check_list = None
                if count == 0:
                    interface_atoms_chain_receptor1_ligand = chain_atoms.select(f"chain {pair[0]} within 5 of chain {pair[1][0]}")
                    interface_atoms_chain_receptor1_receptor = chain_atoms.select(f"chain {pair[1][0]} within 5 of chain {pair[0]}")
                    interface_atoms_chain_receptor2_ligand = chain_atoms.select(f"chain {pair[0]} within 5 of chain {pair[1][1]}")
                    interface_atoms_chain_receptor2_receptor = chain_atoms.select(f"chain {pair[1][1]} within 5 of chain {pair[0]}")

                    if interface_atoms_chain_receptor1_ligand is not None and interface_atoms_chain_receptor1_receptor is not None:
                        interface_atoms_chain_receptor1 = interface_atoms_chain_receptor1_ligand + interface_atoms_chain_receptor1_receptor

                        if interface_atoms_chain_receptor2_ligand is not None and interface_atoms_chain_receptor2_receptor is not None:
                            interface_atoms_chain_receptor2 = interface_atoms_chain_receptor2_ligand + interface_atoms_chain_receptor2_receptor

                            interacting_interface_atom_pairs = measure.contacts.findNeighbors(interface_atoms_chain_receptor1, inter_interface_distance_threshold, interface_atoms_chain_receptor2)  # Maybe threshold can be a bit more loose..?
                            composition_check_list = [pair[0], pair[1][0], pair[1][1]]


                if count > 0:
                    interface_atoms_chain_receptor1_ligand = chain_atoms.select(f"chain {ligand} within 5 of chain {' '.join(pair[0])}")
                    interface_atoms_chain_receptor1_receptor = chain_atoms.select(f"chain {' '.join(pair[0])} within 5 of chain {ligand}")
                    interface_atoms_chain_receptor2_ligand = chain_atoms.select(f"chain {ligand} within 5 of chain {' '.join(pair[1])}")
                    interface_atoms_chain_receptor2_receptor = chain_atoms.select(f"chain {' '.join(pair[1])} within 5 of chain {ligand}")

                    if interface_atoms_chain_receptor1_ligand is not None and interface_atoms_chain_receptor1_receptor is not None:
                        interface_atoms_chain_receptor1 = interface_atoms_chain_receptor1_ligand + interface_atoms_chain_receptor1_receptor

                        if interface_atoms_chain_receptor2_ligand is not None and interface_atoms_chain_receptor2_receptor is not None:
                            interface_atoms_chain_receptor2 = interface_atoms_chain_receptor2_ligand + interface_atoms_chain_receptor2_receptor

                            interacting_interface_atom_pairs = measure.contacts.findNeighbors(interface_atoms_chain_receptor1, inter_interface_distance_threshold, interface_atoms_chain_receptor2)
                            
                            composition_check_list = [char for string in pair for char in string]
                            composition_check_list.append(ligand)

                if interacting_interface_atom_pairs is not None:
                    composition = count_chain_residue_composition_findNeighbor(interacting_interface_atom_pairs)
  
                all_lengths_greater_than_3 = True
                if composition_check_list is not None:
                    for i in composition_check_list:
                        # Check if the receptors is None
                        if i not in composition:
                            all_lengths_greater_than_3 = False
                            break
                        # Check if the length of the receptors is less than 3
                        elif len(composition[i]) <= 3:
                            all_lengths_greater_than_3 = False
                            break
                

                if all_lengths_greater_than_3 and composition_check_list is not None:
                    if count == 0:
                        if pair[1][0] in temp_chain_check_list:  # IF statement to handle the case when some chain in superchains are overlapped (ex. "L" of HL & LM)
                            temp_chain_check_list.remove(pair[1][0])

                        if pair[1][1] in temp_chain_check_list:
                            temp_chain_check_list.remove(pair[1][1])

                        # Update unique_int_chain_dict by delete chainIDs of larger interface from chainID receptors list
                        if pair[1][0] in receptors:
                            receptors.remove(pair[1][0])

                        if pair[1][1] in receptors:
                            receptors.remove(pair[1][1])
                        
                        receptors.append(''.join([pair[1][0], pair[1][1]]))
                    
                    if count > 0:
                        if pair[0] in temp_chain_check_list:  # IF statement to handle the case when some chain in superchains are overlapped (ex. "L" of HL & LM)
                            temp_chain_check_list.remove(pair[0])

                        if pair[1] in temp_chain_check_list:
                            temp_chain_check_list.remove(pair[1])

                        # Update unique_int_chain_dict by delete chainIDs of larger interface from chainID receptors list
                        if pair[0] in receptors:
                            receptors.remove(pair[0])

                        if pair[1] in receptors:
                            receptors.remove(pair[1])
                        
                        receptors.append(''.join([pair[0], pair[1]]))

                else:
                    if count == 0:
                        if pair[1][0] in receptors:
                            receptors.remove(pair[1][0])

                        if pair[1][1] in receptors:
                            receptors.remove(pair[1][1])

                    if count > 0:
                        if pair[0] in receptors:
                            receptors.remove(pair[0])

                        if pair[1] in receptors:
                            receptors.remove(pair[1])

            # Move remaining elements in temp_chain_check_list's value elements as non-interacting interface receptor chain into the final interface chain dictionary.
            for non_interacting_int_chain in temp_chain_check_list:
                if ligand not in super_int_chain_dict:
                    super_int_chain_dict[ligand] = []
                if non_interacting_int_chain not in super_int_chain_dict[ligand]:
                    super_int_chain_dict[ligand].append(non_interacting_int_chain)
            
            if receptors == []:
                break
            
            count += 1
    return super_int_chain_dict


def find_large_interface_across_ligand(super_int_chain_dict, chain_atoms):
    """
    Identifies and constructs larger interface chains across ligands based on interactions.
    
    Parameters:
    super_int_chain_dict (dict): Dictionary with ligands as keys and their corresponding super interface chains as values.
    chain_atoms: Atom selection object containing chain atoms.
    
    Returns an updated dictionary with combined ligands and their corresponding super interface chains.
    """
    # Generate all possible combinations
    combinations = list(itertools.product(*[[(ligand, receptor) for receptor in receptors]
                                            for ligand, receptors in super_int_chain_dict.items()]))

    for combo in combinations:
        ligand1 = combo[0][0]
        receptors1 = combo[0][1]
        ligand2 = combo[1][0]
        receptors2 = combo[1][1]

        interface_atoms_chain_receptor1_ligand = chain_atoms.select(f"chain {' '.join(ligand1)} within 4.5 of chain {' '.join(receptors1)}")
        interface_atoms_chain_receptor1_receptor = chain_atoms.select(f"chain {' '.join(receptors1)} within 4.5 of chain {' '.join(ligand1)}")
        interface_atoms_chain_receptor2_ligand = chain_atoms.select(f"chain {' '.join(ligand2)} within 4.5 of chain {' '.join(receptors2)}")
        interface_atoms_chain_receptor2_receptor = chain_atoms.select(f"chain {' '.join(receptors2)} within 4.5 of chain {' '.join(ligand2)}")

        if interface_atoms_chain_receptor1_ligand is not None and interface_atoms_chain_receptor1_receptor is not None:
            interface_atoms_chain_receptor1 = interface_atoms_chain_receptor1_ligand + interface_atoms_chain_receptor1_receptor

            if interface_atoms_chain_receptor2_ligand is not None and interface_atoms_chain_receptor2_receptor is not None:
                interface_atoms_chain_receptor2 = interface_atoms_chain_receptor2_ligand + interface_atoms_chain_receptor2_receptor

                interacting_interface_atom_pairs = measure.contacts.findNeighbors(interface_atoms_chain_receptor1, inter_interface_distance_threshold, interface_atoms_chain_receptor2)
                
                composition_check_list = [char for tup in combo for element in tup for char in element]
                composition_check_list = set(composition_check_list)

        if interacting_interface_atom_pairs is not None:
            composition = count_chain_residue_composition_findNeighbor(interacting_interface_atom_pairs)

        all_lengths_greater_than_3 = True
        if composition_check_list is not None:
            for i in composition_check_list:
                # Check if the receptors is None
                if i not in composition:
                    all_lengths_greater_than_3 = False
                    break
                # Check if the length of the receptors is less than 3
                elif len(composition[i]) <= 3:
                    all_lengths_greater_than_3 = False
                    break
        
    
        if all_lengths_greater_than_3 and composition_check_list is not None:
            super_int_chain_dict[ligand1].remove(receptors1)
            super_int_chain_dict[ligand2].remove(receptors2)

            super_ligand = ligand1+ligand2
            super_receptor = ', '.join(set(receptors1+receptors2))
            if super_ligand not in super_int_chain_dict:
                super_int_chain_dict[super_ligand] = []
            super_int_chain_dict[super_ligand].append(super_receptor)
        
            super_int_chain_dict = {key: value for key, value in super_int_chain_dict.items() if value}
    
    return super_int_chain_dict



####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #####
#### 2. Final check: confirm the number of interfacing residues is greater than or equal to 20.  ####
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######  ####
                
def find_valid_super_int(chain_selection_obj_dict, super_int_chain_dict):
    """
    Identifies valid super interfaces between ligands and receptors based on interaction thresholds.
    
    Parameters:
    chain_selection_obj_dict (dict): Dictionary with chain selections.
    super_int_chain_dict (dict): Dictionary with ligands as keys and their corresponding super interface receptor chains as values.
    
    Returns a dictionary with valid ligands and their corresponding receptors that meet the interaction criteria ( >= 20 interface resiudes).
    """
    valid_super_int_chain_dict = {}
    for ligand, receptor in super_int_chain_dict.items():
        chainIDs_interface_alignment_pairs = list(itertools.product([ligand], receptor))

        for pair in chainIDs_interface_alignment_pairs:
            ligand = pair[0]
            receptor = pair[1]

            if len(ligand) > 1:
                interface_atoms_chain_ligand = chain_selection_obj_dict[ligand[0]]  # For cases where ligand has more than one chain
                for chain in ligand[1:]:
                    interface_atoms_chain_ligand += chain_selection_obj_dict[chain]

            elif len(ligand) == 1:
                interface_atoms_chain_ligand = chain_selection_obj_dict[pair[0]]


       
            if len(receptor) > 1:
                interface_atoms_chain_receptor = chain_selection_obj_dict[receptor[0]]  # For cases where receptor has more than one chain
                for chain in receptor[1:]:
                    interface_atoms_chain_receptor += chain_selection_obj_dict[chain]

            if len(receptor) == 1:
                interface_atoms_chain_receptor = chain_selection_obj_dict[receptor[0]]
            
            interacting_interface_atom_pairs = measure.contacts.findNeighbors(interface_atoms_chain_ligand, interface_distance_threshold, interface_atoms_chain_receptor)
            int_res_num = count_unique_residue_findNeighbor(interacting_interface_atom_pairs)
            
            if int_res_num >= 20:
                if ligand not in valid_super_int_chain_dict:
                    valid_super_int_chain_dict[ligand] = []
                valid_super_int_chain_dict[ligand].append(receptor)

    return valid_super_int_chain_dict
