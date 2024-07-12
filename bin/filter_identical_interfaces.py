#!/usr/bin/env python3

import itertools

from collections import defaultdict
from prody import *
from prody.proteins import compare


min_in_res = 20                             # Minimum number of protein residues of a interface
interface_distance_threshold = 4.5          # Distance threshold for interfaces
inter_interface_distance_threshold = 5.5    # Distance threshold for larger interfaces
seq_id_threshold = 95                       # 95% sequence identity as a threshold for identical protein chains




####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
#######             1. All-against-all pair-wise neighbor search betweeen receptors chain       #######
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
'''The following blocks of code involve in all-against-all interface pair search returns a list of interface chainID pairs'''

def create_chain_atom_obj(pdb_file, pchainsa, pchainsb):
    '''Takes a pdb file and chain IDs.
       Returns a selection object that store all the atom information of the file.
    '''
    structure = parsePDB(pdb_file)
    chain_atoms = structure.select(f"not hetatm and "
                                   f"(altloc A or altloc _) and "
                                   f"chain {' '.join(pchainsa+pchainsb)} and "
                                   f"not resname nonstdaa and "
                                   f"not name OXT")
    
    return chain_atoms


def filter_int_pair(int_pair_dict):
    '''Filters receptor pairs to extract unique receptor pairs that have the largest number of interacting residues.
       Takes a dictionary that stores interacting receptor pairs and its interacting residue count.
       Returns a list of filtered unique receptor pairs.
    '''
    int_pair_filtered = []
    sorted_pairs = sorted(int_pair_dict.items(), key=lambda item: item[1], reverse=True)

    # To keep track of selected elements
    selected_elements = set()

    # To store the selected unique pairs
    selected_pairs = []

    # Iterate through the sorted pairs and select unique pairs
    for pair, count in sorted_pairs:
        if pair[0] not in selected_elements and pair[1] not in selected_elements:
            selected_pairs.append(pair)
            selected_elements.update(pair)

    for pair in selected_pairs:
        int_pair_filtered.append(pair)
    
    return int_pair_filtered
    
    

def find_interacting_chain(chain_selection_obj_dict, combinations):
    '''Finds interacting receptor chain pairs. Take chain_selection_obj_dict and receptor pair combination list as input.
       Returns a list of interacting receptor chain pairs.
    '''
    int_pair = []
    int_pair_dict = {}
    for pair in combinations:
        atoms_pchaina = chain_selection_obj_dict[pair[0]]
        atoms_pchainb = chain_selection_obj_dict[pair[1]]
        
        interacting_interface_pair_residues = measure.contacts.findNeighbors(atoms_pchaina, 4.5, atoms_pchainb)
        num_interface_residues = count_unique_residue_findNeighbor(interacting_interface_pair_residues)
        
        if num_interface_residues >= 15:
            int_pair_dict[pair] = num_interface_residues
        
    
    int_pair = filter_int_pair(int_pair_dict)
    
    return int_pair


def all_against_all_int_search(chain_selection_obj_dict, pchainsa):
    '''Performs all-against-all interacting chain pair search between given receptors.
       Takes a dictionary that stores chainIDs and their corresponding selection objects & receptor chainIDs.
       Returns a list of a list of unique interacting receptor chain pairs
    '''
    combinations_receptor = list(itertools.combinations(pchainsa, 2))
    int_receptor_pair = find_interacting_chain(chain_selection_obj_dict, combinations_receptor)

    return int_receptor_pair


def find_orphan_receptor(int_receptor_pair, pchainsa):
    '''
    Finds receptor chains that are not interacting with any other receptor chains.
    Takes int_receptor_pair from all_against_all_int_search function and receptor chainIDs.
    Returns a list of receptor chains that are not interacting with any other receptor chains.
    '''
    receptor_chain_track = list(pchainsa)
    for pair in int_receptor_pair:
        if pair[0] in receptor_chain_track:
            receptor_chain_track.remove(pair[0])

        if pair[1] in receptor_chain_track:
            receptor_chain_track.remove(pair[1])
            
    orphan_receptor_list = receptor_chain_track

    return orphan_receptor_list


####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
###    2. Create a dictionary of chains that are interacting with each of antigen/ligand chains     ###
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

def count_unique_residue_findNeighbor(neighbor_pair_list):
    '''
    Takes a list object from findNeighbors() function.
    Returns the number of unique interacting residues found in the list object.
    '''
    list1 = []
    list2 = []

    for neighbor_pair in neighbor_pair_list:

        list1.append(neighbor_pair[0].getResindex())
        list2.append(neighbor_pair[1].getResindex())

    interface_res_num = len(set(list1 + list2))

    return interface_res_num



def create_init_ligand_receptor_dict(chain_selection_obj_dict, int_receptor_pair, pchainsb, min_in_res, interface_distance_threshold):
    '''
    Takes chain_selection_obj_dict, int_receptor_pair (list), ligand chainIDs, min_in_res, and interface_distance_threshold.
    Returns a dictionaty that stores...key: ligand chainIDs, value: receptor chainIDs)
    '''
    int_chain_dict = {}
    for ligand in pchainsb:
        for pair in int_receptor_pair:
            
            atoms_pchaina = chain_selection_obj_dict[ligand]
            atoms_pchainb = chain_selection_obj_dict[pair[0]] + chain_selection_obj_dict[pair[1]]

            interacting_interface_pair_residues = measure.contacts.findNeighbors(atoms_pchaina, interface_distance_threshold, atoms_pchainb)  # interface_distance_threshold = 4.5
            num_interface_residues = count_unique_residue_findNeighbor(interacting_interface_pair_residues)

            if atoms_pchaina is not None and atoms_pchainb is not None:
                if num_interface_residues > min_in_res:  # min_in_res = 20
                    if ligand not in int_chain_dict:
                        int_chain_dict[ligand] = []
                    int_chain_dict[ligand].append(pair[0])
                    int_chain_dict[ligand].append(pair[1])

        if ligand in int_chain_dict:
            int_chain_dict[ligand] = list(set(int_chain_dict[ligand]))
    return int_chain_dict
    

def assign_orphan_receptor(orphan_receptor_list, chain_selection_obj_dict, pchainsb, interface_distance_threshold, int_chain_dict):
    '''
    Takes orphan_receptor_list, chain_selection_obj_dict, ligand chainIDs, interface_distance_threshold, and int_chain_dict.
    Updates int_chain_dict by assigning orphan receptors to appropriate ligand chain as interacting chains.
    '''
    min_in_res = 10
    for ligand in pchainsb:
        for receptor in orphan_receptor_list:
            atoms_pchaina = chain_selection_obj_dict[ligand]
            atoms_pchainb = chain_selection_obj_dict[receptor]

            interacting_interface_pair_residues = measure.contacts.findNeighbors(atoms_pchaina, interface_distance_threshold, atoms_pchainb)  # interface_distance_threshold = 4.5
            num_interface_residues = count_unique_residue_findNeighbor(interacting_interface_pair_residues)

            if num_interface_residues > min_in_res:
                if ligand not in int_chain_dict:
                    int_chain_dict[ligand] = []
                int_chain_dict[ligand].append(receptor)


def find_duplicates(input_dict):
    '''
    Takes input_dict.
    Returns a dictionary that stores a receptor chainID as key and its interacting ligand chainIDs as values when
    a single receptor chain is interacting with multiple ligand chains.
    '''
    element_tracker = {}
    duplicates = {}

    # Iterating through each key-value pair
    for key, values in input_dict.items():
        for value in values:
            # If the element is not already in element_tracker, add it with current key
            if value not in element_tracker:
                element_tracker[value] = {key}
            else:
                # If it's already there but associated with a different key, it's a duplicate
                if key not in element_tracker[value]:
                    element_tracker[value].add(key)
                    # Add or update the duplicates dictionary
                    if value in duplicates:
                        duplicates[value].append(key)
                    else:
                        duplicates[value] = [key]

    # Convert the sets in element_tracker to lists for consistency
    for value, keys in duplicates.items():
        duplicates[value] = list(element_tracker[value])

    return duplicates


def Neighbors_pecentage(selectiona, selectionb, neighbor_pair_list):
    '''
    Takes a pair of selection objects (receptor, ligand) and a interacting atoms pairs list from findNeighbors.
    Returns a percentage of actual interacting atoms over the number of all possible atomic interaction between receptor and ligand chains.
    '''
    percentage = 100*(len(neighbor_pair_list) / (len(selectiona) * len(selectionb)))

    return percentage


def calc_percentage_neigboratoms_for_duplicated_receptor(int_chain_dict, ligand_receptor_pair_duplicate, chain_selection_obj_dict, int_receptor_pair):
    '''
    Takes int_chain_dict, ligand_receptor_pair_duplicate(list), chain_selection_obj_dict, and int_receptor_pair(list).
    Returns the average percentage value of the interaction densities calculated by three different methods.
    '''
    ligand = ligand_receptor_pair_duplicate[1]
    receptor = ligand_receptor_pair_duplicate[0]
    chain_select_list = int_chain_dict[ligand].copy()
    
    # Method 1
    # Calculate the interaction density for all the receptors bound to the ligand versus the ligand itself
    atoms_receptor = chain_selection_obj_dict[receptor]
    atoms_ligand = chain_selection_obj_dict[ligand]
    atoms_receptors = chain_selection_obj_dict[chain_select_list[0]]

    for chain in chain_select_list[1:]:
        atoms_receptors += (chain_selection_obj_dict[chain])

    pair_residues1 = measure.contacts.findNeighbors(atoms_ligand, 4.5, atoms_receptors)
    percentage1 = Neighbors_pecentage(atoms_ligand, atoms_receptors, pair_residues1)
    
    
    # Method 2
    # Calculate the interaction density for a duplicate receptor bound to the ligand versus the rest of the receptors plus the ligand
    chain_select_list.remove(receptor)
    chain_select_list.append(ligand)
    atoms_receptor = chain_selection_obj_dict[receptor]
    atoms_ligand = chain_selection_obj_dict[ligand]
    atoms_other_receptor_plus__ligand = chain_selection_obj_dict[chain_select_list[0]]

    for chain in chain_select_list[1:]:
        atoms_other_receptor_plus__ligand += (chain_selection_obj_dict[chain])

    pair_residues2 = measure.contacts.findNeighbors(atoms_receptor, 7, atoms_other_receptor_plus__ligand)
    percentage2 = Neighbors_pecentage(atoms_receptor, atoms_other_receptor_plus__ligand, pair_residues2)
    

    # Method 3
    # Calculate the interaction density for the ligand-receptor and its partner versus the ligand
    atoms_receptor = chain_selection_obj_dict[receptor]
    atoms_receptor = chain_selection_obj_dict[receptor]

    for pair in int_receptor_pair:
        if receptor in pair:
            if receptor in pair[0]:
                receptor_partner = pair[1]
                atoms_receptor += (chain_selection_obj_dict[receptor_partner])

            elif receptor in pair[1]:
                receptor_partner = pair[0]
                atoms_receptor += (chain_selection_obj_dict[receptor_partner])

    atoms_ligand = chain_selection_obj_dict[ligand]
    pair_residues3 = measure.contacts.findNeighbors(atoms_receptor, 10, atoms_ligand)
    percentage3 = Neighbors_pecentage(atoms_receptor, atoms_ligand, pair_residues3)

    percentage = (percentage1 + percentage2 + percentage3) / 3
    
    return percentage
    


def find_rec_lig_pair_to_delete(duplicates, int_chain_dict, chain_selection_obj_dict, int_receptor_pair):
    '''
    Identifies receptor-ligand pairs to delete based on interaction density percentages.

    Parameters:
    duplicates (dict): A dictionary where keys are shared receptor chains and values are lists of ligand chains sharing those receptors.
    int_chain_dict (dict): A dictionary mapping ligands to their interaction chains.
    chain_selection_obj_dict (dict): A dictionary mapping chain identifiers to their selected atom objects.
    int_receptor_pair (list): A list of interaction receptor pairs.

    Returns a list of receptor-ligand pairs (tuples) that should be deleted based on lower interaction density percentages.
    '''
    delete_pair_list = []

    for key, value in duplicates.items():  # key: shared receptor chain, value: a ligand chain that shares a specific receptor chain
        max_percentage = 0
        max_percentage_object = None
        
        # Generate interface pairs that can be used as keys for interface_pair_Neighbor_result dictionary
        ligand_receptor_combinations = list(itertools.product(key, value))
        
        for ligand_receptor_pair_duplicate in ligand_receptor_combinations:
            percentage = calc_percentage_neigboratoms_for_duplicated_receptor(int_chain_dict, ligand_receptor_pair_duplicate, chain_selection_obj_dict, int_receptor_pair)

            if percentage > max_percentage:
                max_percentage = percentage
                max_percentage_object = ligand_receptor_pair_duplicate
            
        # Collect all other neighbors_objects into delete_pair_list
        for pair in ligand_receptor_combinations:
            if pair != max_percentage_object:
                delete_pair_list.append(pair)

    return delete_pair_list


def filter_duplicate_chain(duplicates, int_chain_dict, chain_selection_obj_dict, int_receptor_pair):
    '''
    Filters out duplicate receptor-ligand pairs from the interaction chain dictionary.

    Parameters:
    duplicates (dict): A dictionary where keys are shared receptor chains and values are lists of ligand chains sharing those receptors.
    int_chain_dict (dict): A dictionary mapping ligands to their interaction chains.
    chain_selection_obj_dict (dict): A dictionary mapping chain identifiers to their selected atom objects.
    int_receptor_pair (list): A list of interaction receptor pairs.
    '''
    delete_pair_list = find_rec_lig_pair_to_delete(duplicates, int_chain_dict, chain_selection_obj_dict, int_receptor_pair)
    for delete_pair in delete_pair_list:
        int_chain_dict[delete_pair[1]].remove(delete_pair[0])
    
    # Delete ligand key if it has an empty value list
    ligand_key_to_remove = []
    for ligand, receptors in int_chain_dict.items():
        if len(receptors) == 0:
            ligand_key_to_remove.append(ligand)
    
    for ligand in ligand_key_to_remove:
        int_chain_dict.pop(ligand)
            


####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
###    3. perform the sequence alignments of chains to identify/filter out identical ligand and receptor chains     ###
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

def create_ligand_chain_res_dict(int_chain_dict, chain_selection_obj_dict):
    ''' Create a dictionary that stores ligand chainIDs as keys and selection objects (used for alignment) as values.'''
    ligand_chain_res_dict = {}  # keys: ligand chainID, values: selection objects (used for alignment)
    for key in int_chain_dict.keys():
        atoms_pchain = chain_selection_obj_dict[key]
        ligand_chain_res_dict[key] = atoms_pchain

    return ligand_chain_res_dict


def group_identical_chains(pairs):
    '''Use graph theory approach with Depth-First Search (DFS) to correctly assign chains to each of identical chain groups
       Returns a list of sublists that represents identical chain groups
    '''
    # Create a graph using adjacency list representation
    graph = defaultdict(list)
    
    # Add edges to the graph
    for u, v in pairs:
        graph[u].append(v)
        graph[v].append(u)
    
    # Function to perform DFS and find a connected component
    def dfs(node, visited, component):
        visited.add(node)
        component.append(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(neighbor, visited, component)
    
    visited = set()
    components = []

    # Find all identical chain groups
    for node in graph:
        if node not in visited:
            component = []
            dfs(node, visited, component)
            components.append(sorted(component))
    
    return components


def calc_seq_identity_ligand(chain_selection_obj_dict, int_chain_dict):
    '''
    Calculates sequence identity for ligand chains in the given interaction chain dictionary.

    Parameters:
    chain_selection_obj_dict (dict): A dictionary mapping chain identifiers to their selected atom objects.
    int_chain_dict (dict): A dictionary mapping ligands to their interaction chains.

    Returns a comparison object containing sequence identity results of the ligand chains.
    '''
    int_chain_dict_keys = list(int_chain_dict.keys())
    ligand_list = int_chain_dict_keys
    chains = chain_selection_obj_dict[ligand_list[0]]

    for chain in ligand_list[1:]:
        chains += (chain_selection_obj_dict[chain])

    chains = chains.toAtomGroup()
    results = compare.matchChains(chains, chains, pwalign=True, overlap=50)
    
    return results


def seq_id_check_ligand(chain_selection_obj_dict, int_chain_dict, seq_id_threshold):
    '''
    Identifies ligand chains with sequence identities above 95%.
    Returns a identical_ligand_list (list) that stores a list of groups of identical ligand chains based on sequence identity.
    '''
    identical_ligand_list = []
    results = calc_seq_identity_ligand(chain_selection_obj_dict, int_chain_dict)
    processed_pairs = []

    if results is not None:
        for result in results:
            ligand1 = result[0].getChids()[0]
            ligand2 = result[1].getChids()[0]
            seq_id = result[2]

            if ligand1 != ligand2:
                pair = tuple(sorted([ligand1, ligand2]))

                # Check if this pair has already been processed
                if pair in processed_pairs:
                    continue  # Skip this pair if it has already been processed
                
                if seq_id > seq_id_threshold:
                    processed_pairs.append(pair)

        identical_chain_groups_set = group_identical_chains(processed_pairs)
        
        for identical_chain_group in identical_chain_groups_set:
            identical_ligand_list.append(identical_chain_group)

    return identical_ligand_list


def calc_seq_identity_receptor(chain_selection_obj_dict, int_chain_dict):
    '''
    Calculates sequence identity for receptor chains in the given interaction chain dictionary.

    Parameters:
    chain_selection_obj_dict (dict): A dictionary mapping chain identifiers to their selected atom objects.
    int_chain_dict (dict): A dictionary mapping ligands to their interaction receptor chains.

    Returns a comparison object containing sequence identity results of the receptor chains.
    '''
    int_chain_dict_values = list(int_chain_dict.values())
    receptor_list = int_chain_dict_values[0].copy()

    for recepor_list in int_chain_dict_values[1:]:
        receptor_list.extend(recepor_list)

    chains = chain_selection_obj_dict[receptor_list[0]]

    for chain in receptor_list[1:]:
        chains += (chain_selection_obj_dict[chain])

    chains = chains.toAtomGroup()
    results = matchChains(chains, chains, pwalign=True, overlap=50)

    return results


def seq_id_check_receptor(chain_selection_obj_dict, int_chain_dict, seq_id_threshold):
    '''
    Identifies receptor chains with sequence identities above a 95%.
    Returns a identical_receptor_dict (dictionary) that stores a list of groups of identical ligand chains based on sequence identity.
    '''
    identical_receptor_dict = {}
    results = calc_seq_identity_receptor(chain_selection_obj_dict, int_chain_dict)
    processed_pairs = []

    if results is not None:
        for result in results:
            receptor1 = result[0].getChids()[0]
            receptor2 = result[1].getChids()[0]
            seq_id = result[2]

            if receptor1 != receptor2:
                pair = tuple(sorted([receptor1, receptor2]))

                # Check if this pair has already been processed
                if pair in processed_pairs:
                    continue  # Skip this pair if it has already been processed
                
                if seq_id > seq_id_threshold:
                    processed_pairs.append(pair)

        identical_chain_groups_set = group_identical_chains(processed_pairs)
        
        for identical_chain_group in identical_chain_groups_set:
            identical_receptor_dict[identical_chain_group[0]] = []
            identical_receptor_dict[identical_chain_group[0]].extend(list(identical_chain_group[1:]))

    return identical_receptor_dict


def map_identical_interactions(int_chain_dict, identical_receptor_comb_list):
    '''
    Maps identical interactions between ligands and receptors based on sequence identity.

    Parameters:
    int_chain_dict (dict): A dictionary mapping ligands to their interaction receptor chains.
    identical_receptor_comb_list (list): A list of lists where each inner list contains a group of receptors that are sequencially identical.

    Returns a dictionary where keys are identifiers of identical interaction combinations and values are lists of these combinations.
    '''
    int_chain_dict_values = list(int_chain_dict.values())
    receptor_track_list = int_chain_dict_values[0].copy()

    for recepor_list in int_chain_dict_values[1:]:
        receptor_track_list.extend(recepor_list)

    def find_equivalent_combinations(receptors, identical_receptor_comb_list, receptor_track_list):
        # Find equivalent receptors for each receptor in the list
        equivalent_groups = []
        
        for receptor in receptors:
            equivalents = next((group for group in identical_receptor_comb_list if receptor in group and receptor in receptor_track_list), [])
            equivalent_groups.append(list(equivalents))
    
        # Generate all possible combinations of equivalent receptors
        equivalent_combinations = set(itertools.product(*equivalent_groups))
        
        return equivalent_combinations
    
    def find_equivalent_receptors(receptor, identical_receptor_comb_list):
        equivalents = set()

        # Find the group that contains the receptor
        receptors_equiv = next((group for group in identical_receptor_comb_list if receptor in group and receptor in receptor_track_list), [])
        
        for eq in receptors_equiv:
            equivalents.add(eq)
       
        return equivalents

    all_ligands = int_chain_dict.keys()
    identical_combinations = {}
    unique_combinations = set()
    comb_counter = 1

    for receptor in int_chain_dict.values():

        if len(receptor) == 1:
            equivalent_receptors = find_equivalent_receptors(receptor, identical_receptor_comb_list)

            identical_comb = set()
            for eq in equivalent_receptors:
                for lig in all_ligands:
                    if eq in int_chain_dict.get(lig, set()):
                        interaction = (lig, eq)
                        identical_comb.add(interaction)
            
            identical_comb_list = sorted(list(identical_comb))
            if identical_comb_list:
                identical_comb_tuple = tuple(identical_comb_list)
                if identical_comb_tuple not in unique_combinations:
                    unique_combinations.add(identical_comb_tuple)
                    identical_combinations[f'identical_comb{comb_counter}'] = identical_comb_list
                    comb_counter += 1
        
        equivalent_combinations = find_equivalent_combinations(receptor, identical_receptor_comb_list, receptor_track_list)
        identical_comb = set()
        
        for eq_comb in equivalent_combinations:
            for lig in all_ligands:
                if all(eq in int_chain_dict.get(lig, set()) for eq in eq_comb):
                    interaction = (lig, *eq_comb)
                    identical_comb.add(interaction)
        
        # Normalize the identical combinations
        normalized_comb = tuple(sorted(tuple(sorted(comb)) for comb in identical_comb))

        if len(normalized_comb) != 0 and normalized_comb not in unique_combinations:
            unique_combinations.add(normalized_comb)
            identical_combinations[f'identical_comb{comb_counter}'] = sorted(list(identical_comb))
            comb_counter += 1

    return identical_combinations


def reorder_identical_interaction_set_dict(identical_interaction_set_dict):
    '''
    Reorders a dictionary of identical interaction sets from map_identical_interaction based on the minimum length of tuples in the values.

    Parameters:
    identical_interaction_set_dict (dict): A dictionary where keys are identifiers and values are lists of tuples representing identical interaction sets.

    Returns a reordered dictionary where items are sorted in a ascending order in terms of the tuples length in the values.
    '''
    # Sort the dictionary items based on the minimum length of tuples in the values
    sorted_items = sorted(identical_interaction_set_dict.items(), key=lambda item: min(len(tup) for tup in item[1]))
    
    # Convert the sorted items back to a dictionary
    reordered_map = {key: value for key, value in sorted_items}
    
    return reordered_map


def refine_identical_interaction_set_dict(identical_interaction_set_dict):
    '''
    Refines a dictionary of identical interaction sets by removing duplicate receptors and filtering out sublists with length 1.
    '''
    seen_list = []
    for key, recepetors in identical_interaction_set_dict.items():
        # Convert each tuple in the list to a list
        identical_interaction_set_dict[key] = [list(tup) for tup in recepetors]

    for value in identical_interaction_set_dict.values():
        for sublist in value:
            for receptor in sublist[1:]:
                if receptor not in seen_list:
                    seen_list.append(receptor)
                else:
                    sublist.remove(receptor)

    for key, receptors in identical_interaction_set_dict.items():
        updated_value = []
    
        # Iterate over each sublist in the current value list
        for sublist in receptors:
            # Check if the sublist has length greater than 1 (to keep it)
            if len(sublist) > 1:
                updated_value.append(sublist)  # Add sublist to updated_value if len(sublist) > 1
                
        # Update the dictionary with the filtered value
        identical_interaction_set_dict[key] = updated_value


def find_best_interfaces(identical_interaction_set, chain_selection_obj_dict):
    ''' For each set of identical interfaces, find a single interface that has the most number of interfacing chains'''
    best_unique_int_dict = {}
    for identical_int_combs in identical_interaction_set.values():
        max_num = 0
        for lig_receptor_comb in identical_int_combs:

            ligand = chain_selection_obj_dict[lig_receptor_comb[0]]
            receptors = chain_selection_obj_dict[lig_receptor_comb[1]]
            for receptor in lig_receptor_comb[2:]:
                receptors += chain_selection_obj_dict[receptor]
            
            interacting_interface_pair_residues = measure.contacts.findNeighbors(ligand, interface_distance_threshold, receptors)
            num_interface_residues = count_unique_residue_findNeighbor(interacting_interface_pair_residues)
            
            if num_interface_residues > max_num:
                max_num = num_interface_residues
                max_num_object = lig_receptor_comb

        if max_num_object[0] not in best_unique_int_dict:
            best_unique_int_dict[max_num_object[0]] = []
        
        for receptor in max_num_object[1:]:
            best_unique_int_dict[max_num_object[0]].append(receptor)
    return best_unique_int_dict


def track_non_unique_receptors(identical_receptor_dict, best_unique_int_dict):
    ''' Keep track of identical receptor chains by creating a "non_unique_receptors" list.'''
    identical_receptor_flatten_list = [receptor_prime for receptor_prime in identical_receptor_dict] + [receptor_secondary for receptor_secondary_list in identical_receptor_dict.values() for receptor_secondary in receptor_secondary_list]
    unique_receptor_flatten_list = [receptor for receptor_list in best_unique_int_dict.values() for receptor in receptor_list]

    for unique_receptor in unique_receptor_flatten_list:
        if unique_receptor in identical_receptor_flatten_list:
            identical_receptor_flatten_list.remove(unique_receptor)
    
    non_unique_receptors = identical_receptor_flatten_list

    return non_unique_receptors


def track_ligands_in_identical_interaction_set_dict(identical_interaction_set_dict):
    '''
    Tracks and returns a list of unique ligands in the identical interaction set dictionary.

    Parameters:
    identical_interaction_set_dict (dict): A dictionary where keys are identifiers and values are lists of lists representing identical interaction sets.

    Returns a list of unique ligand chain IDs found in the input identical interaction set dictionary.
    '''

    ligand_list = []
    for key, sublist in identical_interaction_set_dict.items():
        # Extract the first value from each sublist and add it to first_values_list
        for sublist_item in sublist:
            if sublist_item:  # Check if the sublist item is not empty
                if sublist_item[0] not in ligand_list:
                    ligand_list.append(sublist_item[0])
    
    return ligand_list


def filter_identical_interface(int_chain_dict, identical_ligand_list, unique_int_chain_dict, non_unique_receptors, identical_interaction_set_dict):
    '''
    Filters the identical interfaces by removing non-unique receptors from the unique_int_chain_dict.
    '''
    ligand_list_in_int_map = track_ligands_in_identical_interaction_set_dict(identical_interaction_set_dict)
    for ligand, receptors in int_chain_dict.items():
        for group in identical_ligand_list:
            if ligand not in group and ligand not in ligand_list_in_int_map:

                unique_int_chain_dict[ligand] = receptors

            else:
                unique_receptors = [r for r in receptors if r not in non_unique_receptors]

                if unique_receptors:
                    unique_int_chain_dict[ligand] = unique_receptors
