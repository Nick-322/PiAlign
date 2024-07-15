## PiAlign: a tool for the flexible structural comparison of protein–protein interfaces free of chain number restriction

Proteins play fundamental roles in the living cells of all organisms, responsible for a diverse array of biological functions. ​To elucidate protein functionality, it is crucial to understand their structural interactions, particularly the regions where physical contact occurs between proteins, known as protein-protein interfaces. Comparative analysis of these interfaces facilitates the understanding of biological and evolutionary relationships among proteins. In 2010, iAlign was developed to enable efficient, length-independent structural comparison of interfaces formed by two protein chains. Building upon this foundation, PiAlign is a computational tool that extends iAlign's capabilities by allowing structural comparison of protein-protein interfaces formed by any number of chains, not limited to just two. This expansion of the interface analysis scope enables users to capture more comprehensive structural relationships between large protein complexes, thereby enhancing our understanding of protein interactions and their biological implications.


## Installation

#### First, clone this repository by

	git clone https://github.com/Nick-322/PiAlign

#### To use PiAlign, you will have two choices to create a virtual environment with required dependencies.

"alias pialign='<PATH_TO>/PiAlign/bin/PiAlign.py'"
"alias pialign='<PATH_TO>/PiAlign/bin/PiAlign.py'"
* #### Choice 1. Conda

	*If conda is not installed in your computer, refer to the guide from an official site:
	https://conda.io/projects/conda/en/latest/user-guide/install/index.html


	Once codnda is installed, create an environment by running the following code:
	* For MacOS (ARM64) users:
		> CONDA_SUBDIR=osx-64 conda env create -f <PATH_TO>/PiAlign/environment.yml


	* For Linux (x86-64) users:
		> conda env create -f environment.yml


* #### Choice 2. Python venv
	* Create an environment by running the following code:

			python -m venv <environment_name>


	* And run the code below:

			source <PATH_TO_environment>/bin/activate


	* Lastly, run:

			pip install -r requirements.txt


#### Make sure to make the python/binary file executable by running: 

	chmod +x <PATH_TO>/PiAlign/bin/PiAlign.py
	
	
 * for Linux users:
 	
		chmod +x <PATH_TO>/PiAlign/bin/IS-align_linux
	
* for Mac users:
 			
		chmod +x <PATH_TO>/PiAlign/bin/IS-align_mac


#### In the created environment, you may set the following alias:

	alias pialign='<PATH_TO>/PiAlign/bin/PiAlign.py'

#### and invoke PiAlign by

  	pialign


#### *If you can't run PiAlign successfully, you may need to install dssp in your environment. ###
* For Linux users: run 

		sudo apt-get install dssp

* For Mac users:
1. Install MacPorts
	* a. Install Apple's Command Line Developer Tools: 
	    		
			xcode-select --install
	  
	* b. Install MacPorts for your version of the Mac operating system:
			macOS Sonoma v14
			macOS Ventura v13
			macOS Monterey v12
			macOS Big Sur v11

2. Run 			

		sudo port install dssp

## Examples

Two simple examples illustrating usage of PiAlign

### i. Align two protein-protein interfaces (one chain for each side of the interface)

	../bin/PiAlign.py -p1 1lyl.pdb -c1a A -c1b C -p2 12as.pdb -c2a A -c2b B



### ii. Align two protein-protein interfaces with a multiple interface chain case where we don't know which chains interact with which chains for the first pdb file.

	../bin/PiAlign.py -p1 7e5s.pdb -c1a LEIJPRUVDHKTNOQS -c1b BCA -p2 7uap.pdb -c2a HL -c2b A -searchIntCh1



## Output files

This tool generates a pair of two main files: int.pdb and con.lst files.

* #### int.pdb 

int.pdb file contains atomic information of alpha carbons of residues involved in an interface. Each entry has the same types of information as PDB file except for Occupancy, which is replaced with secondary structure assignment (1: coil, 2: helix, 3: turn, 4: strand)

* #### con.lst

con.lst file provides a list of atomic and residue contacts between receptor chain(s) and ligand chain(s), highlighting the total atomic contacts and residue-residue contacts with specific interacting residues and their contact counts based on a 4.50 Å distance cutoff. An index is assigned to each residue.

### optional output files

* #### parsed.pdb

parsed.pdb contains atomic information of all atoms of all residues in the selected ligand and receptor chains.

* #### vmd

vmd file contains a TCL script for visualizations of an interface with VMD. The molecule visualization styles is NewCartoon.

* #### trans.pdb

trans.pdb is a variant of the first int_pdb file where the atomic coordinates are modified by superposing them onto int_pdb2 using the transformation matrix that generates the optimal interface alignment.

## Stdout

Stdout shows the alignment result of a pair of protein-protein interfaces. It displays, IS-score, P-value, Z-score, Number of aligned residues between interfaces, Number of aligned contacts, RMSD, and Seq identity. "Aligned Interface Residues" section has 12 column values. "Index" has indices of aligned interface residue pairs. "Ch1", "Resid1" and "AA1" are chainIDs, residue IDs, and 3 letter abbreviations of corresponding amino acids for the interface of the first PDB file. "Ch2", "Resid2" and "AA2" are for the interface of the second PDB file. "Distance" provides the distance between interface residues' alpha carbons in angstrom. "NAC" provides the number of atomic contact. "NC1" and "NC2" gives the number of contact residues within interface in the first and second PDB files, respectively. "Note" denotes whether the interface residue pairs are within 5 Ansgtrom and residues are identical.

## Acknowledgment

PiAlign is built upon the pre-existing tool named iAlign, which was developed by Dr. Mu Gao in 2010. We are grateful for the guidance and support from Dr. Gao during the development of this tool. Additionally, this tool utilizes the free and open-source Python package [ProDy](http://prody.csb.pitt.edu/), developed by Bahar Lab at the Laufer Center, Stony Brook University.

## Reference

- Gao, M., & Skolnick, J. (2010). iAlign: A method for the structural comparison of protein–protein interfaces. Bioinformatics, 26(18), 2259–2265. https://doi.org/10.1093/bioinformatics/btq404

- Gao, M., & Skolnick, J. (2011). New benchmark metrics for protein-protein docking methods. Proteins, 79(5), 1623–1634. https://doi.org/10.1002/prot.22987

- Zhang, S., Krieger, J. M., Zhang, Y., Kaya, C., Kaynak, B., Mikulska-Ruminska, K., Doruker, P., Li, H., & Bahar, I. (2021). ProDy 2.0: Increased scale and scope after 10 years of protein dynamics modelling with Python. Bioinformatics, 37(20), 3657–3659. https://doi.org/10.1093/bioinformatics/btab187

## License

PiAlign is available under MIT License. See LICENSE.rst for more details.
