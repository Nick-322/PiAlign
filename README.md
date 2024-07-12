## PiAlign: a tool for the flexible structural comparison of protein–protein interfaces free of chain number restriction

Proteins play fundamental roles in the living cells of all organisms, responsible for a diverse array of biological functions. ​To elucidate protein functionality, it is crucial to understand their structural interactions, particularly the regions where physical contact occurs between proteins, known as protein-protein interfaces. Comparative analysis of these interfaces facilitates the understanding of biological and evolutionary relationships among proteins. In 2010, iAlign was developed to enable efficient, length-independent structural comparison of interfaces formed by two protein chains. Building upon this foundation, PiAlign is a computational tool that extends iAlign's capabilities by allowing structural comparison of protein-protein interfaces formed by any number of chains, not limited to just two. This expansion of the interface analysis scope enables users to capture more comprehensive structural relationships between large protein complexes, thereby enhancing our understanding of protein interactions and their biological implications.


## Installation

#### First, clone this repository by

	git clone https://github.com/

#### To use PiAlign, you will have two choices to create a virtual enviroment with required dependencies.

* #### Choice 1. Conda

	*If conda is not installed in your computer, refer to the guide from an official site:
	https://conda.io/projects/conda/en/latest/user-guide/install/index.html


	Once codnda is installed, create a enironment by running the following code:
	* For MacOS (ARM64) user:
		> CONDA_SUBDIR=osx-64 conda env create -f <PATH_TO>/PiAlign/environment.yml


	* For Linux (x86-64) user:
		> conda env create -f environment.yml



* #### Choice 2. Python venv
	* Create a enironment by running the following code:

			python -m venv <environment_name>


	* And run the code below:

			source <environment_name>/bin/activate


	* Lastly, run:

			pip install -r requirements.txt



#### In the created environment, you may set the following alias:

	alias pialign='<PATH_TO>/PiAlign/bin/PiAlign.py'

#### and invoke PiAlign by

  	pialign



#### To run PiAlign successfuly, you will need to install dssp in your environment. ###
* For Linux user: run 

		sudo apt-get install dssp

* For Mac user:
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

i. Align two protein-protein interfaces (one chain for each side of the interface)

../bin/PiAlign.py -p1 1lyl.pdb -c1a A -c1b C -p2 12as.pdb -c2a A -c2b B



ii. Align two protein-protein interfaces with a multiple interface chain case where we don't know which chains interact with which chains for the first pdb file.

../bin/PiAlign.py -p1 7e5s.pdb -c1a LEIJPRUVDHKTNOQS -c1b BCA -p2 7uap.pdb -c2a HL -c2b A -searchIntCh1



## Output files


## Reference

- Gao, M., & Skolnick, J. (2010). iAlign: A method for the structural comparison of protein–protein interfaces. Bioinformatics, 26(18), 2259–2265. https://doi.org/10.1093/bioinformatics/btq404

- Gao, M., & Skolnick, J. (2011). New benchmark metrics for protein-protein docking methods. Proteins, 79(5), 1623–1634. https://doi.org/10.1002/prot.22987

- Zhang, S., Krieger, J. M., Zhang, Y., Kaya, C., Kaynak, B., Mikulska-Ruminska, K., Doruker, P., Li, H., & Bahar, I. (2021). ProDy 2.0: Increased scale and scope after 10 years of protein dynamics modelling with Python. Bioinformatics, 37(20), 3657–3659. https://doi.org/10.1093/bioinformatics/btab187


## License
