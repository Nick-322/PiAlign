# To use PiAlign, you will have two choices to create a virtual enviroment with required dependencies. 

1. Conda
If conda is not installed in your computer, refer to the guide from an official site:
https://conda.io/projects/conda/en/latest/user-guide/install/index.html


Once codnda is installed, create a enironment by running the following code:
* For MacOS (ARM64) user:
	> CONDA_SUBDIR=osx-64 conda env create -f <PATH_TO>/PiAlign/environment.yml


* For Linux (x86-64) user:
	> conda env create -f environment.yml



2. Python venv
Create a enironment by running the following code:
	> python -m venv <environment_name>


and run the code below
	> source <environment_name>/bin/activate


lastly, 
	> pip install -r requirements.txt



# In the created environment, you may set the following alias:

	> alias pialign='<PATH_TO>/PiAlign/bin/PiAlign.py'

and invoke PiAlign by

  	> pialign



# To run PiAlign successfuly, you will need to install dssp in your environment. ###
* For Linux user: type "sudo apt-get install dssp"

* For Mac user:
1. Install MacPorts
    a. Install Apple's Command Line Developer Tools: xcode-select --install
    b. Install MacPorts for your version of the Mac operating system:
        macOS Sonoma v14
        macOS Ventura v13
        macOS Monterey v12
        macOS Big Sur v11

2. Run "sudo port install dssp"
