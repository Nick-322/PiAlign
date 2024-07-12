Two simple examples illustrating usage of PiAlign

i. Align two protein-protein interfaces (one chain for each side of the interface)

../bin/PiAlign.py -p1 1lyl.pdb -c1a A -c1b C -p2 12as.pdb -c2a A -c2b B



ii. Align two protein-protein interfaces with a multiple interface chain case where we don't know which chains interact with which chains for the first pdb file.

../bin/PiAlign.py -p1 7e5s.pdb -c1a LEIJPRUVDHKTNOQS -c1b BCA -p2 7uap.pdb -c2a HL -c2b A -searchIntCh1

