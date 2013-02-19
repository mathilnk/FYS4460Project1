import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys

project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460Project1/project1-build-desktop-Qt_4_8_1_in_PATH__System__Release/"

if(len(sys.argv)>2):
	filenamebase = sys.argv[1];
	numOfFiles = eval(sys.argv[2]);
	average_energies = np.zeros(numOfFiles);	
	for i in range(numOfFiles):
		filename = project_dir + filenamebase + "%d" %(i) + "_energy.txt";
		infile = open(filename, 'r');
		energies = (infile.readlines());
		numOfEnergies = len(energies);
		pl_energies = zeros(numOfEnergies);
		for j in range(numOfEnergies):
			energy = eval(energies[i]);
			pl_energies[j] = energy;
			average_energies[i]+=energy;
		average_energies[i] = average_energies[i]/numOfEnergies;
		plt.plot(pl_energies);	
		
	print energies;
