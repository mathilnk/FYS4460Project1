import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
def plot_radial():
	plt.figure()
	plt.plot(radial_dist);
	plt.xlabel('r');
	plt.ylabel('Distribution');
	plt.title("The radial distribution function g(r)");
	return 0;
	
	
project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460Project1/project1-build-desktop-Qt_4_8_1_in_PATH__System__Release/"

if(len(sys.argv)>1):
	filename = sys.argv[1];
	infile = open(project_dir + filename, 'r');
	lines = infile.readlines()
	N = len(lines);
	radial_dist = np.zeros(N);
	for i in range(N):
		radial_dist[i] = eval(lines[i]);
	plot_radial();
	plt.<<show();
