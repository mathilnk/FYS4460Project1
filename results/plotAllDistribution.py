import sys, os

filenamebase = sys.argv[1];
number_of_files = eval(sys.argv[2]);

#plot the figures
for i in range(number_of_files):
	line = "python plot_distribution_b.py "+ filenamebase + "%d"%i;
	os.system(line);

#make movie of the figures
line2 = "python makeMovie.py " + sys.argv[1] + " png rm";
os.system(line2);
