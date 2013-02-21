import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys

project_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460Project1/results/"

def plot_energy():
	plt.figure()
	k = plt.plot(kin_energies);
	p = plt.plot(pot_energies);
	t = plt.plot(tot_energies);
	plt.xlabel('Time steps');
	plt.ylabel('Energy');
	plt.title('Energy - dimensionless');
	plt.legend(["Kinetic Energy", "Potential Energy", "Total Energy"]);
	average_energy = sum(tot_energies)/numTimeSteps;
	print "Average energy: ", average_energy;
	return 0;
def plot_temp():
	plt.figure()
	plt.plot(temperature);
	plt.xlabel('Time steps');
	plt.ylabel('Temperature');
	plt.title('Temperature - dimensionless');
	average_temp = sum(temperature)/numTimeSteps;
	print "Average temperature: ", average_temp;
	return 0;
def plot_pressure():
	plt.figure()
	plt.plot(pressure);
	plt.xlabel('Time steps');
	plt.ylabel('Pressure');
	plt.title('Pressure - dimensionless');
	average_pressure = sum(pressure)/numTimeSteps;
	print "Average pressure ", average_pressure;
	return 0;
def plot_pressure_T():
	plt.figure()
	plt.plot(temperature[100:], pressure[100:], '.');
	plt.xlabel('Temperature');
	plt.ylabel('Pressure');
	plt.title('Pressure as a function of Temperature - dimensionless');
	average_pressure = sum(pressure)/numTimeSteps;
	print "Average pressure ", average_pressure;
	return 0;
def plot_displacement():
	plt.figure()
	plt.plot(displacement);
	plt.xlabel('Time steps');
	plt.ylabel('Average displacement');
	plt.title('The average displacement - dimensionless');
	diffusion_constant = (displacement[-1]-displacement[0])/(6*numTimeSteps*dt);
	print "Total displacement: ", displacement[-1];
	print "Diffusion constant: ", diffusion_constant;
	return 0;
def plot_radial():
	plt.figure()
	plt.plot(radial_dist);
	plt.xlabel('r');
	plt.ylabel('Distribution');
	plt.title("The radial distribution function g(r)");
	return 0;

if(len(sys.argv)>3):
	filenamebase = sys.argv[1];
        plot_type = sys.argv[3]; 
	#special case
	#average = 26488.18966;
	numOfFiles = eval(sys.argv[2]);
	average_energies = np.zeros(numOfFiles);	
	for i in range(numOfFiles):
		fluct = 0;
		filename = project_dir + filenamebase +  "_energy.txt";
		infile = open(filename, 'r');
		lines = (infile.readlines());

		numTimeSteps = len(lines)-2;
		tot_energies = np.zeros(numTimeSteps);
		kin_energies = np.zeros(numTimeSteps);
		pot_energies = np.zeros(numTimeSteps);

		temperature = np.zeros(numTimeSteps);
		pressure = np.zeros(numTimeSteps);
		displacement = np.zeros(numTimeSteps);
		
		dt = eval(lines[0].split()[1]);
		print dt;
		for j in range(2,numTimeSteps):
			measurements = lines[j].split();
			kin_energies[j] = measurements[0];
			pot_energies[j] = measurements[1];
			tot_energies[j] = measurements[2];
			temperature[j] = measurements[3];
			pressure[j] = measurements[4];
			displacement[j] = measurements[5];

		#fluct/=numOfEnergies;
		if(plot_type == "E"):
			plot_energy();
		elif(plot_type == "T"):
			plot_temp();
		elif(plot_type == "P"):
			plot_pressure();
		elif(plot_type == "P_T"):
			plot_pressure_T();
		elif(plot_type == "D"):
			plot_displacement();
		elif(plot_type == "all"):
			plot_energy();
			plot_temp();
			plot_pressure();
			plot_pressure_T();
			plot_displacement();
			
		

	#plt.axis([0,5000,0,31000])
	plt.show();
else:
	print "You forgot filname, number of files, and plottype"
