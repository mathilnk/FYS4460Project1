
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
#filename = "unitcell2.xyz";
filename = sys.argv[1] + ".xyz";
infile = open(filename, 'r');
lines = infile.readlines();
length = eval(lines[0]);
vx = np.zeros(length);
vy = np.zeros(length);
vz = np.zeros(length);
v = np.zeros(length);


for i in range(2,length):
    line = lines[i];
    words = line.split()
    vx[i] = eval(words[4])
    vy[i] = eval(words[5])
    vz[i] = eval(words[6])
    v[i] = np.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2);
av_vx = np.sum(vx)/length;
av_vy = np.sum(vy)/length
av_vz = np.sum(vz)/length
av_v = np.sum(v)/length
print av_vx, av_vy, av_vz, av_v
s_vx = np.sqrt(np.sum(vx**2)/length - av_vx**2);
s_vy = np.sqrt(np.sum(vy**2)/length - av_vy**2);
s_vz = np.sqrt(np.sum(vz**2)/length - av_vz**2);
print s_vx,s_vy, s_vz 




# the histogram of the data
#nv, binsv, patchesv = plt.hist(v, 50, normed=1, facecolor='green', alpha=0.75, histtype='step')
#nx, binsx, patchesx = plt.hist(vx, 50, normed=1, facecolor='blue', alpha=0.75, histtype='step')
#ny, binsy, patchesy = plt.hist(vy, 50, normed=1, facecolor='red', alpha=0.75, histtype='step')
nz, binsz, patchesz = plt.hist(vz, 50, normed=1, facecolor='black', alpha=0.75, histtype='step')



plt.xlabel('velocity')
plt.ylabel('Distribution')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
#plt.grid(True)
plt.savefig(sys.argv[1]+".png");

plt.show()
