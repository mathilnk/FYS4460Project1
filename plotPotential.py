import numpy as np;
import matplotlib.pyplot as mpl

r = np.linspace(0.999, 3,100);
U = 4*(1.0/r**12 - 1.0/r**6)
mpl.plot(r,U);
mpl.xlabel('Distance between particles')
mpl.ylabel('Potential');
mpl.title('The Lennard-Jones potential')
mpl.axis([0.999,3,-1,0.2])
mpl.show()
