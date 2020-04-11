import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

data = np.loadtxt('initialConditionsSurfaceDensityBulge.dat',delimiter=',')

plt.semilogx(data[:,0], data[:,1])

#plt.xlabel('distance [kpc]')
#plt.ylabel('circular velocity [km/s]')
#plt.title('Milky Way Potential: Circular Velocity ')
#ax.legend()
plt.show()
