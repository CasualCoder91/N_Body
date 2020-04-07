import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

data = np.loadtxt('frequencyDistribution_z1.dat',delimiter=',')
x = data[:,0]
y = data[:,1]
z = data[:,2]
N = int(len(z)**.5)
z = z.reshape(N, N)
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.title('z = 1.0 [kpc] | [Color] = [$M_{\odot}$]')
plt.imshow(z+10, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
        cmap=cm.hot, norm=LogNorm())
plt.colorbar()
plt.show()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(data[:,0], data[:,1], data[:,2], s=10,c='red')
