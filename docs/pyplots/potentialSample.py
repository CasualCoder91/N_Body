import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

dataDisk = np.loadtxt('potentialDiskPositionsSample5000.dat',delimiter=',')
dataBulge = np.loadtxt('potentialBulgePositionsSample5000.dat',delimiter=',')
fig = plt.figure()
fig.set_size_inches(14, 4.2)
ax = fig.add_subplot(2,2,1)
plt.hist(dataDisk[:,0]+dataBulge[:,0], bins=200)
ax.title.set_text('Distribution along x Axis')

ax = fig.add_subplot(2,2,2)
ax.title.set_text('Distribution along z Axis')
plt.hist(dataDisk[:,2]+dataBulge[:,2], bins=200)


ax = fig.add_subplot(2,2,(3,4), projection='3d')
ax.scatter(dataDisk[:,0], dataDisk[:,1], dataDisk[:,2], s=5,c='red')
ax.scatter(dataBulge[:,0], dataBulge[:,1], dataBulge[:,2], s=5,c='blue')
ax.view_init(elev=5., azim=180)
ax.dist = 8
ax.set_xlim([-40,40])
ax.set_ylim([-40,40])
ax.set_zlim([-7,7])

#plt.autoscale(False)
plt.show()
