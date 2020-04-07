import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

data = np.loadtxt('potentialDiskSample.dat',delimiter=',')
data2 = np.loadtxt('potentialBulgeSample.dat',delimiter=',')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data[:,0], data[:,1], data[:,2], s=5,c='red')
ax.scatter(data2[:,0], data2[:,1], data2[:,2], s=5,c='blue')
ax.view_init(elev=5., azim=180)
ax.dist = 5
plt.show()
