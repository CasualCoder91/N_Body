import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

data = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/potentialVelocity.dat',delimiter=',')
dataDisk = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/potentialVelocityDisk.dat',delimiter=',')
dataBlackHole = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/potentialVelocityBlackHole.dat',delimiter=',')
dataBuldge = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/potentialVelocityBuldge.dat',delimiter=',')
dataHalo = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/potentialVelocityHalo.dat',delimiter=',')
fig, ax = plt.subplots()
ax.plot(data[:,0]*0.001, data[:,1], label='total')
ax.plot(dataDisk[:,0]*0.001, dataDisk[:,1], label='disk')
ax.plot(dataBlackHole[:,0]*0.001, dataBlackHole[:,1], label='black hole')
ax.plot(dataBuldge[:,0]*0.001, dataBuldge[:,1], label='buldge')
ax.plot(dataHalo[:,0]*0.001, dataHalo[:,1], label='halo')
plt.xlabel('distance [kpc]')
plt.ylabel('circular velocity [km/s]')
plt.title('Milky Way Potential: Circular Velocity ')
ax.legend()
plt.show()
