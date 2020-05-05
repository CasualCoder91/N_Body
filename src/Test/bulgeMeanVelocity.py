import numpy as np
import matplotlib.pyplot as plt

#dataMeasured = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/bulgeVelocitiesMeasured.dat',delimiter=',')
data4 = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/bulgeMeanVelocity-4.dat',delimiter=',')
data6 = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/bulgeMeanVelocity-6.dat',delimiter=',')
data8 = np.loadtxt('E:/Master_Thesis/VS_Project/N_Body/src/Test/bulgeMeanVelocity-8.dat',delimiter=',')
plt.figure()
plt.plot(data4[:,0], data4[:,1], 'r-',label='b = -4[°],d=8.5kpc')
plt.plot(data6[:,0], data6[:,1], 'g-',label='b = -6[°],d=8.5kpc')
plt.plot(data8[:,0], data8[:,1], 'b-',label='b = -8[°],d=8.5kpc')
plt.ylabel("radial dispersion km/s")
plt.xlabel("longitute l [°]")
plt.title("Mean velocity bulge")

#for data in dataMeasured:
#    if data[1] ==-4:
#        plt.plot(data[0], data[2], 'ro')
#        plt.errorbar(data[0], data[2],yerr=data[3],fmt='r-o')
#    if data[1] ==-6:
#        plt.plot(data[0], data[2], 'go')
#        plt.errorbar(data[0], data[2],yerr=data[3],fmt='g-o')
#    if data[1] ==-8:
#        plt.plot(data[0], data[2], 'bo')
#        plt.errorbar(data[0], data[2],yerr=data[3],fmt='b-o')
plt.legend()
plt.show()
