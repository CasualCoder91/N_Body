import numpy as np
import matplotlib.pyplot as plt

def bulgeDispersion(dataPath='',showPlot=True,arguments=[]):
    dataMeasured = np.loadtxt(dataPath + '/bulgeVelocitiesMeasured.dat',delimiter=',')
    data4 = np.loadtxt(dataPath + '/bulgeDispersion-4.dat',delimiter=',')
    data6 = np.loadtxt(dataPath + '/bulgeDispersion-6.dat',delimiter=',')
    data8 = np.loadtxt(dataPath + '/bulgeDispersion-8.dat',delimiter=',')
    fig = plt.figure()
    plt.plot(data4[:,0], data4[:,1], 'r-',label='b = -4[°],d=8.5kpc')
    plt.draw()
    plt.plot(data6[:,0], data6[:,1], 'g-',label='b = -6[°],d=8.5kpc')
    plt.draw()
    plt.plot(data8[:,0], data8[:,1], 'b-',label='b = -8[°],d=8.5kpc')
    plt.draw()

    fig.ylabel = "radial dispersion km/s"
    fig.xlabel = "longitute l [°]"
    plt.title("Radial velocity dispersion bulge")

    for data in dataMeasured:
        if data[1] ==-4:
            plt.plot(data[0], data[2], 'ro')
            plt.draw()
            plt.errorbar(data[0], data[2],yerr=data[3],fmt='r-o')
        if data[1] ==-6:
            plt.plot(data[0], data[2], 'go')
            plt.draw()
            plt.errorbar(data[0], data[2],yerr=data[3],fmt='g-o')
        if data[1] ==-8:
            plt.plot(data[0], data[2], 'bo')
            plt.draw()
            plt.errorbar(data[0], data[2],yerr=data[3],fmt='b-o')
    plt.legend()

    if showPlot:
        plt.show()

if __name__ == "__main__":
    bulgeDispersion()
