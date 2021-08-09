import numpy as np
import matplotlib.pyplot as plt

from database import Database
from point import Point

def plot_points(observed_points,simulated_points):
    #loop timesteps
    fig = plt.figure()

    fig.set_size_inches(9,9)
    #print(maxDist[0],maxDist[1])
    #plotDist = 0.00002345 #np.minimum(maxDist[0],maxDist[1])
    #
    #delta = 0.0075
    #plt.xlim(20.4921513908745-delta, 20.4921513908745+delta)
    #plt.ylim(-5.09226524786853-delta, -5.09226524786853+delta)
    #ax1 = fig.add_subplot(111)

    plt.legend(loc='upper left')
    plt.scatter(observed_points[:,0], observed_points[:,1], s=1, c='b', marker="s", label='observed')
    plt.scatter(simulated_points[:,0], simulated_points[:,1], s=1, c='r', marker="o", label='simulated')

    #plt.scatter(observed_points[:,0], observed_points[:,1], c='r',s=1)#c=timestepData[:,5]
    #plt.scatter(simulated_points[:,0], simulated_points[:,1], c='b',s=1)#c=timestepData[:,5]
    plt.xlabel('ascension [arcsec]', fontsize=16)
    plt.ylabel('declination [arcsec]', fontsize=16)
    #print(plotDist)
    #name = output+'\starPositions'+str(int(i))
    plt.show()
    #plt.savefig(name+'.jpg')
    #plt.close(fig)

def main():
    db = Database()
    observed_points = db.select_points(0, True)
    simulated_points = db.select_points(0, False)

    print(len(observed_points))
    print(len(simulated_points))

    op_arr = np.vstack(observed_points[:]).astype(float)
    sp_arr = np.vstack(simulated_points[:]).astype(float)

    plot_points(op_arr,sp_arr)
    #print(float_arr)

if __name__ == '__main__':
    main()