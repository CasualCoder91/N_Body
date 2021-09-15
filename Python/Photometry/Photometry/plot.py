import numpy as np
import matplotlib.pyplot as plt

from database import Database
from point import Point
from config import n_pixel

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

    plt.scatter(observed_points[:,0], observed_points[:,1], s=1, c=observed_points[:,4], marker="s", label='observed')
    plt.scatter(simulated_points[:,0], simulated_points[:,1], s=1, c=simulated_points[:,4], marker="o", label='simulated')

    plt.xlabel('ascension [arcsec]', fontsize=16)
    plt.ylabel('declination [arcsec]', fontsize=16)
    plt.legend(loc='upper left')
    plt.show()
    #name = output+'\starPositions'+str(int(i))

    #plt.savefig(name+'.jpg')
    #plt.close(fig)

def plot_magnitude_hist(observed_points,simulated_points):
    fig = plt.figure()
    fig.set_size_inches(9,9)

    plt.hist(observed_points[:,5]-15.5, density=True, bins=50, alpha=0.5, label='observed')
    plt.hist(simulated_points[:,5], density=True, bins=50, alpha=0.5, label='simulated')

    plt.legend(loc='upper left')
    plt.show()


def plot_points_velocity(observed_points,simulated_points):
    fig = plt.figure()

    fig.set_size_inches(9,9)
    #print(maxDist[0],maxDist[1])
    #plotDist = 0.00002345 #np.minimum(maxDist[0],maxDist[1])
    #
    #delta = 0.0075
    #plt.xlim(20.4921513908745-delta, 20.4921513908745+delta)
    #plt.ylim(-5.09226524786853-delta, -5.09226524786853+delta)

    plt.scatter(observed_points[:,2], observed_points[:,3], s=1, c='r', marker="s", label='observed')
    plt.scatter(simulated_points[:,2], simulated_points[:,3], s=1, c='g', marker="o", label='simulated')

    plt.xlabel('v_asc [arcsec/dt]', fontsize=16)
    plt.ylabel('v_dec [arcsec/dt]', fontsize=16)
    plt.legend(loc='upper left')
    plt.show()
    #name = output+'\starPositions'+str(int(i))

    #plt.savefig(name+'.jpg')
    #plt.close(fig)


def main():
    db = Database()
    #observed_points = db.select_points(0, True)
    #simulated_points = db.select_points(0, False)
    observed_points = db.select_cluster(0, True)
    simulated_points = db.select_cluster(0, False)

    #print(len(observed_points))
    #print(len(simulated_points))

    op_arr = np.vstack(observed_points[:]).astype(float)
    sp_arr = np.vstack(simulated_points[:]).astype(float)

    plot_points_velocity(op_arr,sp_arr)
    #plot_points(op_arr,sp_arr)
    #plot_magnitude_hist(op_arr,sp_arr)
    #print(float_arr)

if __name__ == '__main__':
    main()