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

    plt.scatter(observed_points[:,0], observed_points[:,1], s=1, c='b', marker="s", label='observed')
    plt.scatter(simulated_points[:,0], simulated_points[:,1], s=1, c='r', marker="o", label='simulated')

    plt.xlabel('ascension [arcsec]', fontsize=16)
    plt.ylabel('declination [arcsec]', fontsize=16)
    plt.legend(loc='upper left')
    #name = output+'\starPositions'+str(int(i))

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

    avg_asc_obs = np.median(op_arr[:,0])
    avg_dec_obs = np.median(op_arr[:,1])

    avg_asc_sim = np.median(sp_arr[:,0])
    avg_dec_sim = np.median(sp_arr[:,1])

    print("avg ascension observed: ", avg_asc_obs)
    print("avg dec observed: ", avg_dec_obs)

    print("avg ascension simulated: ", avg_asc_sim)
    print("avg dec simulated: ", avg_dec_sim)

    fig = plt.figure()

    fig.set_size_inches(9,9)
    #print(maxDist[0],maxDist[1])
    #plotDist = 0.00002345 #np.minimum(maxDist[0],maxDist[1])
    #
    #delta = 0.0075
    #plt.xlim(20.4921513908745-delta, 20.4921513908745+delta)
    #plt.ylim(-5.09226524786853-delta, -5.09226524786853+delta)

    plt.scatter(op_arr[:,0], op_arr[:,1], s=1, c='b', marker="s", label='observed')
    plt.scatter(sp_arr[:,0]-0.002, sp_arr[:,1]-0.002, s=1, c='r', marker="o", label='simulated')

    plt.xlabel('ascension [arcsec]', fontsize=16)
    plt.ylabel('declination [arcsec]', fontsize=16)
    plt.legend(loc='upper left')

    #plt.hist(op_arr[:,1],bins=100)
    #plt.hist(sp_arr[:,1],bins=100)
    #plt.axvline(np.median(op_arr[:,1]))
    #plt.axvline(np.median(sp_arr[:,1]))

    plt.scatter(avg_asc_obs,avg_dec_obs,s=100,label='obs_avg',zorder=100)
    plt.scatter(avg_asc_sim,avg_dec_sim,s=100,label='sim_avg',zorder=100)

    plt.show()

    #print(float_arr)

if __name__ == '__main__':
    main()