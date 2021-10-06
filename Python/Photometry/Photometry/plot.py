import numpy as np
import matplotlib.pyplot as plt

from database import Database
from point import Point
from config import n_pixel

def plot_points(b_observed_points=True,b_simulated_points=True):

    db = Database()
    fig = plt.figure()
    fig.set_size_inches(9,9)

    if(b_observed_points):
        observed_points = db.select_points(0, True)
        op_arr = np.vstack(observed_points[:]).astype(float)
        plt.scatter(op_arr[:,0], op_arr[:,1], s=1, c='r', marker="s", label='observed')
    if(b_simulated_points):
        simulated_points = db.select_points(0, False)
        sp_arr = np.vstack(simulated_points[:]).astype(float)
        plt.scatter(sp_arr[:,0], sp_arr[:,1], s=1, c='g', marker="o", label='simulated')

    plt.xlabel('ascension [arcsec]', fontsize=16)
    plt.ylabel('declination [arcsec]', fontsize=16)
    plt.legend(loc='upper left')
    plt.show()


def plot_magnitude_hist():
    db = Database()
    fig = plt.figure()
    fig.set_size_inches(9,9)

    observed_points = db.select_points(0, True)
    op_arr = np.vstack(observed_points[:]).astype(float)

    simulated_points = db.select_points(0, False)
    sp_arr = np.vstack(simulated_points[:]).astype(float)

    plt.hist(op_arr[:,5]-12.5, density=True, bins=50, alpha=0.5, label='observed')
    plt.hist(sp_arr[:,5], density=True, bins=50, alpha=0.5, label='simulated')

    plt.legend(loc='upper left')
    plt.show()


def plot_points_velocity(b_observed_points=True,b_simulated_points=True):

    db = Database()
    fig = plt.figure()
    fig.set_size_inches(9,9)
    #print(maxDist[0],maxDist[1])
    #plotDist = 0.00002345 #np.minimum(maxDist[0],maxDist[1])
    #
    #delta = 0.0075
    #plt.xlim(20.4921513908745-delta, 20.4921513908745+delta)
    #plt.ylim(-5.09226524786853-delta, -5.09226524786853+delta)

    if(b_observed_points):
        observed_points = db.select_points(0, True)
        op_arr = np.vstack(observed_points[:]).astype(float)
        op_cluster = op_arr[op_arr[:,4] > -1]
        plt.scatter(op_cluster[:,2], op_cluster[:,3], s=1, c='r', marker="s", label='observed_cluster')
        op_fs = op_arr[op_arr[:,4] == -1]
        plt.scatter(op_fs[:,2], op_fs[:,3], s=1, c='black', marker="s", label='observed_fs')
    if(b_simulated_points):
        simulated_points = db.select_points(0, False)
        sp_arr = np.vstack(simulated_points[:]).astype(float)
        sp_cluster = sp_arr[sp_arr[:,4] > -1]
        plt.scatter(sp_cluster[:,2], sp_cluster[:,3], s=1, c='g', marker="s", label='simulated_cluster')
        sp_fs = sp_arr[sp_arr[:,4] == -1]
        plt.scatter(sp_fs[:,2], sp_fs[:,3], s=1, c='blue', marker="s", label='simulated_fs')

    plt.xlabel('v_asc [arcsec/dt]', fontsize=16)
    plt.ylabel('v_dec [arcsec/dt]', fontsize=16)
    plt.legend(loc='upper left')
    plt.show()

#def plot_points_velocity(points):
#    arr = np.vstack(points[:]).astype(float)
 #   fig = plt.figure()
#
#    fig.set_size_inches(9,9)
#
#    plt.scatter(arr[:,2], arr[:,3], s=1, c='r', marker="s")
#
#    plt.xlabel('v_asc [arcsec/dt]', fontsize=16)
#    plt.ylabel('v_dec [arcsec/dt]', fontsize=16)
#    plt.show()

def plot_map():
    masses =np.array([640,1600,4000,10000,25000])
    x = [a+0.5 for a in range(len(masses))]
    angles = np.array([180,25,10,5,0])
    y = [a+0.5 for a in range(len(angles))]
    #z = [0.999831366,0.999470002,0.999053414,0.99863088,0.993632332,0.961191213,0.965806713,0.975989894,0.980775459,0.985457224,0.855673786,
    #    0.873354153,0.9015785,0.92996853,0.950666177,0.926183906,0.939546943,0.956671385,0.963135931,0.972534242,0.625338221,0.697084463,0.771201291,
    #    0.823656487,0.866523918] #total
    z = [0.999755501,0.999253485,0.998640674,0.998111424,0.992187177,0.953405272,0.960113393,0.971606404,0.97693511,0.982810825,0.808736986,
         0.846920703,0.878962645,0.916749616,0.942218371,0.870929796,0.912485069,0.938408026,0.951847771,0.964182264,0.486505394,0.612119595,0.698490662,0.771850264,0.821191887] #C P 0.5 - 0.08
    z = np.reshape(z, (5, 5))
    #plt.xscale('log')
    #plt.ylim(0,60)
    plt.xlabel('cluster mass [$M_{\odot}$]', fontsize=16)
    plt.ylabel('angle [$^\circ$]', fontsize=16)
    plt.title('confidence 0.5 - 0.08 [$M_{\odot}$]')
    plt.pcolor(z)
    plt.xticks(x, masses, fontsize=16)
    plt.yticks(y, angles, fontsize=16)
    #plt.contourf(x,y,z)
    plt.colorbar()
    plt.show();

def main():
    db = Database()
    plot_map()
    #print(len(observed_points))
    #print(len(simulated_points))

    #op_arr = np.vstack(observed_points[:]).astype(float)
    #sp_arr = np.vstack(simulated_points[:]).astype(float)
    #plot_points()
    #plot_magnitude_hist()
    #plot_points_velocity(True,True)
    #plot_points(sp_arr)
    #plot_magnitude_hist(op_arr,sp_arr)
    #print(float_arr)

    #observed_t0 = db.select_observed_cluster_stars(0)
    #observed_t1 = db.select_observed_cluster_stars(1)
    #observed_t0_arr = np.vstack(observed_t0[:]).astype(float)
    #observed_t1_arr = np.vstack(observed_t1[:]).astype(float)

    #simulated_t0 = db.select_cluster(0,False)
    #simulated_t1 = db.select_cluster(1,False)
    #simulated_t0_arr = np.vstack(observed_t0[:]).astype(float)
    #simulated_t1_arr = np.vstack(observed_t1[:]).astype(float)

    #fig = plt.figure()

    #fig.set_size_inches(9,9)

    #plt.scatter(observed_t0_arr[0,0], observed_t0_arr[0,1], s=2, c='r', marker="s", label='observed_t0')
    #plt.scatter(observed_t1_arr[:,0], observed_t1_arr[:,1], s=2, c='g', marker="o", label='observed_t1')

    #plt.scatter(simulated_t0_arr[0,0], simulated_t0_arr[0,1], s=1, c='b', marker="x", label='simulated_t0')
    #plt.scatter(simulated_t1_arr[:,0], simulated_t1_arr[:,1], s=1, c='black', marker="+", label='simulated_t1')

    #data = db.select_distace()

    #distance = print(sum(data[:,0] - data[:,2]),sum(data[:,1] - data[:,3]))

    #plt.xlabel('ascension [arcsec]', fontsize=16)
    #plt.ylabel('declination [arcsec]', fontsize=16)
    #plt.legend(loc='upper left')
    #plt.show()

    #plot_points_velocity(points)


if __name__ == '__main__':
    main()