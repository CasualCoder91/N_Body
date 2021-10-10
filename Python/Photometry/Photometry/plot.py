import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from database import Database
from point import Point
import config

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

def plot_cluster(b_observed_points=True,b_simulated_points=True):
    db = Database()
    fig = plt.figure()
    fig.set_size_inches(9,9)

    if(b_observed_points):
        observed_points = db.select_cluster(0,True)
        op_arr = np.vstack(observed_points[:]).astype(float)
        plt.scatter(op_arr[:,0], op_arr[:,1], s=1, c='r', marker="s", label='observed')
    if(b_simulated_points):
        simulated_points = db.select_cluster(0,False)
        sp_arr = np.vstack(simulated_points[:]).astype(float)
        plt.scatter(sp_arr[:,0], sp_arr[:,1], s=1, c='g', marker="o", label='simulated')

    false_negative = db.select_false_negative()
    fn_arr = np.vstack(false_negative[:]).astype(float)
    plt.scatter(fn_arr[:,0], fn_arr[:,1], s=1, c='b', marker="x", label='false negative')

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


def plot_points_velocity(b_observed_points=True,b_simulated_points=True,b_false_negative=False):

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
        #plt.scatter(op_fs[:,2], op_fs[:,3], s=1, c='black', marker="s", label='observed_fs')
    if(b_simulated_points):
        simulated_points = db.select_points(0, False)
        sp_arr = np.vstack(simulated_points[:]).astype(float)
        sp_cluster = sp_arr[sp_arr[:,4] > -1]
        plt.scatter(sp_cluster[:,2], sp_cluster[:,3], s=1, c='g', marker="s", label='simulated_cluster')
        sp_fs = sp_arr[sp_arr[:,4] == -1]
        plt.scatter(sp_fs[:,2], sp_fs[:,3], s=1, c='blue', marker="s", label='simulated_fs')
    if(b_false_negative):
        false_negative = db.select_false_negative()
        fn_arr = np.vstack(false_negative[:]).astype(float)
        plt.scatter(fn_arr[:,2], fn_arr[:,3], s=1, c='yellow', marker="x", label='false negative')

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

def plot_map(z,title,cmap):
    masses =np.array([640,1600,4000,10000,25000])
    x = [a+0.5 for a in range(len(masses))]
    angles = np.array([180,25,10,5,0])
    y = [a+0.5 for a in range(len(angles))]
    z = np.reshape(z, (5, 5))
    #plt.xscale('log')
    #plt.ylim(0,60)
    plt.clf()
    plt.xlabel('cluster mass [$M_{\odot}$]', fontsize=16)
    plt.ylabel('angle [$^\circ$]', fontsize=16)
    plt.title(title+' [$M_{\odot}$]')
    plt.pcolor(z, cmap=cmap)
    plt.xticks(x, masses, fontsize=16)
    plt.yticks(y, angles, fontsize=16)
    #plt.contourf(x,y,z)
    plt.colorbar()
    fig = plt.gcf()
    fig.set_size_inches(9, 9)
    fig.savefig(config.output_base_path+'\\Clustering\\'+title+'.png', dpi=100)
    #plt.show();

def plot_maps():
    #C P 0.5 - 0.08
    z = [0.99976,0.99926,0.99865,0.99812,0.9922,0.954,0.961,0.9717,0.977,0.9829,0.809,0.847,0.879,0.9168,0.9423,0.871,0.913,0.9385,0.9519,0.9642,0.487,0.613,0.6985,0.7719,0.8212]
    plot_map(z,'confidence 0.5 - 0.08','Reds_r')
    z = [1,1,1,1,0.9968,0.9697,0.9733,0.9817,0.9874,0.9897,0.914,0.905,0.933,0.9484,0.9616,0.9683,0.9685,0.978,0.9789,0.9839,0.809,0.817,0.8722,0.899,0.9261]
    plot_map(z,'confidence 2 - 0.5','Oranges_r')
    z = [1,1,1,1,0.99961,1,1,0.9994,0.99974,0.99989,0.992,0.9951,0.9975,0.9981,0.9982,1,0.9987,0.9985,0.9991,0.9993,1,0.9976,1,0.9996,0.99968]
    plot_map(z,'confidence 100 - 2','Blues_r')

def plot_clustering_map():
    z=np.array([0,0.888889,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0,0,0,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0,0,0,0,0,0,0,0,0,0.893333,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0,0,0,0,0,0,0,0,0.893333,0.878981,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0.854031,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0,0,0,0,0,0,0.846715,0.830313,0.824,0.81085,0.797521,0.780612,0.770964,0.745724,0.745209,0.720807,0.682171,0.620837,0.55102,0.519401,0.478501,0,0,0,0])
    z = np.ma.masked_where(z < 0.05, z)
    n = np.arange(50,550,25)
    x = [a+0.5 for a in range(len(n))]
    eps = np.round(np.arange(0.000006,0.0000299,0.000001)*10e4,4)
    y = [a+0.5 for a in range(len(eps))]
    z = np.reshape(z, (20, 24))
    #plt.xscale('log')
    #plt.ylim(0,60)
    plt.clf()

    #plt.gca().ticklabel_format(axis='both', style='plain', useOffset=False)
    cmap = mpl.cm.get_cmap("autumn_r").copy()
    cmap.set_bad(color='black')

    plt.xlabel('eps * 10^4', fontsize=16)
    plt.ylabel('min pts', fontsize=16)
    plt.title('DBSCAN parameter space')
    plt.pcolor(z, cmap=cmap)
    plt.xticks(y, eps, fontsize=16)
    plt.yticks(x, n, fontsize=16)
    #plt.contourf(x,y,z)
    plt.colorbar()
    fig = plt.gcf()
    fig.set_size_inches(9, 9)
    #fig.savefig(config.output_base_path+'\\Clustering\\'+title+'.png', dpi=100)
    plt.show();

def main():
    db = Database()
    #plot_maps()
    #print(len(observed_points))
    #print(len(simulated_points))

    #op_arr = np.vstack(observed_points[:]).astype(float)
    #sp_arr = np.vstack(simulated_points[:]).astype(float)
    #plot_points(False,True)
    #plot_clustering_map()
    #plot_cluster(True,True)
    #plot_magnitude_hist()
    plot_points_velocity(False,True)
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