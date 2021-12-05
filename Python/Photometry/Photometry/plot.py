import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import pandas as pd

from database import Database
from point import Point
import ErrorAnalysis
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

def plot_precision_maps(simulated=False):

    if simulated:
        z = [[1,1,1,1,0.998,0.9943,0.9931,0.9947,0.99523,0.99591,0.9689,0.9749,0.9751,0.9801,0.98332,0.9878,0.9909,0.9912,0.99224,0.9931,0.893,0.916,0.9239,0.9405,0.9475],
             [1,1,1,1,0.99847,0.9946,0.9931,0.995,0.99663,0.99654,0.978,0.9742,0.9787,0.9841,0.9861,0.9925,0.9918,0.9926,0.9943,0.99459,0.933,0.9399,0.9508,0.9636,0.968],
             [1,1,1,1,0.99981,1,1,0.9996,0.99984,0.99994,1,0.9981,0.9992,0.9989,0.99937,1,1,1,0.99984,0.99975,1,1,1,0.99984,0.99969]]
    else:
        z = [[0.99976,0.99926,0.99865,0.99812,0.9922,0.954,0.961,0.9717,0.977,0.9829,0.809,0.847,0.879,0.9168,0.9423,0.871,0.913,0.9385,0.9519,0.9642,0.487,0.613,0.6985,0.7719,0.8212],
             [1,1,1,1,0.9968,0.9697,0.9733,0.9817,0.9874,0.9897,0.914,0.905,0.933,0.9484,0.9616,0.9683,0.9685,0.978,0.9789,0.9839,0.809,0.817,0.8722,0.899,0.9261],
             [1,1,1,1,0.99961,1,1,0.9994,0.99974,0.99989,0.992,0.9951,0.9975,0.9981,0.9982,1,0.9987,0.9985,0.9991,0.9993,1,0.9976,1,0.9996,0.99968]]

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(12,5))

    masses =np.array([640,1600,4000,10000,25000])
    x = [a+0.5 for a in range(len(masses))]
    angles = np.array([180,25,10,5,0])
    y = [a+0.5 for a in range(len(angles))]

    plt.setp(axes, xticks=x, xticklabels=masses, yticks=y, yticklabels=angles)

    c = axes[0].pcolor(np.reshape(z[0], (5, 5)), cmap='Reds_r')
    fig.colorbar(c, ax=axes[0], pad=0.01)
    axes[0].set_title('0.5 - 0.08 [$M_{\odot}$]')
    axes[0].set_ylabel("angle [$^\circ$]", fontsize=14)

    c = axes[1].pcolor(np.reshape(z[1], (5, 5)), cmap='Oranges_r')
    fig.colorbar(c, ax=axes[1], pad=0.01)
    axes[1].set_title('2 - 0.5 [$M_{\odot}$]')
    axes[1].set_xlabel("cluster mass [$M_{\odot}$]", fontsize=14)

    c = axes[2].pcolor(np.reshape(z[2], (5, 5)), cmap='Blues_r')
    fig.colorbar(c, ax=axes[2], pad=0.01)
    axes[2].set_title('100 - 2 [$M_{\odot}$]')

    fig.tight_layout()

    plt.show()


def plot_velocity_hist():

    config.simulation_id=2
    db = Database()

    fig, axes = plt.subplots(1, 2, figsize=(10,5))

    scalefactor = 0.030866666667*2.419e+6

    #data
    db.connect(config.output_base_path + r"\Database\Default_0_640_ext.db")
    cluster = db.select_cluster(0,False)
    cluster_arr = np.vstack(cluster[:]).astype(float)[:,2:4]*scalefactor
    cluster_arr = ErrorAnalysis.remove_outliers(cluster_arr)
    fs = db.select_fs(0,False)
    fs_arr = np.vstack(fs[:]).astype(float)[:,2:4]*scalefactor

    axes[0].scatter(fs_arr[:,0],fs_arr[:,1], marker='o', alpha=0.5, label='FS')
    axes[0].scatter(cluster_arr[:,0],cluster_arr[:,1], marker='x', alpha=0.6, label='CS')


 #   axes[0].hist2d(cluster_arr[:, 0], 
 #          cluster_arr[:, 1],
 #          bins = 100, 
 #          cmap = "RdYlGn_r",
 #          norm = mpl.colors.LogNorm())

    cluster_mean = np.mean(cluster_arr, axis=0)
    cluster_std = np.std(cluster_arr, axis=0)
    dx = np.max(cluster_std)*4

    axes[0].set_xlim(cluster_mean[0]-dx,cluster_mean[0]+dx)
    axes[0].set_ylim(cluster_mean[1]-dx,cluster_mean[1]+dx)
    print(axes[0].get_xlim())
    print(axes[1].get_xlim())

    axes[1].set_xlim(axes[0].get_xlim())
    axes[1].set_ylim(axes[0].get_ylim())

    from matplotlib.patches import Circle
    radius = np.max(cluster_std)*2
    circle = Circle(cluster_mean, radius=radius, color='r', fill=False)
    norm = np.linalg.norm(cluster_arr-cluster_mean,axis=1)
    print('points inside circle: ', len(norm[norm < radius]))


    circle2 = Circle(cluster_mean, radius=radius, color='r', fill=False)

    axes[0].add_patch(circle)
    axes[1].add_patch(circle2)

    db.connect(config.output_base_path + r"\Database\Default_0_10000_ext.db")

    cluster = db.select_cluster(0, False)
    cluster_arr = np.vstack(cluster[:]).astype(float)[:,2:4]*scalefactor
    norm = np.linalg.norm(cluster_arr-cluster_mean,axis=1)
    print('points inside circle: ', len(norm[norm < radius]))


    fs = db.select_fs(0, False)
    fs_arr = np.vstack(fs[:]).astype(float)[:,2:4]*scalefactor

    axes[1].scatter(fs_arr[:,0],fs_arr[:,1], marker='o', alpha=0.5, label='FS')
    axes[1].scatter(cluster_arr[:,0],cluster_arr[:,1], marker='x', alpha=0.6, label='CS')

 #   axes[1].hist2d(cluster_arr[:, 0], 
 #          cluster_arr[:, 1],
 #          bins = 5000, 
 #          cmap = "RdYlGn_r",
 #          norm = mpl.colors.LogNorm(),
 #          alpha = 0.5)

    axes[1].set_xlim(axes[0].get_xlim())
    axes[1].set_ylim(axes[0].get_ylim())
    axes[1].set_title('10 kM[$_{\odot}$]')


    #plt.hist(cluster_arr[:,i], density=True, bins=100, alpha=0.5, label='CS')
    #plt.hist(fs_arr[:,i], density=True, bins=100, alpha=0.5, label='FS')
    axes[0].set_xlabel('v_asc [km/s]', fontsize=16)
    axes[0].set_ylabel('v_dec [km/s]', fontsize=16)
    axes[0].set_title('0.64 kM[$_{\odot}$]')
    axes[0].legend(loc="upper left")

    plt.show()

def plot_number_hist():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ## the data
    df = pd.read_excel(config.output_base_path+r'\25_observations.xlsx','PlotData', usecols = 'A:BJ')
    masses = df['Mass'].unique()
    columns = ['SNCS Tot','SNFS Tot','MNCS Tot','MNFS Tot']

    #mass_ranges = ['Tot','> 2','2 - 0.5','0.5 - 0.08']
    mass_range = '0.5 - 0.08'

    query = 'Angle==10'
    SNCSTotdf = df.query(query)['SNCS '+mass_range]
    SNCSTotdferr = df.query(query)['SNCS '+mass_range+' Err']
    SNFSotdf = df.query(query)['SNFS '+mass_range]
    SNFSTotdferr = df.query(query)['SNFS '+mass_range+' Err']
    MNCSTotdf = df.query(query)['MNCS '+mass_range]
    MNCSTotdferr = df.query(query)['MNCS '+mass_range+' Err']
    MNFSTotdf = df.query(query)['MNFS '+mass_range]
    MNFSTotdferr = df.query(query)['MNFS '+mass_range+' Err']

    N = len(masses)

    ## necessary variables
    ind = np.arange(N)                # the x locations for the groups
    width = 0.35                      # the width of the bars

    ## the bars
    rects1 = ax.bar(ind, SNCSTotdf, width, alpha=0.5,
                    color='black',
                    yerr=SNCSTotdferr,
                    error_kw=dict(elinewidth=2,ecolor='red'))

    rects2 = ax.bar(ind+width, SNFSotdf, width, alpha=0.5,
                        color='red',
                        yerr=SNFSTotdferr,
                        error_kw=dict(elinewidth=2,ecolor='black'))

    rects3 = ax.bar(ind, MNCSTotdf, width, alpha=0.5,
                    color='black',
                    yerr=MNCSTotdferr,
                    error_kw=dict(elinewidth=2,ecolor='red'))

    rects4 = ax.bar(ind+width, MNFSTotdf, width, alpha=0.5,
                        color='red',
                        yerr=MNFSTotdferr,
                        error_kw=dict(elinewidth=2,ecolor='black'))



    # axes and labels
    ax.set_xlim(-width,len(ind)+width)
    #ax.set_ylim(0,45)
    ax.set_ylabel('#Stars')
    #ax.set_title('Amount of stars by clustersize')
    xTickMarks = masses
    ax.set_xticks(ind+width/2)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)

    ## add a legend
    ax.legend( (rects1[0], rects2[0]), ('CS', 'FS') )

    plt.show()

def plot_avg_2D_vel():
    #640
    disp_2D_cluster = [0.00310743,0.00316207,0.0029074,0.0026192,0.0021763]
    disp_2D_cluster_error = [0.00000084,0.00000076,0.0000027,0.0000066,0.0000090]
    avg_vel_2D_cluster = [0.00310743,0.00316207,0.0029074,0.0026192,0.0021763]
    avg_vel_2D_cluster_error = [0.00000084,0.00000076,0.0000027,0.0000066,0.0000090]
    avg_vel_2D_fs = [0.00826,0.002124,0.003193,0.005075,0.0036446]
    avg_vel_2D_fs_error = [0.00022,0.000015,0.000014,0.000014,0.0000056]

    x = [180,25,10,5,0]
    a = np.arange(min(x),max(x),np.abs((x[0]-x[-1])/len(x)))
    plt.xticks(a,x)

    plt.scatter(a,avg_vel_2D_cluster,label='cluster stars')

    plt.scatter(a,avg_vel_2D_fs,marker='x',label='field stars')

    # Plot error bar

    plt.errorbar(a, avg_vel_2D_cluster, yerr = avg_vel_2D_cluster_error,fmt='o',ecolor = 'cyan')

    plt.xlabel('angle [$^\circ$]', fontsize=16)
    plt.ylabel('velocity [arcsec/dt]', fontsize=16)

    plt.legend()

    plt.show()


def plot_f1_maps():

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(12,5))

    df = pd.read_excel(config.output_base_path+r'\25_observations.xlsx',
                   'extinction', usecols = 'A:AG')

    masses = df['Mass'].unique()
    x = [a+0.5 for a in range(len(masses))]
    angles = df['Angle'].unique()
    y = [a+0.5 for a in range(len(angles))]

    plt.setp(axes, xticks=x, xticklabels=masses, yticks=y, yticklabels=angles)

    z = [0.635,0.8125,0.8476,0.925,0.9403,0.626,0.747,0.8098,0.8719,0.9163,0.45,0.621,0.7224,0.7979,0.8667,0.338,0.486,0.6061,0.7196,0.8117,0.39,0.518,0.619,0.7053,0.7637]
    z = np.reshape(z, (5, 5))
    c = axes[0].pcolor(z, cmap='Reds_r')
    fig.colorbar(c, ax=axes[0], pad=0.01)
    axes[0].set_title('0.5 - 0.08 [$M_{\odot}$]')
    axes[0].set_ylabel("angle [$^\circ$]", fontsize=14)

    z = [0.732,0.9068,0.924,0.9531,0.9582,0.811,0.893,0.9355,0.9476,0.9547,0.751,0.875,0.9079,0.9318,0.9467,0.782,0.865,0.9233,0.9410,0.9539,0.786,0.828,0.8837,0.8978,0.9174]
    z = np.reshape(z, (5, 5))
    c = axes[1].pcolor(z, cmap='Oranges_r')
    fig.colorbar(c, ax=axes[1], pad=0.01)
    axes[1].set_title('2 - 0.5 [$M_{\odot}$]')
    axes[1].set_xlabel("cluster mass [$M_{\odot}$]", fontsize=14)

    z = [0.712,0.944,0.9075,0.9340,0.9492,0.897,0.924,0.9278,0.9511,0.9556,0.744,0.9384,0.9507,0.9594,0.9671,0.838,0.907,0.9448,0.9585,0.9705,0.926,0.9622,0.9536,0.9641,0.9657]
    z = np.reshape(z, (5, 5))
    c = axes[2].pcolor(z, cmap='Blues_r')
    fig.colorbar(c, ax=axes[2], pad=0.01)
    axes[2].set_title('100 - 2 [$M_{\odot}$]')

    fig.tight_layout()

    plt.show()

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

    plt.xlabel(r'$\epsilon * 10^5$', fontsize=16)
    plt.ylabel('nPoints', fontsize=16)
    plt.title('DBSCAN parameter space')
    plt.pcolor(z, cmap=cmap)
    plt.xticks(y, eps, fontsize=16)
    plt.yticks(x, n, fontsize=16)
    #plt.contourf(x,y,z)
    plt.colorbar()
    fig = plt.gcf()
    #fig.set_size_inches(9, 9)
    #fig.savefig(config.output_base_path+'\\Clustering\\'+title+'.png', dpi=100)
    plt.show();

def main():

    #plot_number_hist()

    plot_velocity_hist()

    #db = Database()
    #plot_precision_maps(True)
    #plot_f1_maps()

    #ErrorAnalysis.Vel2D_error()
    #ErrorAnalysis.Precision_error('0.5 - 0.08',False)

    #Vel2D_error()
    #plot_avg_2D_vel()

    #mass_ranges = ['Tot','> 2','2 - 0.5','0.5 - 0.08']

    #ErrorAnalysis.F1_error('0.5 - 0.08')

    #print(len(observed_points))
    #print(len(simulated_points))

    #op_arr = np.vstack(observed_points[:]).astype(float)
    #sp_arr = np.vstack(simulated_points[:]).astype(float)
    #plot_points(False,True)
    #plot_clustering_map()
    #plot_cluster(True,True)
    #plot_magnitude_hist()
    #plot_points_velocity(False,True)
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