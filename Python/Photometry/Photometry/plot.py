import os # for relative paths
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
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

def plot_3d_fs():
    from matplotlib.colors import LogNorm
    from scipy.interpolate import interpn
    from matplotlib.colors import Normalize 
    config.database_path = os.path.join(config.output_base_path,r"Database\Cone_at_GC.db")
    db = Database()
    test = db.select_3d_stars(0, False)

    x = test[:,0]
    y = test[:,1]
    z = test[:,2]
    sort = True
    bins = 20

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    fig.set_size_inches(9.5, 5)

    ax1.scatter(x, y, z, c='r', s=0.1)
    ax1.set_zlabel('z$_{GCA}$ [pc]', rotation = 0)
    ax1.set_ylabel('y$_{GCA}$ [pc]', rotation = 0)
    ax1.set_xlabel('x$_{GCA}$ [pc]', rotation = 0)
    #plt.setp(ax1, xticks=[298,299,300,301,302],yticks=[-2,-1,0,1,2],zticks=[25,26,27,28,29])
    #ax1.set_xlim(298,302)
    #ax1.set_ylim(-2,2)
    #ax1.set_zlim(25,29)

    ax2 = fig.add_subplot(1, 2, 2)
    #ax2.set_xlim(299,301)
    #ax2.set_ylim(-1,1)

    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax2.scatter( x, y, c=z)
    ax2.set_ylabel('y$_{GCA}$ [pc]')
    ax2.set_xlabel('x$_{GCA}$ [pc]', rotation = 0)
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    #cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax2)
    #cbar.ax.set_ylabel('Density')


    #ax.title.set_text('Bulge')
    #pcm = ax.imshow(z+10, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #        cmap=cm.hot, norm=LogNorm())
    #divider = make_axes_locatable(ax2)
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    #cbar = fig.colorbar(pcm,cax=cax,ax=ax2)
    #cbar.ax.set_title('[$M_{\odot}$]')

    #ax.grid(b = True, color ='grey',
    #    linestyle ='-.', linewidth = 0.3,
    #    alpha = 0.2)

    #ax.scatter(test[:,0], test[:,1],test[:,2], s=1, c='r', marker="s", label='observed_cluster')
    fig.tight_layout(pad=3.0)
    plt.show()

def plot_3d_cluster():
    from matplotlib.colors import LogNorm
    from scipy.interpolate import interpn
    from matplotlib.colors import Normalize 
    config.database_path = os.path.join(config.output_base_path,r"Database\Default.db")#os.path.join(output_base_path,r"Database\Default_0_10000_ext.db")
    db = Database()
    test = db.select_3d_stars(0)

    x = test[:,0]
    y = test[:,1]
    z = test[:,2]
    sort = True
    bins = 20

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    fig.set_size_inches(9.5, 5)

    ax1.scatter(x, y, z, c='r', s=0.1)
    ax1.set_zlabel('z$_{GCA}$ [pc]', rotation = 0)
    ax1.set_ylabel('y$_{GCA}$ [pc]', rotation = 0)
    ax1.set_xlabel('x$_{GCA}$ [pc]', rotation = 0)
    ax1.set_xlim(298,302)
    plt.setp(ax1, xticks=[298,299,300,301,302],yticks=[-2,-1,0,1,2],zticks=[25,26,27,28,29])
    #ax1.set_xticks(298,299,300,301,302)
    ax1.set_ylim(-2,2)
    ax1.set_zlim(25,29)

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_xlim(299,301)
    ax2.set_ylim(-1,1)
    ax2.set_xticks([299,299.5,300,300.5,301])

    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax2.scatter( x, y, c=z)
    ax2.set_ylabel('y$_{GCA}$ [pc]')
    ax2.set_xlabel('x$_{GCA}$ [pc]', rotation = 0)
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    #cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax2)
    #cbar.ax.set_ylabel('Density')


    #ax.title.set_text('Bulge')
    #pcm = ax.imshow(z+10, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #        cmap=cm.hot, norm=LogNorm())
    #divider = make_axes_locatable(ax2)
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    #cbar = fig.colorbar(pcm,cax=cax,ax=ax2)
    #cbar.ax.set_title('[$M_{\odot}$]')

    #ax.grid(b = True, color ='grey',
    #    linestyle ='-.', linewidth = 0.3,
    #    alpha = 0.2)

    #ax.scatter(test[:,0], test[:,1],test[:,2], s=1, c='r', marker="s", label='observed_cluster')
    fig.tight_layout(pad=3.0)
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

    plt.rcParams.update({'font.size': 14})

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

    plt.setp(axes, xticks=x, xticklabels=masses*0.001, yticks=y, yticklabels=angles)

    c = axes[0].pcolor(np.reshape(z[0], (5, 5)), cmap='Reds_r')
    fig.colorbar(c, ax=axes[0], pad=0.01)
    axes[0].set_title('0.5 - 0.08 [$M_{\odot}$]')
    axes[0].set_ylabel("angle [$^\circ$]", fontsize=14)

    c = axes[1].pcolor(np.reshape(z[1], (5, 5)), cmap='Oranges_r')
    fig.colorbar(c, ax=axes[1], pad=0.01)
    axes[1].set_title('2 - 0.5 [$M_{\odot}$]')
    axes[1].set_xlabel("cluster mass [$kM_{\odot}$]", fontsize=14)

    c = axes[2].pcolor(np.reshape(z[2], (5, 5)), cmap='Blues_r')
    fig.colorbar(c, ax=axes[2], pad=0.01)
    axes[2].set_title('100 - 2 [$M_{\odot}$]')

    fig.tight_layout()

    plt.show()


def plot_velocity_hist():

    config.simulation_id=2
    db = Database()

    fig, axes = plt.subplots(1, 2, figsize=(10,5))

    scalefactor = 28

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
    axes[1].set_title('10 [kM$_{\odot}$]')
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axes[1].ticklabel_format(axis="x", style="sci", scilimits=(0,0))

    #plt.hist(cluster_arr[:,i], density=True, bins=100, alpha=0.5, label='CS')
    #plt.hist(fs_arr[:,i], density=True, bins=100, alpha=0.5, label='FS')
    axes[0].set_xlabel('v_asc [arcsec/day]', fontsize=16)
    axes[0].set_ylabel('v_dec [arcsec/day]', fontsize=16)
    axes[0].set_title('0.64 [kM$_{\odot}$]')
    axes[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axes[0].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
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
    mass_range = 'Tot'

    query = 'Angle==10'
    SNCSTotdf = df.query(query)['SNCS '+mass_range]
    SNCSTotdferr = df.query(query)['SNCS '+mass_range+' Err']
    SNFSotdf = df.query(query)['SNFS '+mass_range]
    SNFSTotdferr = df.query(query)['SNFS '+mass_range+' Err']
    MNCSTotdf = df.query(query)['MNCS '+mass_range]
    MNCSTotdferr = df.query(query)['MNCS '+mass_range+' Err']
    MNFSTotdf = df.query(query)['MNFS '+mass_range]
    MNFSTotdferr = df.query(query)['MNFS '+mass_range+' Err']

    print(SNCSTotdf)
    print(MNCSTotdf)

    N = len(masses)

    ## necessary variables
    ind = np.arange(N)                # the x locations for the groups
    width = 0.35                      # the width of the bars

    ## the bars
    rects1 = ax.bar(ind, SNCSTotdf, width, alpha=0.5,
                    color='blue',
                    yerr=SNCSTotdferr,
                    error_kw=dict(elinewidth=2,ecolor='red'))

    rects2 = ax.bar(ind, MNCSTotdf, width, alpha=1,
                    color='blue',
                    yerr=MNCSTotdferr,
                    error_kw=dict(elinewidth=2,ecolor='red'))

    rects3 = ax.bar(ind+width, SNFSotdf, width, alpha=0.5,
                    color='black',
                    yerr=SNFSTotdferr,
                    error_kw=dict(elinewidth=2,ecolor='red'))

    rects4 = ax.bar(ind+width, MNFSTotdf, width, alpha=1,
                    color='black',
                    yerr=MNFSTotdferr,
                    error_kw=dict(elinewidth=2,ecolor='red'))



    # axes and labels
    ax.set_xlim(-width,len(ind)+width*2)
    ax.set_xlabel('cluster mass [kM$_{\odot}$]')
    ax.set_ylabel('#Stars')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax.set_title('Amount of stars by clustersize')
    xTickMarks = masses/1000
    ax.set_xticks(ind+width/2)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=0, fontsize=10)

    ## add a legend
    ax.legend( (rects1[0], rects2[0], rects3[0], rects4[0]), ('SCS','MCS', 'SFS','MFS') )

    plt.show()

def plot_avg_2D_vel():
    #640

    scale_factor = 28 #dt in days

    disp_2D_cluster = np.array([0.00310743,0.00316207,0.0029074,0.0026192,0.0021763]) * scale_factor
    disp_2D_cluster_error = np.array([0.00000084,0.00000076,0.0000027,0.0000066,0.0000090]) * scale_factor
    avg_vel_2D_cluster = np.array([0.00310743,0.00316207,0.0029074,0.0026192,0.0021763]) * scale_factor
    avg_vel_2D_cluster_error = np.array([0.00000084,0.00000076,0.0000027,0.0000066,0.0000090]) * scale_factor
    avg_vel_2D_fs = np.array([0.00826,0.002124,0.003193,0.005075,0.0036446]) * scale_factor
    avg_vel_2D_fs_error = np.array([0.00022,0.000015,0.000014,0.000014,0.0000056]) * scale_factor



    x = [180,25,10,5,0]
    a = np.arange(min(x),max(x),np.abs((x[0]-x[-1])/len(x)))
    plt.xticks(a,x)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    plt.scatter(a,avg_vel_2D_cluster,label='cluster stars')

    plt.scatter(a,avg_vel_2D_fs,marker='x',label='field stars')

    # Plot error bar

    plt.errorbar(a, avg_vel_2D_cluster, yerr = avg_vel_2D_cluster_error,fmt='o',ecolor = 'cyan')

    plt.xlabel('angle [$^\circ$]', fontsize=16)

    plt.ylabel('velocity [arcsec/day]', fontsize=16)

    plt.legend()

    plt.show()


def plot_title_img():
    # Create main container
    fig = plt.figure()
    ax = plt.subplot(111) #whole path

    config.database_path = os.path.join(config.output_base_path,r"Database\Default_0_10000_ext.db")
    db = Database()

    simulated_points = db.select_points(0, False)
    sp_arr = np.vstack(simulated_points[:]).astype(float)
    sp_cluster = sp_arr[sp_arr[:,4] > -1]
    #plt.scatter(sp_cluster[:,2], sp_cluster[:,3], s=1, c='g', marker="s", label='simulated_cluster')
    sp_fs = sp_arr[sp_arr[:,4] == -1]
    #plt.scatter(sp_fs[:,2], sp_fs[:,3], s=1, c='blue', marker="s", label='simulated_fs')

    #plt.legend(loc='upper left')
    ax.scatter(sp_cluster[:,2], sp_cluster[:,3], s = 5, c = 'g')
    ax.scatter(sp_fs[:,2], sp_fs[:,3], s = 5, c = 'blue')

    # Create zoom-out plot
    rect_x = 0.002145
    rect_y = 0.0001765
    rect_width = 0.00049
    rect_height = 0.000587
    #rect = mpl.patches.Rectangle((rect_x, rect_y), rect_width, rect_height, linewidth=1, edgecolor='r', facecolor='none')
    # Add the patch to the Axes
    #plt.gca().add_patch(rect)
    plt.xlim(-0.019, 0.019)
    plt.ylim(-0.017, 0.017)
    plt.xlabel('v$_{x,\mathrm{HTP}}$ [arcsec/dt]', fontsize=16)
    plt.ylabel('v$_{y,\mathrm{HTP}}$ [arcsec/dt]', fontsize=16)
    plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
    plt.ticklabel_format(axis='y',style='sci', scilimits=(0,0))

    # Create zoom-in plot
    axins = zoomed_inset_axes(ax,25,loc='lower right', 
                              axes_kwargs={"facecolor" : "lightgray"})

    #axins.plot(random_walk)

    #x1,x2,y1,y2 = 1000,2000, -60,-15
    axins.set_xlim(rect_x-rect_width/2, rect_x+rect_width/2)
    axins.set_ylim(rect_y-rect_height/2, rect_y+rect_height/2)

    #ax_new = fig.add_axes([0.6, 0.6, 0.4, 0.4]) # the position of zoom-out plot compare to the ratio of zoom-in plot 
    #axins.axis('off')
    axins.get_xaxis().set_visible(False)
    axins.get_yaxis().set_visible(False)
    #ax_new.spines['bottom'].set_visible(True)
    #ax_new.spines['left'].set_visible(True)
    #plt.xlim(rect_x-rect_width/2, rect_x+rect_width/2)
   # plt.ylim(rect_y-rect_height/2, rect_y+rect_height/2)
    axins.scatter(sp_cluster[:,2], sp_cluster[:,3], s = 1, c = 'g')
    axins.scatter(sp_fs[:,2], sp_fs[:,3], s = 1, c = 'blue')


    pp,p1,p2 = mark_inset(ax,axins,loc1=1,loc2=3)
    pp.set_fill(True)
    pp.set_facecolor("lightgray")
    pp.set_edgecolor("k")


    # Save figure with nice margin
    #plt.savefig('zoom.png', dpi = 300, bbox_inches = 'tight', pad_inches = .1)
    plt.show()

def plot_f1_maps():

    plt.rcParams.update({'font.size': 14})

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(12,5))

    df = pd.read_excel(config.output_base_path+r'\25_observations.xlsx',
                   'extinction', usecols = 'A:AG')

    masses = df['Mass'].unique()
    x = [a+0.5 for a in range(len(masses))]
    angles = df['Angle'].unique()
    y = [a+0.5 for a in range(len(angles))]

    plt.setp(axes, xticks=x, xticklabels=masses*0.001, yticks=y, yticklabels=angles)

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
    axes[1].set_xlabel("cluster mass [$kM_{\odot}$]", fontsize=14)

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

def plot_fits(plot_observed = False, plot_simulated = False, gc=True):

    from astropy.visualization import astropy_mpl_style
    plt.style.use(astropy_mpl_style)
    from astropy.utils.data import get_pkg_data_filename
    from astropy.io import fits
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle
    if gc:
        image_file = get_pkg_data_filename(r'scopesim_t0.fits')
    else:
        image_file = get_pkg_data_filename(r'scopesim_t0_ac.fits')#scopesim_t0_ac   scopesim_t0.fits
    image_data = fits.getdata(image_file, ext=0)

    from numpy import unravel_index
    max_index = unravel_index(image_data.argmax(), image_data.shape)
    max_value = image_data[max_index[0],max_index[1]]
    min_index = unravel_index(image_data.argmin(), image_data.shape)
    min_value = image_data[min_index[0],min_index[1]]
    print(max_index)

    print(image_data)

    print(max_value)
    print(min_value)
    #return

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6), tight_layout = True)
    ax1.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    ax2.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    ax3.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

    cmap = 'hot'
    norm = LogNorm(max_value/20000,max_value,clip=True)

    if gc:
        x_min_1 = 6800
        y_min_1 = 6800
        dxy_1 = 700
        x_min_2 = 6826
        y_min_2 = 7203
        dxy_2 = 233
    else:
        x_min_1 = 6400
        y_min_1 = 7000
        dxy_1 = 700
        x_min_2 = 6700
        y_min_2 = 7203
        dxy_2 = 233

    #ax1 = fig.add_subplot(1, 3, 1)
    ax1.set_title('Full Image')
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax1.imshow(image_data, norm=norm, cmap=cmap, origin='lower')
    #ax1.imshow(image_data, norm=LogNorm(), cmap=cmap, origin='lower')
    ax1.add_patch(Rectangle((x_min_1,y_min_1),dxy_1,dxy_1, facecolor="none", ec='w', lw=2))

    #ax2 = fig.add_subplot(1, 3, 2)
    ax2.set_title('21x Zoom')
    ax2.set_xlim(x_min_1,x_min_1+dxy_1)
    ax2.set_ylim(y_min_1,y_min_1+dxy_1)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax2.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax2.imshow(image_data, norm=norm, cmap=cmap, origin='lower')
    ax2.add_patch(Rectangle((x_min_2,y_min_2),dxy_2,dxy_2, facecolor="none", ec='w', lw=2))

    #ax3 = fig.add_subplot(1, 3, 3)
    ax3.set_title('64x Zoom')
    ax3.set_xlim(x_min_2, x_min_2 + dxy_2)
    ax3.set_ylim(y_min_2, y_min_2 + dxy_2)
    im = ax3.imshow(image_data, norm=norm, cmap=cmap , origin='lower')

    if plot_simulated or plot_observed:
        config.database_path = os.path.join(config.output_base_path,r"Database\Default_0_10000_ext.db")
        db = Database()

    if plot_simulated:
        simulated_points = db.select_points(0, False)
        sim_arr = np.vstack(simulated_points[:]).astype(float)
        origin = 14976/2.
        x = sim_arr[:,0]/config.pixelfactor + origin
        y = sim_arr[:,1]/config.pixelfactor + origin
        ax1.scatter(x, y, s=0.1, c='cyan', marker="x", label='simulated')
        ax2.scatter(x, y, s=10, c='cyan', marker="x", label='simulated')
        ax3.scatter(x, y, s=50, c='cyan', marker="x", label='simulated')

    if plot_observed:
        observed_points = db.select_points(0, True)
        op_arr = np.vstack(observed_points[:]).astype(float)
        origin = 14976/2.
        x = op_arr[:,0]/config.pixelfactor + origin
        y = op_arr[:,1]/config.pixelfactor + origin
        ax1.scatter(x, y, s=0.1, facecolors='none', edgecolors='chartreuse', marker="o", label='observed')
        ax2.scatter(x, y, s=10, facecolors='none', edgecolors='chartreuse', marker="o", label='observed')
        ax3.scatter(x, y, s=50, facecolors='none', edgecolors='chartreuse', marker="o", label='observed')

    if plot_simulated and plot_observed:
        lgnd = ax1.legend(loc="upper right", scatterpoints=1, fontsize=10)
        for lh in lgnd.legendHandles: 
            lh._sizes = [30]
            lh.set_alpha(1)

    #fig.supxlabel()
    #ax2.set_xlabel('x$_{HTP}$ [px]')
    #fig.supylabel('y$_{HTP}$ [px]')
    fig.tight_layout()
    plt.show()

def main():

    #plot_3d_cluster()
    plot_3d_fs()
    #plot_fits(False,False,False)
    #plot_fits(True,True)

    #plot_title_img()

    #plot_number_hist()
    #plot_avg_2D_vel()
    #plot_velocity_hist()

    #plot_title_img()

    #db = Database()
    #plot_precision_maps(False)
    #plot_f1_maps()

    #ErrorAnalysis.Vel2D_error()
    #ErrorAnalysis.Precision_error('0.5 - 0.08',False)

    #Vel2D_error()
    

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