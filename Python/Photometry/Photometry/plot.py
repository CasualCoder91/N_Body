import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import pandas as pd
from uncertainties import ufloat

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

def plot_precision_maps():

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(12,5))

    masses =np.array([640,1600,4000,10000,25000])
    x = [a+0.5 for a in range(len(masses))]
    angles = np.array([180,25,10,5,0])
    y = [a+0.5 for a in range(len(angles))]

    plt.setp(axes, xticks=x, xticklabels=masses, yticks=y, yticklabels=angles)

    z = [0.99976,0.99926,0.99865,0.99812,0.9922,0.954,0.961,0.9717,0.977,0.9829,0.809,0.847,0.879,0.9168,0.9423,0.871,0.913,0.9385,0.9519,0.9642,0.487,0.613,0.6985,0.7719,0.8212]
    z = np.reshape(z, (5, 5))
    c = axes[0].pcolor(z, cmap='Reds_r')
    fig.colorbar(c, ax=axes[0], pad=0.01)
    axes[0].set_title('0.5 - 0.08 [$M_{\odot}$]')
    axes[0].set_ylabel("angle [$^\circ$]", fontsize=14)

    z = [1,1,1,1,0.9968,0.9697,0.9733,0.9817,0.9874,0.9897,0.914,0.905,0.933,0.9484,0.9616,0.9683,0.9685,0.978,0.9789,0.9839,0.809,0.817,0.8722,0.899,0.9261]
    z = np.reshape(z, (5, 5))
    c = axes[1].pcolor(z, cmap='Oranges_r')
    fig.colorbar(c, ax=axes[1], pad=0.01)
    axes[1].set_title('2 - 0.5 [$M_{\odot}$]')
    axes[1].set_xlabel("cluster mass [$M_{\odot}$]", fontsize=14)

    z = [1,1,1,1,0.99961,1,1,0.9994,0.99974,0.99989,0.992,0.9951,0.9975,0.9981,0.9982,1,0.9987,0.9985,0.9991,0.9993,1,0.9976,1,0.9996,0.99968]
    z = np.reshape(z, (5, 5))
    c = axes[2].pcolor(z, cmap='Blues_r')
    fig.colorbar(c, ax=axes[2], pad=0.01)
    axes[2].set_title('100 - 2 [$M_{\odot}$]')

    fig.tight_layout()

    plt.show()

def RoundToSigFigs( x, sigfigs ):

    #temporarily suppres floating point errors
    errhanddict = np.geterr()
    np.seterr(all="ignore")

    matrixflag = False
    if isinstance(x, np.matrix): #Convert matrices to arrays
        matrixflag = True
        x = np.asarray(x)

    xsgn = np.sign(x)
    absx = xsgn * x
    mantissas, binaryExponents = np.frexp( absx )

    decimalExponents = 3.010299956639811952137388947244930267681898814621085413104274611e-1 * binaryExponents
    omags = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - omags)

    if type(mantissas) is float or isinstance(mantissas, np.floating):
        if mantissas < 1.0:
            mantissas *= 10.0
            omags -= 1.0
        
    else: #elif np.all(np.isreal( mantissas )):
        fixmsk = mantissas < 1.0, 
        mantissas[fixmsk] *= 10.0
        omags[fixmsk] -= 1.0

    result = xsgn * np.around( mantissas, decimals=sigfigs - 1 ) * 10.0**omags
    if matrixflag:
        result = np.matrix(result, copy=False)

    np.seterr(**errhanddict)
    return result

def F1_error(mass_range):
    df = pd.read_excel(config.output_base_path+r'\25_observations.xlsx',
                   'extinction', usecols = 'A:AG')

    masses = df['Mass'].unique()
    x = [a+0.5 for a in range(len(masses))]
    angles = df['Angle'].unique()
    y = [a+0.5 for a in range(len(angles))]

    for angle in angles:
        for mass in masses:
            UPdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['UP '+mass_range]
            UP = ufloat(UPdf.mean(), UPdf.std())
            CFPdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['CFP '+mass_range]
            CFP = ufloat(CFPdf.mean(), CFPdf.std())
            CTPdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['CTP '+mass_range]
            CTP = ufloat(CTPdf.mean(), CTPdf.std())
            CFNdf = df.query('Mass=='+str(mass)+'&'+'Angle=='+str(angle))['CFN '+mass_range]
            CFN = ufloat(CFNdf.mean(), CFNdf.std())
            F1 = CTP/(CTP+0.5*(CFP+UP+CFN))
            #print(uncertainty.s)
            digits = np.abs(int(np.log10(abs(F1.s))-2))
            F1_round = np.round(F1.n,digits)
            format_string = '{:.'+str(digits)+'f}'
            print(format_string.format(F1_round))

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
    #db = Database()
    #plot_precision_maps()
    plot_f1_maps()

    #F1_error('> 2')

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