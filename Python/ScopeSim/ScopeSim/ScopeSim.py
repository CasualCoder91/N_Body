import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import scopesim as sim
import pyckles
import scopesim_templates as sim_tp
#import astropy.table as table
from astropy.table import Table
from astropy import units
from matplotlib.colors import LogNorm
import os # for relative paths

def createConnection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except:
        print("Connection failed!")

    return conn

def select2dPositions(conn, simulationID, timestep):
    cur = conn.cursor()
    cur.execute("""SELECT position.rH, position.aHTP, position.dHTP, star.mass 
       FROM star
       INNER JOIN position on position.id_star = star.id
       where star.id_simulation = ?1
       and position.timestep = ?2""", (simulationID,timestep)) #and star.isCluster=1 LIMIT 10000000 OFFSET 10000000
    rows = np.array(cur.fetchall())
    return rows

#and star.isCluster = 1

def makeSource(data):
    """
    Generate a source object from np.array (radius,ascension,declination,mass)

    Returns
    -------
    src: scopesim.Source
    """

    masses = data[:,3]
    distances = data[:,0]

    # 2. get spec_types for masses
    spec_types = sim_tp.utils.cluster_utils.mass2spt(masses)
    spec_types = sim_tp.utils.cluster_utils.closest_pickles(spec_types)

    # 3. get spectra from pyckles
    pickles = pyckles.SpectralLibrary("pickles", return_style="synphot")
    unique_spts = np.unique(spec_types)
    spectra = [pickles[spt] for spt in unique_spts]

    # 4. scale all spectra to V=0
    spectra = [sim_tp.rc.ter_curve_utils.scale_spectrum(spec, "V", 0*units.mag) for spec in spectra]

    # 5. make ref list for spec_types from spectra
    ref = [np.where([spt == u_spt for u_spt in unique_spts])[0][0]
           for spt in spec_types]

    # 6. make weight list from Mv + dist_mod(distance)
    Mvs = np.array(sim_tp.utils.cluster_utils.mass2Mv(masses))
    dist_mod = 5 * np.log10(distances) - 5
    weight = 10 ** (-0.4 * (Mvs + dist_mod))


    # 8. make table with (x,y,ref,weight)
    tbl = Table(names=["x", "y", "ref", "weight", "masses", "spec_types"],
                data= [ data[:,1],   data[:,2],   ref,   weight,   masses,   spec_types ])

    # 9. make Source with table, spectra
    src = sim_tp.rc.Source(table=tbl, spectra=spectra)

    return src

def main():

    simulationID = 1
    bPlot = False #If True then Plot result

    #paths
    outputBasePath = os.path.join(os.path.abspath(__file__ + r"\..\..\..\.."), r"Output")
    databasePath = os.path.join(outputBasePath,r"Database\Default.db")
    outputPath = os.path.join(outputBasePath, "Simulation" + str(simulationID))

    conn = createConnection(databasePath) # create a database connection
    #sim.download_package(["instruments/MICADO_Sci",
    #                      "telescopes/ELT", 
    #                      "locations/Armazones", 
    #                      "instruments/MICADO"]) #"MAORY"

    with conn:

        cmd = sim.UserCommands(use_instrument="MICADO_Sci")
        cmd["!DET.width"] = 4096
        cmd["!DET.height"] = 4096

        opt = sim.OpticalTrain(cmd)

        data = select2dPositions(conn, simulationID, 0)
        source = makeSource(data)
        opt.observe(source)
        opt.readout(filename=outputPath+r"\scopesim.fits")

        if bPlot:
            plt.figure(figsize=(12,12))
            plt.imshow(opt.image_planes[0].image, norm=LogNorm())
            plt.show()


if __name__ == '__main__':
    main()

