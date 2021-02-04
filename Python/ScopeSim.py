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

def createConnection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn

def select2dPositions(conn, simulationID, timestep):
    cur = conn.cursor()
    cur.execute("""SELECT position.rHEQ,position.aScope,position.dScope,star.mass FROM star
       INNER JOIN position on position.id_star = star.id
       where star.id_simulation = ?1
       and position.timestep = ?2""", (simulationID,timestep)) #and star.isCluster=0 LIMIT 10000000 OFFSET 10000000 and star.isCluster=1
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
    #print(8000 * np.tan(0.5*0.002))
    simulationID = 1
    database = r"E:\Master_Thesis\VS_Project\N_Body\Output\Database\Default.db"
    output = r"E:\Master_Thesis\VS_Project\N_Body\Output\NGC2244"# + str(simulationID)
    # create a database connection
    conn = createConnection(database)
    with conn:
        data = select2dPositions(conn, simulationID, 0)
        #data[:,1] -= np.average(data[:,1]) #average a around 0
        #data[:,2] -= np.average(data[:,2]) #average d around 0
        print(min(data[:,1]),max(data[:,1]))
        print(min(data[:,2]),max(data[:,2]))
        source = makeSource(data)
        #vega = sim.source.source_templates.vega_spectrum(mag=10)
        #ab_spec = sim.source.source_templates.ab_spectrum(mag=20)
        #tbl = Table(names=["x",   "y",    "ref",  "weight"],
        #                  data=data)
        #print(,np.std(tbl["x"]))
        #tbl["x"]-=np.average(tbl["x"])
        #tbl["y"]-=np.average(tbl["y"])
        #tbl["ref"]=np.array(tbl["ref"]).astype(int)
        #table_source = sim.Source(table=tbl, spectra=[vega]) # why spectra=[vega, ab_spec] not working?

        #print(table_source)
        #print(source.fields)
        #print(source.spectra)

        #sim.server.download_package("telescopes/LFOA")
        if 1==0: #micado
            sim.server.download_package(["locations/armazones",
                                 "telescopes/ELT",
                                 "instruments/MICADO",
                                 "instruments/MAORY"])
            micado = sim.OpticalTrain("MICADO")
            micado.cmds["!OBS.dit"] = 36000 # grid integration time -> länge der beobachtung in s
            #micado.cmds["!OBS.ndit"] = 1   # #beobachtungen
            #micado.cmds["!INST.filter_name"] = "Ks" #-> optischer filter
            #micado["scope_vibration"].include = False
            #micado["detector_window"].include=False #nur ausschnitt
            #micado["full_detector_array"].include=False

            micado.update()

            micado.observe(source)
            hdus_micado = micado.readout()[0]
            print(hdus_micado)
            #hdus_micado.writeto(output+"test.fits")

            #for i in range(1,10):
            #    plt.subplot(3,3,i)

            im = hdus_micado[1].data
        elif 1==0:
            sim.server.download_package(["locations/Paranal","telescopes/VLT", "instruments/HAWKI"])
            cmds = sim.UserCommands(use_instrument="HAWKI")
            hawki = sim.OpticalTrain(cmds)
            hawki.cmds["!OBS.dit"] = 360 # grid integration time -> länge der beobachtung in s
            hawki.cmds["!OBS.ndit"] = 1   # #beobachtungen
            hawki.update()
            hawki.observe(source)
            hdus_hawki = hawki.readout()[0]
            print(hdus_hawki)
            im = hdus_hawki[1].data

        else: #LFOA
            sim.server.download_package("telescopes/LFOA")
            lfoa = sim.OpticalTrain("LFOA")
            lfoa.cmds["!OBS.dit"] = 360 # grid integration time -> länge der beobachtung in s
            lfoa.cmds["!OBS.ndit"] = 2   # #beobachtungen
            lfoa.update()
            lfoa.observe(source)
            lfoa.readout()[0].writeto(output+"/my_image.fits")
            hdus_lfoa = lfoa.readout()[0]
            print(hdus_lfoa)
            im = hdus_lfoa[1].data

        #plt.imshow(hdus_micado[1].data, norm=LogNorm())#, origin="lower", vmax=1E6

        plt.figure(figsize=(12,12))
        plt.imshow(im, vmin=0.9*np.median(im), vmax=2*np.median(im),  norm=LogNorm())

        plt.show()


        #np.savetxt(output+r"\NGC2244_pos.dat", data2save, delimiter=',')

        #for i in np.unique(data[:,0]):
        #    timestepData = data[data[:,0] == i]
            #plt.scatter(timestepData[:,1], timestepData[:,2])
            #plt.show()
        #    if(i<9 or (i%4==0 and i<52) or (i%13==0 and i>52)):
        #        a = numpy.append(a, a[0])
        #        np.savetxt(output+r"\NGC2244_pos_"+str(int(i))+'.dat', timestepData[:,1:], delimiter=',')   # X is an array

if __name__ == '__main__':
    main()
