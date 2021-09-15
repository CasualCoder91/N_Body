import sqlite3
import numpy as np #used in makeSource
import matplotlib.pyplot as plt
import scopesim as sim
sim.rc.__config__['!SIM.sub_pixel.flag'] = True
import pyckles #used for spectra in makeSource
import scopesim_templates as sim_tp
from astropy.table import Table
from astropy import units
from matplotlib.colors import LogNorm

from config import simulation_id, database_path, fits_path, output_path, timestep, n_pixel, save_img, exposure_time #my "global variables"
from database import Database

def make_source(data):
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

def ss_all():

    print("generating fits")

    save_file = True #If True Output fits file

    db = Database() # create a database connection
    #sim.download_package(["instruments/MICADO_Sci",
    #                      "telescopes/ELT", 
    #                      "locations/Armazones", 
    #                      "instruments/MICADO"]) #"MAORY"

    cmd = sim.UserCommands(use_instrument="MICADO_Sci")
    cmd["!DET.width"] = n_pixel
    cmd["!DET.height"] = n_pixel
    cmd["!DET.dit"] = exposure_time #seconds | 3600s=1h

    opt = sim.OpticalTrain(cmd)
    opt["scao_const_psf"].meta["convolve_mode"] = "same"
    opt["scao_const_psf"].meta["rotational_blur_angle"] = 15*exposure_time/3600 # depends on time  ~15°/h (complex function)
    #opt["scao_const_psf"].meta["psf_side_length"] = 1024 #size of diameter in sechseck hinter hellen sternen

    for timestep in [0,1]:
        data = db.select_2d_stars(timestep)

        source = make_source(data)
        opt.observe(source)

        if save_file:
            opt.readout(filename=output_path +"/scopesim_t"+str(timestep)+".fits")

        if save_img:
            fig = plt.figure(figsize=(n_pixel/960, n_pixel/960), dpi=96)
            plt.imshow(opt.image_planes[0].image, norm=LogNorm())
            #plt.colorbar()
            fig.savefig(output_path + "/scopesim_t"+str(timestep)+".png",dpi=960)


if __name__ == '__main__':
    ss_all()

