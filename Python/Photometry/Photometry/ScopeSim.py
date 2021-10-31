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

import config #my "global variables"
from database import Database

import pyckles

def make_source(data):
    """
    Generate a source object from np.array (radius,ascension,declination,mass)

    Returns
    -------
    src: scopesim.Source
    """

    masses = data[:,3]
    distances = data[:,0]
    extinction = data[:,4]

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
    weight = 10 ** (-0.4 * (Mvs + dist_mod + extinction))


    # 8. make table with (x,y,ref,weight)
    tbl = Table(names=["x", "y", "ref", "weight", "masses", "spec_types"],
                data= [ data[:,1],   data[:,2],   ref,   weight,   masses,   spec_types ])

    # 9. make Source with table, spectra
    src = sim_tp.rc.Source(table=tbl, spectra=spectra)

    return src

def ss_all():

    print("generating fits")

    db = Database() # create a database connection
    #sim.download_package(["instruments/MICADO_Sci",
    #                      "telescopes/ELT", 
    #                      "locations/Armazones", 
    #                      "instruments/MICADO"]) #"MAORY"

    for timestep in [0,1]:
        data = db.select_2d_stars(timestep)
        make_fits(data, timestep, True, config.save_img, config.n_pixel)

def make_fits(data,timestep=0,save_file=True,save_img=False,n_pixel=config.n_pixel):

    cmd = sim.UserCommands(use_instrument="MICADO_Sci")
    cmd["!DET.width"] = n_pixel
    cmd["!DET.height"] = n_pixel
    cmd["!DET.dit"] = config.exposure_time #seconds | 3600s=1h

    opt = sim.OpticalTrain(cmd)
    opt["scao_const_psf"].meta["convolve_mode"] = "same"
    opt["scao_const_psf"].meta["rotational_blur_angle"] = 15*config.exposure_time/3600 # depends on time  ~15°/h (complex function)
    #opt["scao_const_psf"].meta["psf_side_length"] = 1024 #size of diameter in sechseck hinter hellen sternen

    source = make_source(data)
    opt.observe(source)

    if save_file:
        opt.readout(filename=config.output_path +"/scopesim_t"+str(timestep)+".fits")

    if save_img:
        fig = plt.figure(figsize=(n_pixel/960, n_pixel/960), dpi=96)
        plt.imshow(opt.image_planes[0].image, norm=LogNorm())
        fig.savefig(config.output_path + "/scopesim_t"+str(timestep)+".png",dpi=960)


def test_make_source():
    """SELECT position.rH, position.aHTP, position.dHTP, star.mass, star.extinction  
           FROM star
           INNER JOIN position on position.id_star = star.id
           where star.id_simulation = ?1
           and position.timestep = ?2
           and position.rH NOT NULL"""
    data = np.array([[8000,0,0,0.5,0],[8100,1,1,1.5,3]])
    make_source(data)

if __name__ == '__main__':
    #test_make_source()
    ss_all()
    #print(sim.__file__)