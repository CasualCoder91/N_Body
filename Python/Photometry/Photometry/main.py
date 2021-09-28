import numpy as np
import sys

from Photutils import pu_all, find_mask_borders
from ScopeSim import ss_all, make_fits
import config # import simulation_id, database_path, fits_path, output_path, timestep, n_pixel, save_img, exposure_time #my "global variables"


def all():
    if len(sys.argv) > 1:
        config.simulation_id = int(sys.argv[1])
        config.save_img = False
    ss_all()
    pu_all()


def fake_star_test():
    for m in np.linspace(1,0,50,False):
        data = np.array([[7000,0,0,m],[9000,0,0,0.01]]) #todo: loop for different masses
        make_fits(data,n_pixel=1000)
        flux, dist = find_mask_borders()
        print(flux, dist)


all()