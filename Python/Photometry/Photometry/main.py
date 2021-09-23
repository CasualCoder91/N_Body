import numpy as np

from Photutils import pu_all, find_mask_borders
from ScopeSim import ss_all, make_fits


def all():
    ss_all()
    pu_all()


def fake_star_test():
    for m in np.linspace(1,0,50,False):
        data = np.array([[7000,0,0,m],[9000,0,0,0.01]]) #todo: loop for different masses
        make_fits(data,n_pixel=1000)
        flux, dist = find_mask_borders()
        print(flux, dist)


all()