import numpy as np

from Photutils import pu_all, find_mask_borders
from ScopeSim import ss_all, make_fits


def all():
    ss_all()
    pu_all()


def fake_star_test():
    array = np.empty(2)
    for m in np.linspace(10,0,100,False):
        data = np.array([[7000,0,0,m],[9000,0,0,0.01]]) #todo: loop for different masses
        make_fits(data,n_pixel=1000)
        flux, dist = find_mask_borders()
        print(flux, dist)
        array = np.append(array,[flux, dist],axis=0)
    print(array)

fake_star_test()