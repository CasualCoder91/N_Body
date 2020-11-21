import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.coordinates import Angle

ra=Angle('6d31m55.0s')
dec=Angle('04d56m30.0s')

#NGC 6611
print('NGC 6611')
gc = SkyCoord(l=16.94*u.degree, b=0.8*u.degree, frame='galactic', distance=1.719*u.kpc)
print('Galactic: ',gc)
gc = gc.transform_to(coord.Galactocentric) #-6.65563328, 0.50081482, 0.04565517
print('Galactocentric: ', gc)
gc = gc.transform_to('icrs')
print('ICRS: ', gc.to_string()) #gc.to_string(style='hmsdms')

#NGC 2244
print('NGC 2244')
gc = SkyCoord(l=206.30*u.degree, b=-2.08*u.degree, frame='galactic', distance=1.445*u.kpc)
print('Galactic: ',gc)
gc = gc.transform_to(coord.Galactocentric) #-6.65563328, 0.50081482, 0.04565517
print('Galactocentric: ', gc)
gc = gc.transform_to('icrs')
print('ICRS: ', gc.to_string()) #gc.to_string(style='hmsdms')
