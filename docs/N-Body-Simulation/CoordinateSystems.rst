==================
Coordinate Systems
==================

In the simulation, locations and velocities of stars are stored in galactocentric cartesian coordinates (GCA).
However, the observer/scopesim expect heliocentric equatorial polar coordinates.
Moreover, observational data is usually given in heliocentric galactic polar coordinates.
This data is used for initial cluster positions in the simulation and to compare results.
The implementation of transformations is therefor inevitable. The code has been adapted from GalPot (:cite:`McMillan_2016`).

Galactocentric Cartesian (GCA)
------------------------------

GCA is a right-handed coordinate system with the galactic center in its origin.
The projection of :math:'\hat{e}_x' onto the galactic equator (or midplane) points to the initial location of the sun and
'\hat{e}_z' towards the galactic north pole. Therefor, the direction of galactic rotation at the location of the sun is the negative z axis.

Units: positions [pc], velocities [km/s]

Galactocentric Cylindrical (GCY)
--------------------------------

GCY is similar to GCA. Positions and velocities between the two systems goes as follows.




Local Standard of Rest (LSR)
----------------------------

Heliocentric Cartesian (HCA)
----------------------------

Heliocentric Galactic Polar (HGP)
---------------------------------

Heliocentric Equatorial Polar (HEQ)
-----------------------------------
