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

Galactocentric Cylindrical (GCY)
--------------------------------

Local Standard of Rest (LSR)
----------------------------

Heliocentric Cartesian (HCA)
----------------------------

Heliocentric Galactic Polar (HGP)
---------------------------------

Heliocentric Equatorial Polar (HEQ)
-----------------------------------
