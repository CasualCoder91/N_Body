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

Galactocentric Polar (GCP)
--------------------------------

GCP is a spherical coordinate system and similar to GCA.
Position and velocity transformation between the two systems goes as follows.

.. math::
    \rho  = \sqrt{x^2+y^2+z^2} \\
    \theta = \textup{atan2}{\left ( y,x \right )}\\
    \varphi  = \arcsin\left ( \frac{z}{\sqrt{x^2+y^2+z^2}} \right ) \\
    \dot{\rho} =  \frac{x\dot{x}+y\dot{y}+z\dot{z}}{\sqrt{x^2+y^2+z^2}}\\
    \dot{\theta} = \frac{\dot{x}y-x\dot{y}}{x^2+y^2} \\
    \dot{\varphi} = \frac{z(x\dot{x}+y\dot{y})-\dot{z}(x^2+y^2)}{(x^2+y^2+z^2)\sqrt{x^2+y^2}}


Local Standard of Rest (LSR)
----------------------------

Heliocentric Cartesian (HCA)
----------------------------

Heliocentric Galactic Polar (HGP)
---------------------------------

Heliocentric Equatorial Polar (HEQ)
-----------------------------------
