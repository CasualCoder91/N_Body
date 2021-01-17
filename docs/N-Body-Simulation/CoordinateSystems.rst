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
The projection of :math:`\hat{e}_x` onto the galactic equator (or midplane) points to the initial location of the sun and
:math:`\hat{e}_z` towards the galactic north pole. Therefor, the direction of galactic rotation at the location of the sun is the negative z axis.

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

Like GCA, LSR is a right-handed coordinate system.
The origin of positions is the location of the sun
and the origin of velocity is the velocity of a star on a circular orbit with mean velocity of stars in the solar neighborhood.
:math:`\hat{e}_x` points towards the galactic center, :math:`\hat{e}_y` towards the direction of galactic rotation
and :math:`\hat{e}_z` roughly towards the galactic north pole.

The location of the sun is approximately given by (:cite:`McMillan_2016_2`, :cite:`Binney_1997`)

.. math::
    \vec{x}_{sun,GCA} \cong  (8.20,0,0.014)[kpc]

and the mean velocity (:cite:`McMillan_2016_2`)

.. math::
    \vec{v}_{mean,GCA} \cong  (0,-232.8,0)[km/s]

Since the sun is generally not in the galactic midplane, there is an angle between the planes spanned by :math:`(x,y)_{GCA}` and :math:`(x,y)_{LSR}`.
This angle can be expressed as

.. math::
    \sin(\alpha) = \frac{z_{sun,GCA}}{\sqrt{x^2+z^2}} \\
    \cos(\alpha) = \frac{x_{sun,GCA}}{\sqrt{x^2+z^2}}

The transformation of position and velocity vectors from GCA to LSR has to contain a rotation by :math:`-\alpha` about the y axis.

.. math::
    x_{LSR} = \cos(\alpha)( x_{sun,GCA} - x_{GCA} ) - \sin(\alpha)(z_{GCA}-z_{sun,GCA}) \\
    y_{LSR} = -y_{GCA} \\
    z_{LSR} = \sin(\alpha)(x_{sun,GCA} - x_{GCA}) + \cos(\alpha)( z_{GCA} - z_{sun,GCA}) \\ \\
    u_{LSR} = -\cos(\alpha)u_{GCA} - \sin(\alpha)w_{GCA} \\
    v_{LSR} = v_{sun,GCA}-v_{GCA} \\
    w_{LSR} = -\sin(\alpha)u_{GCA} + \cos(\alpha)w_{GCA} \\

Heliocentric Cartesian (HCA)
----------------------------

Heliocentric Galactic Polar (HGP)
---------------------------------

Heliocentric Equatorial Polar (HEQ)
-----------------------------------
