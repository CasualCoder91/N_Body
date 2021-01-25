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

The only distinction between HCA and LSR is the origin of the velocity. In HCA the origin is the velocity of the sun.
The difference between the average velocity in the solar neighborhood and the sun itself, also called the peculiar motion of the sun,
is given by (:cite:`2010MNRAS.403.1829S`)

.. math::
    \vec{v}_{sun,LSR} \cong  (11.1, 12.24, 7.25)[km/s]

Transformation from LSR to HCA can be written as

.. math::
    \vec{x}_{HCA} = \vec{x}_{LSR} \\
    \vec{v}_{HCA} = \vec{v}_{LSR} - \vec{v}_{sun,LSR}

Heliocentric Galactic Polar (HGP)
---------------------------------

HGP is a spherical coordinate system with origins for position and velocity identical to those defined in HCA.
Coordinates given in this system are often called galactic coordinates.
The galactic longitude :math:`l` is the angular distance in the galactic midplane with :math:`l=0^{\circ}` towards the galactic center.
The galactic latitude :math:`b` denotes the angle below and above the galactic midplane ranging from :math:`-90^{\circ}` to :math:`90^{\circ}`.
and :math:`r` the radial distance.

The transformation from HCA to HGP is identical to the transformation from GCA to GCP (?) with :math:`l \equiv \varphi`, :math:`b \equiv \theta` and :math:`\rho \equiv r`

Heliocentric Equatorial Polar (HEQ)
-----------------------------------

HEQ, like HGP, is a spherical coordinate system having the same origins for position and velocity.
However, angles are given in and normal to the celestial equator which is not parallel to the galactic midplane.
The right ascension (:math:`a`) is the angular distance in the equator with :math:`a=0^{\circ}` towards the northward equinox.
The equinox is the intersection of the ecliptic - the plane in which the earth orbits the sun - and the celestial equator.
The declination (:math:`d`) is the angular distance above or below in the equator.

Since the ecliptic and the equator are in motion, a reference frame is needed.
A reference frame consists of quantities defining the coordinate system at a specific time as well as methods to
calculate those quantities for any other date. A commonly used reference frame is defined for the J2000.0 epoch (:math:`\epsilon_0`).

In order to transform between HCA and HEQ at :math:`\epsilon_0`,
the direction of the north Galactic pole (NGP) and the galactic center (GC) are needed in both basis.

In HCA the NGP is simply :math:`\vec{x}_{NGP,HCA}=(0, 0, 1)`.
In HGP, since the direction is normal to the fundamental plane, :math:`b=90^{\circ}_{GC,HGP}`.
In HEQ at :math:`\epsilon_0` the direction is

.. math::
    a_{NGP} = 12^h51^m26.28^s \\
    d_{NGP} = 27^{\circ}7^\prime41.7^{\prime\prime}

The GC defines the x axis of HCA: :math:`\vec{x}_{GC,HCA}=(1, 0, 0)`. In GC the same direction is

.. math::
    a_{GC,HEQ} = 17^h45^m40.0409^s \\
    d_{GC,HEQ} = -29^{\circ}0^\prime28.118^{\prime\prime}

To express these basis vectors in HCA basis, they can to be transformed as follows

.. math::
    x_{HCA} = \cos(d)\cos(a) \\
    y_{HCA} = \cos(d)\sin(a) \\
    z_{HCA} = \sin(d) \\

The third basis vector is the cross product of :math:`\vec{x}_{NGP}` and :math:`\vec{x}_{GC}`.
With these basis vectors the change of basis matrix is

.. math::
    M = \left [\hat{e}_x,\hat{e}_y,\hat{e}_z\right ]

The full transformation from HCA to HEQ consists of the two steps: the multiplication with :math:`M` followed by
the transformation from cartesian to spherical as given in section (?) GCP.

For the transformation between HGP and HEQ the direction of the north celestial pole (NCP) is required .
NCP is perpendicular to the celestial equator, hence :math:`d_{NGP} = 90^{\circ}`.
In HGP at :math:`\epsilon_0`, NCP is

.. math::
    l_{NCP} = 123^{\circ}55^\prime55.2^{\prime\prime}\\
    b_{NCP} = 27^{\circ}7^\prime41.7^{\prime\prime}

Using NGP and NCP the transformation from HGP to HEQ at :math:`\epsilon_0` is

.. math::
    \sin(d) = \sin(d_{NGP})\sin(b) + \cos(d_{NGP})\cos(b)\cos(l_{NCP}-l) \\
    \cos(d)\sin(a-a_{NGP}) = \cos(b)\sin(l_{NCP}-l) \\
    \cos(d)\cos(a-a_{NGP}) = \cos(d_{NGP})\sin(b)-sin(d_{NGP})\cos(b)\cos(l_{NGP}-l)

Three angles describe the precision of both planes between a epoch :math:`\epsilon_0` and the date of observation :math:`\epsilon_D`.

(...)

With these three rotations a precession matrix :math:`P` as well as its inverse can be formalized.

(...)
