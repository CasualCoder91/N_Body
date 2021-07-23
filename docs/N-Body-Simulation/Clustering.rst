==================
Cluster analysis
==================

The goal is to detect all clusters present in the observed data with a certain accuracy.

Velocity approximation
----------------------

The observed data points are positions (galactic coordinates: l and b) and mass at certain timestamps.
The proper motion of a star can be approximated from its position at two different timestamps via

.. math::
    v_{i} = \frac{x_{i}\left ( t+dt \right )-x_{i}\left ( t \right )}{dt}

However, this is only possible if the positions at two different timesteps can be accurately associated with one specific star.
The assumtion that :math:`\vec{x}(t+dt)` belongs to the same star as :math:`\vec{x}(t)` if their Euclidian distance is smaller than the distance of :math:`\vec{x}(t)` to any other :math:`\vec{x}(t+dt)`,
only hold true for small :math:`dt`. However, :math:`dt` has to be large enough so that the change in position is detectable between timesteps. 
The pixel scale - the ratio of arcsec to pixel - dictates a lower bound for :math:`dt`.

Since the data results from a simulation, it is trivial to verify whether or not the attribution of the two positions is indeed correct.
In fact, tests using this assumtion with :math:`dt = 1 day` lead to a small but significant amount of wrong assignments.

An additional condition was introduced. Stars at :math:`\vec{x}(t)` and :math:`\vec{x}(t+dt)` have to have an similar apparent magnitude :math:`m`.

.. math::
    \left | 1-\frac{m_{i}\left ( t \right )}{m_{j}\left ( t+dt \right )} \right | < \varepsilon_{m}

where :math:`\varepsilon_{m}` is the maximum relative difference in apparent magnitude.

With this constraint, all positions where correctly assigned during further tests.

DBSCAN
------

The clustering algorithm of choice is DBSCAN since it can detect clusters of arbitrary shape and the amount of clusters is not a parameter.
In addition to the usual condition :math:`\epsilon_{x}` regarding the spatial distance between points a additional condition :math:`\epsilon_{v}` is introduced.
The difference in velocity between two stars has to be smaller than :math:`\epsilon_{v}` to be deemed neighbors:

.. math::
    \left \|\vec{v}_{1}-\vec{v}_{2}  \right \|_{2}< \epsilon_{v}
