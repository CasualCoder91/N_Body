==================
Cluster analysis
==================

The goal is to detect all clusters present in the observed data with a certain accuracy.
The observed data points are positions (galactic coordinates: l and b) and mass at certain timestamps.
The proper motion of a star can be approximated from its position at two different timestamps via

.. math::
    v_{i} = \frac{x_{i}\left ( t+dt \right )-x_{i}\left ( t \right )}{dt}

However, this is only possible if the positions at two different timesteps can be accurately associated with one specific star.
:math:`\vec{x}(t+dt)` belongs to the same star as :math:`\vec{x}(t)` if their Euclidian distance is smaller than the distance of :math:`\vec{x}(t)` to any other :math:`\vec{x}(t+dt)`.
Since the data results from a simulation, it is trivial to verify whether or not the attribution of the two positions is indeed correct.

The clustering algorithm of choice is DBSCAN since it can detect clusters of arbitrary shape and the amount of clusters is not a parameter.
In addition to the usual condition :math:`\epsilon_{x}` regarding the spatial distance between points a additional condition :math:`\epsilon_{v}` is introduced.
The difference in velocity between two stars has to be smaller than :math:`\epsilon_{v}` to be deemed neighbors:

.. math::
    \left \|\vec{v}_{1}-\vec{v}_{2}  \right \|_{2}< \epsilon_{v}
