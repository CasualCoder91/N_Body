==================
Cluster analysis
==================

The goal is to detect all clusters present in the observed data with a certain accuracy.
The observed data points are positions (galactic coordinates: l and b) and mass at certain timestamps.
The proper motion of a star can be approximated from its position at two different timestamps via

.. math::
    v_{i} = \frac{x_{i}\left ( t+dt \right )-x_{i}\left ( t \right )}{dt}

However, this is only possible if the positions at two different timesteps can be accurately associated with one specific star.
Since the data results from a simulation, it is trivial to verify whether or not the attribution of the two positions is indeed correct.
:math:`\vec{x}(t+dt)` belongs to the same star as :math:`\vec{x}(t)` if their Euclidian distance is smaller than the distance of :math:`\vec{x}(t)` to any other :math:`\vec{x}(t+dt)`.
