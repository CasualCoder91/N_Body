=================
Cluster analysis
=================

The objective is to detect all cluster stars present in the observed data with a certain accuracy.
One image contains a set of stars :math:`s_{i}` with positions :math:`\vec{x}_{i}\left ( t \right )` in pixel coordinates and apparent magnitudes :math:`m_{i}`.
This information is not enough to detect cluster memberships. The number density of cluster stars in position space  is not high enough relative to the field stars.
Moreover, depending on the distance of the cluster to the observer and the cluster size, cluster stars often cover a significant area of if not the whole image.
The proper motion :math:`\vec{v_{i}}(t) density of cluster stars on the other hand is very high relative to that of field stars.

Velocity approximation
----------------------

The proper motion of a star can be approximated from its position at two different timestamps via

.. math::
    v_{i} = \frac{x_{i}\left ( t+dt \right )-x_{i}\left ( t \right )}{dt}

However, this is only possible if the positions at two different timesteps can be accurately associated with the same star :math:`s_{i}`.
The assumtion that :math:`\vec{x}(t+dt)` belongs to the same star as :math:`\vec{x}(t)` if their Euclidian distance is smaller than the distance of :math:`\vec{x}(t)` to any other :math:`\vec{x}(t+dt)`,
only hold true for small :math:`dt`. However, due to the discreteness of images, :math:`dt` has to be large enough so that the change in position is detectable between timesteps. 
The pixel scale - the ratio of arcsec to pixel - dictates a lower bound for :math:`dt`.

Since the data results from a simulation, it is trivial to verify whether or not the attribution of the two positions is indeed correct.
In fact, tests using only the Euclidian distance as a metric with :math:`dt = 1 day` lead to a small but significant amount of wrong assignments.

An additional condition was introduced. Stars at :math:`\vec{x}(t)` and :math:`\vec{x}(t+dt)` have to have an similar apparent magnitude :math:`m` so be considered the same star.

.. math::
    \left | 1-\frac{m_{i}\left ( t \right )}{m_{j}\left ( t+dt \right )} \right | < \varepsilon_{m}

where :math:`\varepsilon_{m}` is the maximum relative difference in apparent magnitude.

With this constraint, all positions where correctly assigned during further tests.

DBSCAN
------

The clustering algorithm of choice is DBSCAN since it is density based and able to detect clusters of arbitrary shape. 
Additionaly and contrary to other algorithms the amount of clusters to find is not a parameter. DBSCAN detects any clusters present in the data based on two parameters:

1. :math:`\epsilon`: the maximum distance between points to be considered neighbors
1. nPoints: the minimum amount of neighbors for a point to be classified as core point.

During excecution all datapoints are classified as one of the following:

1. core point: a point with at least nPoints within :math:`\epsilon`:
2. border point: a point having at least one core point but less than nPoints within :math:`\epsilon`:
3. noise/outlier: any other point

The implementation of DBSCAN can be summarized as follows: Iterate the list of points. If the current point is not already classified, check if it meets the requirements to be classified as core point.
Once a core point has been found, the neighboring points of that point are tested. If they too have enough neigbors the recursion continues untill all neighbors are classified as eighter core or border points. 

In a naive implementation, the distance of each point to every other point is checked. The time complexity of such an implementation is :math:`O(n^2)`.
Moreover, for large datasets the recursion can lead to stack overflow.

R*-tree and other data structures can be used to improve the performance to an average of :math:`O(n\log{n})` (:cite:`Ester1996`)
The library mlpack (:cite:`Curtin2018`) includes an implementation of DBSCAN supporting R*-tree and many other trees.

Initially and in addition to the standard condition :math:`\epsilon_{x}` regarding the spatial distance between points an supplementary condition :math:`\epsilon_{v}` was introduced.
The difference in velocity between two stars has to be smaller than :math:`\epsilon_{v}` to be deemed neighbors:

.. math::
    \left \|\vec{v}_{1}-\vec{v}_{2}  \right \|_{2}< \epsilon_{v}

Large :math:`\epsilon_{x}` lead to more accurate membership detection. It turned out, the spatial distance condition does not benefit the results at all and was dropped.
