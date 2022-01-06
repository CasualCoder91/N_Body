=================
Cluster analysis
=================

The objective is to detect all cluster stars present in the observed data with a certain accuracy.
One image contains a set of stars :math:`s_{i}` with positions :math:`\vec{x}_{i}\left ( t \right )` in pixel coordinates and apparent magnitudes :math:`m_{i}`.
This information is not enough to detect cluster memberships. The number density of cluster stars in position space is similar to that of the field stars.
Furthermore, cluster and field star positions overlap significantly.
Moreover, depending on the distance of the cluster to the observer, the cluster size and the field of view angle, cluster stars often cover a significant area of if not the whole image.
The proper motion :math:`\vec{v_{i}}(t)` density of cluster stars on the other hand is highly relative to that of field stars and hardly any field stars have :math:`\vec{v_{i}}(t)` close to cluster stars.

Velocity approximation
----------------------

The proper motion of a star can be approximated from its position at two different timestamps via first order Taylor expansion

.. math::
    v_{i} = \frac{x_{i}\left ( t+dt \right )-x_{i}\left ( t \right )}{dt}

However, this is only possible if the positions at two different timesteps can be accurately associated with the same star :math:`s_{i}`.
The assumption that :math:`\vec{x}(t+dt)` belongs to the same star as :math:`\vec{x}(t)` if their Euclidian distance is smaller than the distance of :math:`\vec{x}(t)` to any other :math:`\vec{x}(t+dt)`,
only hold true for small :math:`dt`. However, due to the discreteness of images, :math:`dt` has to be large enough so that the change in position is detectable between timesteps. 
The pixel scale - the ratio of arcsec to pixel - dictates a lower bound for :math:`dt`.

Since the data results from a simulation, it is trivial to verify whether or not the attribution of the two positions is indeed correct.
In fact, tests using only the Euclidian distance as a metric with :math:`dt = 1 day` lead to a small but significant amount of wrong assignments.

An additional condition was introduced. Stars at :math:`\vec{x}(t)` and :math:`\vec{x}(t+dt)` have to have a similar apparent magnitude :math:`m` to be considered the same star.

.. math::
    \left | 1-\frac{m_{i}\left ( t \right )}{m_{j}\left ( t+dt \right )} \right | < \varepsilon_{m}

where :math:`\varepsilon_{m}` is the maximum relative difference in apparent magnitude.

With this constraint, all positions where correctly assigned during further tests.

DBSCAN
------

The clustering algorithm of choice is DBSCAN since it is density based and able to detect clusters of arbitrary shape. 
Additionally and contrary to other algorithms the amount of clusters to find is not a parameter. DBSCAN detects any clusters present in the data based on two parameters:

#. :math:`\epsilon`: the maximum distance between points to be considered neighbors
#. nPoints: the minimum amount of neighbors for a point to be classified as core point.

During execution all data points are classified as one of the following:

#. core point: a point with at least nPoints within :math:`\epsilon`:
#. border point: a point having at least one core point but less than nPoints within :math:`\epsilon`:
#. noise/outlier: any other point

The implementation of DBSCAN can be summarized as follows: Iterate the list of points. If the current point is not already classified, check if it meets the requirements to be classified as core point.
Once a core point has been found, the neighboring points of that point are tested. If they too have enough neighbors the recursion continues until all neighbors are classified as either core or border points. 

In a naive implementation, the distance of each point to every other point is checked. The time complexity of such an implementation is :math:`O(n^2)`.
Moreover, for large datasets the recursion can lead to stack overflow.

R*-tree or similar data structures can be used to improve the performance to an average of :math:`O(n\log{n})` (:cite:`Ester1996`)
The library mlpack (:cite:`Curtin2018`) was used, which includes an implementation of DBSCAN supporting R*-tree and many other trees.

Initially and in addition to the standard condition :math:`\epsilon_{x}` for the spatial distance between points a supplementary condition :math:`\epsilon_{v}` was introduced.
The difference in velocity between two stars has to be smaller than :math:`\epsilon_{v}` for them to be classified as neighbors:

.. math::
    \left \|\vec{v}_{1}-\vec{v}_{2}  \right \|_{2}< \epsilon_{v}

Large :math:`\epsilon_{x}` lead to more accurate membership detection. It turned out, the spatial distance condition does not benefit the results at all and it was dropped.
For larger areas than used here, constraining the spatial distance, for instance via subdivision, should be beneficial.

Performance
-----------

Observed stars are mapped to simulated stars via their proximity in order to measure the performance. If a star is not the closest observed star to any simulated star it remains not mapped.
Not mapped stars exist due to :ref:`background-label` and parts of the PSFs of bright stars being detected as stars. Observed stars can then be attributed one of the following types:

#. True Positive (TP): correctly classified as cluster star.
#. False Positive (FP): wrongly classified as cluster star.
#. True Negative (TN): correctly classified as field star.
#. False Negative (FN): wrongly classified as field star.
#. Unconfirmed Positive (UP): not mapped star classified as cluster star. Treated as FP.
#. Unconfirmed Negative (UN): not mapped star classified as field star. Not taken into account since FN is unlikely.

Accuracy
^^^^^^^^

.. math::
    A = \frac{TP+TN}{TP+TN+FP+FN}

With a large amount of field stars relative to cluster stars, this metric is not ideal as it will give a good rating even if most cluster stars are FN.

Precision and Recall
^^^^^^^^^^^^^^^^^^^^

When FPs are more problematic than FNs, the precision :math:`P` should be high

.. math::
    P = \frac{TP}{TP+FP+UP}

On the flip side, if FNs are a big concern, but FPs tolerable, the recall :math:`R` is a good metric

.. math::
    R = \frac{TP}{TP+FN}

F1 Score
^^^^^^^^

This metric is a balance between :math:`P` and :math:`R`. Contrary to :math:`A` TN is not taken into account.

.. math::
    F_1 = 2 \frac{P*R}{P+R} = \frac{TP}{TP+0.5(FP+UP+FN)}





