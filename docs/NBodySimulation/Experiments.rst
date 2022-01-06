===========
Experiments
===========

Parameter optimization
----------------------

.. _fig-DBSCAN:
.. figure:: Images/DBSCAN_parameter_space.*
    :align: center

    Precision depending on DBSCAN parameters

The quality of the cluster analysis with DBSCAN depends on the choice of its parameters.
:numref:`fig-DBSCAN` depicts the precision for a :math:`1 [kM_\odot]` cluster :math:`300 \textup{[pc]}` from the GC.
Decreasing :math:`\epsilon` leads to a higher precission until it is too small, at which point multiple clusters are detected.
Based on this observation :math:`\epsilon` was set to :math:`1.1*10^{-5}` and nPoints to 200 for all simulations.

Setup
-----

A total of 25 simulations with varying focus and cluster mass were carried out in order to study the effects of these parameters on the cluster detection performance.
The focus was set to :math:`l=0^{\circ}` and :math:`b \in \left \{0,5,10,25,180 \right \}[^\circ]` HGP and the cluster mass :math:`\in \left \{0.64, 1.6, 4.0, 10, 25 \right \} [\textup{kM}_\odot]`.
Each simulation was repeated 10 times for error estimation.

.. _25_parameters-label:

Parameters
^^^^^^^^^^

The following parameters remained unchanged between the simulations.

**General simulation parameters**

* | FOV angle: :math:`54 \textup{[arcsec]}`
  | The field of view angle, chosen large enough to cover most of the cluster stars. See :math:`\alpha` in :numref:`fig-cone`
* | View distance: :math:`9 \textup{[kpc]}`
  | The height of the COV or line of sight distance reaching behind the GC when looking towards it.
* | Cluster distance: :math:`8 \textup{[kpc]}`
  | The distance between the observer and the cluster. When looking straight at the GC the distance between the cluster and the GC is :math:`300 \textup{[pc]}` 
  | The mean cluster velocity is set to the circular velocity at this location.
* | View point: :math:`\begin{bmatrix}8300 & 0 & 27\end{bmatrix}^{T}_{GCA} \textup{[pc]}`
  | The position of the observer in GCA.
* | Timestep: :math:`28 \textup{[day]}`
  | Amount of time between the two recorded snapshots. The time per integration is :math:`7 \textup{[day]}`. Therefor snapshots are taken every 4 integrations.

**Cluster specific parameters (McLuster)**

* (P) Profile: 0 Plummer density profile
* | (R) Radius: -1
  | With this setting the radius is calculated by mcluster via a mass to half-mass radius relation as described in :cite:`Marks2012`
* | (Q) Virial ratio: 0.5 
  | The cluster is in virial equilibrium.
* (f) IMF: 1 Kroupa ranging from 0.08 Msun to 100 Msun
* (C) Output: 3 The resulting mass, position and velocity for each star is written into a file.

Results and Interpretation
--------------------------

.. _fig-25_n_stars:
.. figure:: Images/25_n_stars.*
    :align: center

    Number of simulated/detected cluster and field stars at :math:`10^\circ` depending on the cluster masses

Depending on the total cluster mass the amount of simulated CS (SCS) ranges from :math:`1.3e3` to :math:`40.4e3` while the amount of mapped CS (MCS) lies between :math:`1.0e3` and :math:`16.4e3`.
The decrease in detectability of CS is due to increasing CS density and has no direct impact on the clustering performance.
The difference between observed stars (OS) and mapped stars (MS) is negligible, in fact the amount of not mapped stars decreases with increasing number of SCS.

While the number of simulated FS (SFS) remains constant, the number of mapped FS (MFS) decreases with increasing number of CS because bright CS hide FS.
This inverse relationship does impact the clustering performance to some degree, less MFS means less potential FPs as well as TNs, the latter are not relevant for Precision and F1 score.

.. _fig-25_precision:
.. figure:: Images/25_precision.*
    :align: center

    Precision for different mass bins, angles and cluster masses

.. _fig-25_precision_sim:
.. figure:: Images/25_precision_sim.*
    :align: center

    Precision using accurate velocities

:numref:`fig-25_precision` displays the precision using the velocity of OS and :numref:`fig-25_precision_sim` for the velocity of SS.
The only relevant difference between simulated and observed HTP proper motion is the accuracy of position and consequently of velocity.
While both figures show the same relationships, the drop in overall performance due to inaccuracies introduced during observation and source detection are painfully apparent.

.. _fig-25_F1:
.. figure:: Images/25_F1.*
    :align: center

    F1 score for different mass bins, angles and cluster masses


As is discernable in :numref:`fig-25_precision` and :numref:`fig-25_precision_sim` the precision, with one exception, is correlated with the angle.
Curiously for the same cluster mass the precision is lower at :math:`10^\circ` than at :math:`5^\circ`.

.. _fig-25_avg_vel_640:
.. figure:: Images/25_avg_vel_640.*
    :align: center

    Average cluster and field star velocity at different angles

:numref:`fig-25_avg_vel_640` provides the explanation for this outlier. 
At :math:`10^\circ` the average field star velocity is closer to the average cluster star velocity than at any other angle, making it harder to differentiate between cluster and noise.

The bigger the cluster mass, the higher the cluster star velocity density, which implies the second correlation - precision with cluster mass - presented in :numref:`fig-25_precision`.

.. _fig-25_vel_scatter:
.. figure:: Images/25_vel_scatter.*
    :align: center

    2D HTP velocity of simulated clusters

:numref:`fig-25_vel_scatter` displays examples for the 2D HTP velocity space of two simulated clusters near the GC.

In this example the :math:`0.64 [kM_{\odot}]` cluster only has 1143 stars inside the circle while :math:`10 [kM_{\odot}]` has 3158.
In both cases statistically the same amount of field stars fall within that area, leading to a higher ratio of FPs and therefore a lower precision for the lower mass cluster.

