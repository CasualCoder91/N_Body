===========
Experiments
===========

Setup
-----

A total of 25 simulations with varying focus and cluster mass were carried out and in order to study the effects of these parameters on the cluster detection performance.
Each simulation was repeated 10 times for error estimation.

The following parameters remained unchanged:

Simulation

* FOV angle: :math:`54 \textup{[arcsec]}` \\
The field of view angle, chosen large enough to cover most of the most massive cluster. See :math:`\alpha` in :numref:`fig-cone`
* View distance: :math:`9 \textup{[kpc]}` \\ The height of the COV,
* Cluster distance: :math:`8 \textup{[kpc]}` \\ The distance between the observer and the cluster.
* View point: :math:`\begin{bmatrix}8300 & 0 & 27\end{bmatrix}^{T}_{GCA} \textup{[kpc]}` \\ The position of the observer
* Timestep: :math:`28 \textup{[day]}` \\ Amount of time between the two recorded smapshots. The time per integration is :math:`7 \textup{[day]}`. Therefor snapshots are taken every 4 integations.

Cluster specific parameters (McLuster)

* (P) Profile: 0 Plummer density profile
* (R) Radius: -1 
* (Q) Virial ratio: -0.5 
* (f) IMF: 1 Kroupa
* (C) Output: 3 


five different FOV angles (180°,25,10, 5, 0  ) and cluster masses (0.64, 1.6, 4.0, 10, 25)

ToDo: Continue here (25 observations)