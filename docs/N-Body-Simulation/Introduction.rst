============
Introduction
============

The Extremely Large Telescope (ELT) is currently under construction. 
This first next generation 40m class telescope will have the resolution and sensitivity needed to study the motions of individual stars in the galactic centre (GC).
These position and motions can be used to detect star clusters using some type of clustering algorithm and consequently estimate their Initial Mass Function (IMF).
Whether the IMF of star clusters is universal is subject of discussion :cite:`Bastian2010` at this time.
The study of Young Massive Clusters (YMCs) :cite:`Zwart2010` near the GC will hopefully give more insight into this hotly debated topic.
YMCs are tightly grouped clusters of stars, generally containing >10000 members. 
All original members are still present, the IMF is well sampled in all mass regimes. These facts combine to give a good picture of the end product of a star formation event.

Given detected star clusters the question remains, how reliable these results are.
If the true classification of the studied set of stars is known, the performance of the clustering algorithm and hence the reliability of the results can be calculated.
However, this is generally not the case.

In this master thesis an N-body simulation containing cluster stars and field stars under the influence of the milky way potential is performed.
Snapshots are taken at different timesteps and fed to ScopeSim :cite:`Leschinski2020` to create mock observations.
With the help of Photutils :cite:`Bradley2021` positions are extracted and the DBSCAN algorithm used to detect cluster and field stars.
Finally, the reliability of the results can be determined by comparing its results with the initially simulated stars.

