==============
Data Reduction
==============

With the help of Photutils :cite:`Bradley2021` the 2D HTP positions and fluxes of stars are detected from the FITS files generated during observation.
In a first step all stars are stored in the database, each having exactly one location corresponding to the timestep of the FITS file.

The following method, while yielding decent results, is certainly not optimal.
Due to hardware and time constraints, options such as Image Segmentation were not feasible. Moreover, the parameters of the chosen method were not optimized beyond some spot checks.

Detecting Stars
---------------

PSFs of bright stars were wrongfully detected as stars. Increasing the detection threshold would have meant not detecting real faint stars in different areas.
Hence the decission was made to mask square areas around bright stars leading to only one detection within those areas.

The side length :math:`l_i` of the square depends on the flux :math:`F_i` of the stars: :math:`l_i = f(F_i)`.

To find the approriate function, FITS for single stars at a given distance and with varying mass were generated, 
the sources in each file detected using DAOStarFinder and their flux and maximum distance between the correct and any wrong sources calculated.
A linear fit of the resulting dataset :ref:`appendix-label` lead to the "empirical" function

..  math::
    l_i = \begin{cases}
     & 0\text{ if } F_i < 100\\ 
     & 0.01*F_i+28\text{ if } F_i \geqslant 100\\ 
    \end{cases}

ToDo: Add 3 images "without mask" "drawn mask" "after mask"

The DAOStarFinder method is called twice:
1. To find the bright stars and generate the mask. The resulting table contains one row for each source. 
This table is sorted by the flux column in descending order and iterated from top to bottom until the current entry has :math:`F_i < 100`. 
Elements of the mask - a 2D boolean array with the same size as the image - is updated according to function (?) and the current entry stored in a new table if located outside a masked area.
2. passing the mask parameter generated in the previous step and returning sources outside the masked areas.

Both the bright sources recorded after the first and the faint sources returned from second call are stored in the database.
