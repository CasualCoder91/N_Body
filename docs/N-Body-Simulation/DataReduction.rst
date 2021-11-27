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
Hence the decission was made to mask square areas around bright stars and hence only detect one star within those areas.

The side length :math:`l_i` of the square depends on the flux :math:`F_i` of the stars: :math:`l_i = f(F_i)`.

To find the approriate function, FITS for single stars at a given distance and with varying mass were generated, 
the sources in each file detected using DAOStarFinder and their flux and maximum distance between the correct and any wrong sources calculated.
A linear fit of the resulting dataset :ref:`appendix-label` lead to the "empirical" function

..  math::
    l_i = \begin{cases}
     & 0\text{ if } F_i < 100\\ 
     & 0.01*F_i+28\text{ if } F_i \geqslant 100\\ 
    \end{cases}


