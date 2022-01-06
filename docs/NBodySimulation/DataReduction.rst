==============
Data Reduction
==============

With the help of Photutils :cite:`Bradley2021` the 2D HTP positions and fluxes of stars are detected from the FITS files generated during observation.
In a first step, all stars are stored in the database, each having exactly one location corresponding to the timestep of the FITS file.

The following method, while yielding decent results, is certainly not optimal.
Due to hardware and time constraints, options such as Image Segmentation were not feasible. Moreover, the parameters of the chosen method were not optimized beyond some spot checks.

Detecting Stars
---------------

.. _background-label:

Background
^^^^^^^^^^

Testing FITS files generated with no input sources still yield some amount of detected sources.
With the :ref:`25_parameters-label` used for the 25 simulations, 125 sources were detected from an empty background.
The FITS files contain a raster of 64 images, sometimes overlapping and sometimes separated by one pixel due to rounding.
This leads to wrong detections at the corners. However, this effect only explains a fraction of the detections.

.. _masking_bright_stars-label:

Masking bright stars
^^^^^^^^^^^^^^^^^^^^

PSFs of bright stars were wrongfully detected as stars. Increasing the detection threshold would have meant not detecting real faint stars in different areas.
Hence the decision was made to mask square areas around bright stars leading to only one detection within those areas.

The side length :math:`l_i` of the square depends on the flux :math:`F_i` of the stars: :math:`l_i = f(F_i)`.

To find the appropriate function, FITS for single stars at a given distance and with varying mass were generated, 
the sources in each file detected, using DAOStarFinder, and their flux and maximum distance between the correct and any wrong sources calculated.
A linear fit of the resulting dataset (:ref:`appendix-label`) lead to the "empirical" function

..  math::
    l_i = \begin{cases}
     & 0\text{ if } F_i < 100\\ 
     & 0.01*F_i+28\text{ if } F_i \geqslant 100\\ 
    \end{cases}
    :label: side_length

ToDo?: Add 3 images "without mask" "drawn mask" "after mask"

The DAOStarFinder method is called twice:

#. | To find the bright stars and generate the mask. The resulting table contains one row for each source. 
   | This table is sorted by the flux column in descending order and iterated from top to bottom until the current entry has :math:`F_i < 100`.
   | Elements of the mask - a 2D Boolean array with the same size as the image - are updated.
   | All elements inside the box with side length :eq:`side_length` are set to true and the current table entry stored in a new table if located outside a masked area.
#. For passing the mask parameter generated in the previous step and returning sources outside the masked areas.

Both the bright sources recorded after the first and the faint sources returned from second call are stored in the database.
