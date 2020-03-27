Tasks
=====
An overview of the tasks that we aim to complete:

Observations
------------
1. star catalogue of simulated open clusters
    - N-body simulations of cluster in free space
    - N-body simulations of cluster in galactic potential
    - Determine realistic initial conditions
    - Set realistic boundary conditions
    - Determine time span to look at - [0.1, 1, 10] years?
    - Adjust the mass and velocity dispersion of the cluster to fit the observations

2. star catalogue for fore- and background stars in the galaxy
    - function for generating number of stars in a volume (dx, dy, dz)
      at any galactic position (x,y,z)
    - use a present-day mass function for the field stars
    - find plausible values for velocity and dispersion and orbital parameter distribution

3. Images of clusters and on-sky scenes
    - Images/Movies of clusters evolving over time in free space
    - Images/Movies of field star population (evolving over time)
    - Images/Movies of clusters evolving over time with contaminating fore- and
      background stars

4. Mock observations of clusters
    - "Observations" of clusters evolving over time including all optical
      artifacts from atmosphere, telescope, instrument, etc
    - "Observe" this cluster with the ELT using SimCADO
    - Extra: Make several observations with different telescopes (ELT, VLT, Hubble, JWST, etc)

Data Reduction
--------------
1. Extract the stars and determine relative position of the stars for each timestep
    - By using the known position
    - By using a blind extraction algorithm (e.g. starfind)
    - By using aperture and/or PSF photometry


2. Track cluster members
    - Write function that tracks motion of extracted stars

Analysis
--------
1. Observational limits
    - Extract the star properties (luminosity, position, velocity)
    - Determine the percentage of stars that have detectable motions in each snapshot
    - **Find the time (snapshot) where >90% of the stars have detectable motions with respect to t=0**
    - **Find the lowest recoverable mass for a given set of observing conditions**

2. Cluster membership
    - Define a criteria and/or method for determining probability of cluster membership
    - Determine cluster members. Compare this list to the input list
    - **Determine percentage of false positives, false negatives based on velocity**

3. Get IMF for cluster members
    - Isolate only the cluster members
    - **Find mass bins where 90% of input stars have cluster membership**
    - **Find mass bins where 90% of input stars are detected**
    - Reconstruct IMF for cluster based on stars with >50%, >85%, >95% membership probabilites
    - **Compare reconstructed IMF to input IMF**
