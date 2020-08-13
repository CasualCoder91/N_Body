J/A+A/468/151       Nearby open clusters radii and masses      (Piskunov+, 2007)
================================================================================
Towards absolute scales for the radii and masses of open clusters.
    Piskunov A.E.,  Schilbach E., Kharchenko N.V., Roeser S., Scholz R.-D.
   <Astron. Astrophys. 468, 151 (2007)>
   =2007A&A...468..151P
================================================================================
ADC_Keywords: Clusters, open ; Morphology
Keywords: Galaxy: open clusters and associations: general -
          solar neighbourhood - Galaxy: stellar content

Description:
    The table presents tidal parameters and masses of 236 open clusters in
    the nearest kiloparsecs around the Sun. The parameters are derived
    from a fitting of three-parameter King profiles to the the observed
    density distribution. The clusters are sub-sample of the Catalogue of
    Open Cluster Data and their members are selected from high precision
    homogeneous all sky catalogue ASCC-2.5 (Kharchenko, 2001, Cat.
    <I/280>). Up to four cluster membership samples in wide cluster area
    are considered for every cluster and the best solution is presented in
    the table. Clusters in the table are sorted by numbers in the COCD
    order.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table2.dat    110      236   King parameters and tidal masses for 236 open
                              clusters
--------------------------------------------------------------------------------

See also:
       I/280 : All-sky Compiled Catalogue of 2.5 million stars (Kharchenko 2001)
 J/A+A/438/1163 : Catalogue of Open cluster Data (COCD) (Kharchenko+, 2005)
 J/A+A/440/403  : 109 new Galactic open clusters (COCD1) (Kharchenko+, 2005)

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units      Label  Explanations
--------------------------------------------------------------------------------
   1-  4  I4    ---        COCD   COCD number (>1000 for COCD1)
   6- 22  A17   ---        Name   Cluster designation (NGC, IC or other)
  24- 29  F6.2  deg        GLON   Galactic longitude
  31- 36  F6.2  deg        GLAT   Galactic latitude
  38- 41  I4    pc         Dist   Heliocentric distance
  43- 46  F4.2  mag        E(B-V) Colour-excess
  48- 51  F4.2  [yr]       logt   Logarithm of age
  53- 56  F4.2  deg        r1     Empirical angular radius of the core
  58- 61  F4.2  deg        r2     Empirical angular radius of the cluster
      63  I1    ---        Sol    [1-4] Number of the acceptable solutions of
                                    King's parameters from 4 membership groups
      65  I1    ---        Sam    [1-4] Number of the membership sample
                                    providing the best solution
  67- 69  I3    ---        N2     Number of the cluster members of this sample
                                    within r2
  71- 74  F4.1  pc         rc     Core radius (King)
  76- 78  F3.1  pc       e_rc     rms error of King core radius
  80- 83  F4.1  pc         rt     Tidal radius (King)
  85- 88  F4.1  pc       e_rt     rms error of King tidal radius
  90- 93  F4.1  ---        k      King's normalizing factor
  95- 98  F4.1  ---      e_k      rms error of King normalizing factor
 100-104  F5.3 [solMass]   logM   Logarithm of cluster mass
 106-110  F5.3 [solMass] e_logM   rms error of logarithm of cluster mass
--------------------------------------------------------------------------------

Acknowledgements: Anatoly Piskunov
================================================================================
(End)                                        Patricia Vannier [CDS]  23-Mar-2007
