J/PAZh/21/13    New list of OB associations              (Melnik+, 1995)
================================================================================
New list of OB associations of our Galaxy.
     Melnik A.M., Efremov Yu.N.
    <Pis'ma Astron. Zh. 21, 13 (1995)>
    =1995PAZh...21...13M
================================================================================
ADC_Keywords: Milky Way ; Associations, stellar ; Stars, OB ; Stars, early-type

Abstract:
    A new partition into associations of OB stars within 3kpc of the Sun
    is derived using Battinelli's modification of the cluster analysis
    method. We have found 58 associations, of which only 10% may be
    considered as a random clusters of field stars.

Description:
    We identified 88 associations with at least 5 members which average
    characteristics are shown in Table 1 of the paper by Mel'nik and
    Efremov (1995PAZh...21...13M). The new table (list.dat) exhibits the
    list of 1392 stars -- members of 88 associations of the new partition.
    For each star we give its name, spectral type, V magnitude, galactic
    coordinates l and b, distance from the Sun rBH. All these values are
    taken from the catalog by Blaha and Humphreys (1989AJ..98...1598B)
    which includes two parts: the catalog of young stars -- members of
    associations or clusters and the catalog of field stars. The last two
    columns exhibit the name of association to which the star is included
    in the partition by Blaha and Humphreys (Assc_BH) and in the partition
    by Mel'nik and Efremov (Assc_ME). The letter F in the last but one
    column means that it is a field star in the catalog by Blaha and
    Humphreys (1989AJ..98...1598B). There are some evidences that the
    distance scale by Blaha and Humphreys (r_BH) should be shorten by
    10-20% (Dambis et al. 2001AstL..27...58D; Mel'nik and Dambis 2009,
    Cat. J/MNRAS/400/518).

File Summary:
--------------------------------------------------------------------------------
 FileName  Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe        80        .   This file
table.dat     75       88   Table of OB-associations
list.dat      67     1392   List of stars -- members of OB associations of
                            the new partition
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
   Bytes Format Units  Label Explanations
--------------------------------------------------------------------------------
   1-  9  A9    ---    AssME Association name
  11- 16  F6.2  deg    GLON  Galactic longitude
  18- 23  F6.2  deg    GLAT  Galactic latitude
  26- 29  F4.2  kpc    D     Distance
  32- 36  F5.1  km/s   RV    ? Average radial velocity of the association
                               relative to the Sun
  40- 43  F4.1  km/s e_RV    ? Dispersion of RV
  46- 47  I2    ---    NRV   Number of stars with measured RV
  50- 52  I3    ---    Ntot  Total number of stars in the association
  57- 61  F5.1  pc     Dl    Size of association along longitude
  65- 68  F4.1  pc     Db    Size of association along latitude
  71      I1    ---    p     [0/1] p=1 for the association with 90% reliability
                                   p=0  for reliability below 90%
  74- 75  I2    ---    Nsg   Number of K and M supergiants
--------------------------------------------------------------------------------

Byte-by-byte Description of file: list.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 14  A14   ---     Name      Designation of the star
  16- 19  A4    ---     SpType    MK spectral type
  21- 25  F5.2  mag     Vmag      Magnitude in V band
  27- 32  F6.2  deg     GLON      Galactic longitude
  34- 39  F6.2  deg     GLAT      Galactic latitude
  42- 45  F4.2  kpc     rBH       Distance from the Sun derived by BH
  48- 55  A8    ---     AssBH     Association of Blaha and Humphreys,
                                   1989AJ..98...1598B
  59- 67  A9    ---     AssME     Association of Mel'nik and Efremov, as
                                   in table.dat file
--------------------------------------------------------------------------------

History:
   * 25-Sep-1997: Original
   * 06-Oct-2011: list.dat added, from Anna Mel'nik, anna(at)sai.msu.ru

================================================================================
(End)                                       Veta Avedisova [INASAN]  10-Jul-1997
