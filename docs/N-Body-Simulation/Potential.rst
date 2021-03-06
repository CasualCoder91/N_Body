==================
Galactic Potential
==================

.. https://iopscience.iop.org/article/10.1088/0004-637X/794/1/4
   https://academic.oup.com/mnras/article/487/3/4025/5511782
   (MWPotential2014) https://arxiv.org/pdf/1412.3451.pdf
   (more on dark halo) https://www.researchgate.net/publication/316334286_Mass_Distribution_and_Gravitational_Potential_of_the_Milky_Way
   (another galactic potential) http://www.astro.utu.fi/~cflynn/galdyn/lecture4.html
   (halo from here if needed: https://iopscience.iop.org/article/10.1088/0004-637X/794/1/4)
   (Measurements of circular velocity) https://arxiv.org/pdf/1810.09466.pdf
   (circular velocity, using derivation) https://www.researchgate.net/publication/316680117_GalRotpy_an_educational_tool_to_understand_and_parametrize_the_rotation_curve_and_gravitational_potential_of_disk-like_galaxies

The potential consists of four parts: black hole, bulge, disc and the dark matter halo.

The **black hole** is represented by a Keplerian potential:

..  math::
    \Phi_{bh}\left ( r \right ) = -\frac{G*M_{bh}}{r}

with :math:`r` being the spherical radius.

The **disk** can be modeled via a Miyamoto Nagai potential :cite:`1975PASJ..27..533M`

.. math::
    \Phi_{disk}\left ( R,z \right ) = -\frac{G*M_{disk}}{\sqrt{R^{2+\left ( a_{disk}+\sqrt{z^{2}+b_{disk}^{2}} \right )^{2}}}}

    \rho_{disk}(R,z)=\frac{b_{disk}^{2}M_{disk}}{4\pi}\frac{a_{disk}R^{2}+\left[a_{disk}+3(z^{2}+b_{disk}^{2})^{\frac{1}{2}}\right]\left[a_{disk}+(z^{2}+b_{disk}^{2})^{\frac{1}{2}}\right]^{2}}{\left \{ R^{2}+\left[a_{disk}+(z^{2}+b_{disk}^{2})^{\frac{1}{2}}\right]^{2} \right \}^{\frac{5}{2}}\left(z^{2}+b_{disk}^{2}\right)^{\frac{3}{2}}}

with :math:`R` the cylindrical radius and with :math:`z` the distance to the reference plane.

For the **bulge** the Hernquist potential :cite:`1990ApJ..356..359H` is used

.. math::
    \Phi_{bulge}\left ( r \right ) = -\frac{G*M_{bulge}}{\left ( r+a_{bulge} \right )}

    \rho_{bulge}(r)= \frac{M_{bulge}}{2\pi }\frac{a_{bulge}}{r}\frac{1}{\left ( r+a_{bulge} \right )^{3}}

:math:`a` is the scale-length of the spheroid potential

and NFW potential is used for the **dark matter halo** potential

.. math::
    \Phi_{halo}\left ( r \right ) = \frac{-4\pi G\rho _{s}r_{s}^{3}\ln\left ( 1+\frac{r}{r_{s}}\right )}{r}

where :math:`\rho _{s}` is the characteristic density and :math:`r_{s}` the scale length.

Parameters for bulge and disk taken from :cite:`Price_Whelan_2014` and the radius :math:`r_{s}` for the halo from :cite:`Bovy_2015`:

.. (rs: https://arxiv.org/pdf/1412.3451.pdf)
   (others: https://iopscience.iop.org/article/10.1088/0004-637X/714/1/229)

.. csv-table:: Parameters
   :header: "Parameter", "Value", "Unit"
   :widths: 20, 20, 10

   :math:`M_{bh}`, :math:`4*10^{6}`, :math:`M_\odot`
   :math:`M_{disk}`, :math:`10^{11}`, :math:`M_\odot`
   :math:`a_{disk}`, :math:`6.5`, :math:`kpc`
   :math:`b_{disk}`, :math:`0.26`, :math:`kpc`
   :math:`M_{bulge}`, :math:`3.4*10^{10}`, :math:`M_\odot`
   :math:`a_{bulge}`, :math:`0.70`, :math:`kpc`
   :math:`r_{s}`, :math:`16`, :math:`kpc`

:math:`\rho _{s}` can be determined by imposing

.. math::
   v_{c}\left ( R_{0},z=0 \right ) = 220\left [ \frac{km}{s} \right ]

   R_{0} = 8 \left [ kpc \right ]

The circular velocity :math:`v_{c}` is defined via

.. math::
   v_{c}\left ( R,z \right ) = \sqrt{R\frac{\partial \Phi \left (R,z  \right )}{\partial R}   }

with in the presented case total potential being

.. math::
   \Phi \left ( R,z \right ) = \Phi_{bh} \left ( R,z \right )+\Phi_{disk} \left ( R,z \right )+\Phi_{bulge} \left ( R,z \right )+\Phi_{halo} \left ( R,z \right )

therefor

.. math::
   v_{c}^{2} = v_{c,bh}^{2}+v_{c,disk}^{2}+v_{c,bulge}^{2}+v_{c,halo}^{2}

   v_{c,bh}^{2} = \frac{G M_{bh} R^2}{\left(R^2+z^2\right)^{3/2}}

   v_{c,disk}^{2} = \frac{GM_{disk}R^2}{\left(\left(a_{disk}+\sqrt{b_{disk}^2+z^2}\right)^2+R^2\right)^{3/2}}

   v_{c,bulge}^{2} = \frac{G M_{bulge} R^2}{\sqrt{R^2+z^2} \left(a_{bulge}+\sqrt{R^2+z^2}\right)^2}

   v_{c,halo}^{2} = \frac{4 \pi G \rho _{s} R^2 r_{s}^3 \log \left(\frac{\sqrt{R^2+z^2}}{r_{s}}+1\right)}{\left(R^2+z^2\right)^{3/2}}-\frac{4\pi G\rho_{s} R^2 {r_{s}}^2}{\left(R^2+z^2\right) \left(\frac{\sqrt{R^2+z^2}}{r_{s}}+1\right)}

.. math::
   v_{c}^{2}=R\left ( \frac{G M_{bh} R}{\left(R^2+z^2\right)^{3/2}}+
   \frac{GM_{disk}R}{\left(\left(a_{disk}+\sqrt{b_{disk}^2+z^2}\right)^2+R^2\right)^{3/2}}+
   \frac{G M_{bulge} R}{\sqrt{R^2+z^2} \left(a_{bulge}+\sqrt{R^2+z^2}\right)^2}+
   \frac{4 \pi G \rho _{s} R r_{s}^3 \log \left(\frac{\sqrt{R^2+z^2}}{r_{s}}+1\right)}{\left(R^2+z^2\right)^{3/2}}-\frac{4\pi G\rho_{s} R {r_{s}}^2}{\left(R^2+z^2\right) \left(\frac{\sqrt{R^2+z^2}}{r_{s}}+1\right)} \right )

plugging in all the parameters (and :math:`G\approx 4.302*10^{-6}\left [ \frac{kpc}{M_\odot}\frac{km^{2}}{s^{2}} \right ]`) results in :math:`\rho_{s}\approx 4.5*10^{6} \left [ \frac{M_\odot}{kpc^{3}} \right ]`

Circular Velocity
-----------------

.. plot:: pyplots/potentialCircularVelocity.py

Angular Velocity
----------------

.. math::
    \Omega ^{2}\left ( R \right ) = \frac{1}{R}\frac{\partial \Phi \left ( R,0 \right )}{\partial r}

    \Omega ^{2}\left ( R \right ) = \frac{G}{R} \left\{-\frac{M_{bulge}}{(a_{bulge}+R)^2}+\frac{2 M_{disk} R^3}{\left[\left(a_{disk}+b_{disk}\right)^2+R^4\right]^{1.5}}+\frac{M_{bh}}{R^2}-\frac{4 \pi  p_{s} r_{s}^3}{R^2+R r_{s}}+\frac{4 \pi  p_{s} r_{s}^3 \ln \left(\frac{R+r_{s}}{r_{s}}\right)}{R^2}\right\}

Mass Distribution
-----------------

The mass inside a volume is calculated by numerical integration of the density.
GSL implementation of Monte Carlo Integration is used. For further details refer to the GSL documentation_.

.. _documentation: https://www.gnu.org/software/gsl/doc/html/montecarlo.html

Example at z=1pc

.. plot:: pyplots/massDistribution.py

Surface Mass Density (SMD)
--------------------------

The SMD is defined by

.. math::
    \Sigma \left ( R \right )=2\int_{0}^{\infty}\rho \left ( r \right )\mathrm{d}z

.. doxygenfunction:: Potential::surfaceDensityDisk

.. doxygenfunction:: Potential::surfaceDensityBulge

.. plot:: pyplots/potentialSurfaceDensity.py

"Junk"

Halo potential given by (https://iopscience.iop.org/article/10.1088/0004-637X/714/1/229)

.. math::
    \Phi_{halo}\left ( x,y,z \right ) = v_{halo}^{2}*\ln\left ( C_{1}*x^{2}+C_{2}*y^{2}+C_{3}*x*y +\left (\frac{z}{q_{z}}  \right )^{2}+r_{halo}^{2}\right )

    \Phi_{halo}\left ( r \right ) = \frac{1}{2}v_{halo}^{2}\ln\left ( r^{2}+r_{halo}^{2}\right )

.. bibliography:: bibtex.bib
