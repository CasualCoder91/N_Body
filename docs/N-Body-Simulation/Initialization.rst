==============
Initialization
==============

Initial Mass Function (IMF)
---------------------------

Salpeter (1955)
^^^^^^^^^^^^^^^

.. doxygenfunction:: initialMassSalpeter

.. plot:: pyplots/initialConditionsMassSalpeter.py

Background

When it comes to sampling an IMF, one possible approach is called inverse transform sampling. Here one must integrate the IMF over the mass, yielding a cumulative probability function (cpf), and normalize it (ncpf).
Then one has to take the inverse of the ncpf. Since Salpeter is a power-law distribution function all this can be done analytically:

 https://www.usna.edu/Users/cs/crabbe/SI475/current/particleFilter/particleFilter.pdf
 https://local.strw.leidenuniv.nl/events/phdtheses/haas/05.pdf

.. math::
    p(m)=\frac{dN}{dm}=A*m^{-\alpha }

    P(m)=\int_{m_{min}}^{m}A*m^{-\alpha } dm = \frac{A}{-\alpha +1}\left ( m^{-\alpha +1} -m_{min}^{-\alpha +1}\right )

A is defined by normalization:

.. math::
    P({m_{max}})\equiv 1\rightarrow A=\frac{-\alpha +1}{m_{max}^{-\alpha +1} -m_{min}^{-\alpha +1} }

Inserting this into P(m) yields:

.. math::
    P(m)=\frac{m^{-\alpha +1} -m_{min}^{-\alpha +1}}{m_{max}^{-\alpha +1} -m_{min}^{-\alpha +1}}

Inverting this function and some simplifications lead to:

.. math::
    m_{rand} = m_{min}*\left [ 1+x*\left ( \left ( \frac{m_{max}}{m_{min}} -1\right )^{-\alpha +1} \right ) \right ]^{\frac{1}{-\alpha +1}}

where x is a random number in range [0,1]


Present Day Mass Function (PDMF)
--------------------------------

Disk Stellar Mass Function
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: InitialConditions::massDisk

Stars belonging to the disk are given a mass by rejection sampling the PDMF as given by :cite:`2003PASP..115..763C`.

For :math:`m<1` the PDMF reads

.. math::
    \xi\left(\mathrm{log}(m)\right) = \frac{\mathrm{dN}}{\mathrm{dlog}(m))} = A*\mathrm{exp}[ \frac{-( \mathrm{log}(m) -\mathrm{log}( m_{c} ))^{2} }{2 \sigma^{2}}]

or equivalently (this version is sampled)

.. math::
    \xi\left(m\right) = \frac{\mathrm{dN}}{\mathrm{dm}} = \frac{A}{m\mathrm{ln}(10)}*\mathrm{exp}[ \frac{-( \mathrm{log}(m) -\mathrm{log}( m_{c} ))^{2} }{2 \sigma^{2}}]

For :math:`m>1` the PDMF has the form

.. math::
    \xi\left(\mathrm{log}(m)\right) = \frac{\mathrm{dN}}{\mathrm{dlog}(m))} = A m^{-x}

or depending on :math:`m` rather than :math:`\mathrm{log}(m)`

.. math::
    \xi\left(m\right) = \frac{\mathrm{dN}}{\mathrm{dm}} = \frac{A}{m\mathrm{ln}(10)} m^{-x}

.. plot:: pyplots/initialConditionsMassDisk.py

Positions
---------

Miyamoto & Nagai
^^^^^^^^^^^^^^^^

mass distribution (from Poisson equation)(Source: http://articles.adsabs.harvard.edu/pdf/1975PASJ...27..533M)

.. math::
    \rho \left ( R,z \right )=\frac{b_{disk}^{2}M_{disk}}{4\pi }\frac{a_{disk}R^{2}+\left ( a+3\left (z^{2}+b_{disk}^{2}  \right )^{\frac{1}{2}} \right )\left ( a+\left ( z^{2}+b_{disk}^{2} \right )^{\frac{1}{2}} \right )^{2}}{\left ( R^{2}+\left ( a+\left ( z^{2}+b_{disk}^{2} \right )^{\frac{1}{2}} \right )^{2} \right )^{\frac{5}{2}}\left (  z^{2}+b_{disk}^{2}\right )^{\frac{3}{2}}}

**Radial distribution**

normalized cumulative distribution function

.. math::
    E^{*}\left ( R \right )=\frac{1}{\Sigma _{0}}\int_{0}^{R}\Sigma ^{*}\left ( R' \right )dR'

with probability distribution function in cylindrical coordinates for the radial component:

.. math::
    \Sigma^{*}\left(R\right )=R\int_{-\infty}^{\infty}\rho\left ( R,z \right )dz

and normalization constant

.. math::
    \Sigma_{0}=\int_{0}^{\infty}\Sigma^{*}\left(R'\right )dR'


.. bibliography:: bibtex.bib
