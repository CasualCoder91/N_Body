==============
Initialization
==============

Initial Mass Function (IMF)
---------------------------

Salpeter (1955)
^^^^^^^^^^^^^^^

.. doxygenfunction:: InitialConditions::initialMassSalpeter

.. plot:: pyplots/initialConditionsMassSalpeter.py

Background

When it comes to sampling an IMF, one possible approach is called inverse transform sampling. Here one must integrate the IMF over the mass, yielding a cumulative probability function (cpf), and normalize it (ncpf).
Then one has to take the inverse of the ncpf. Since Salpeter is a power-law distribution function all this can be done analytically:

 https://www.usna.edu/Users/cs/crabbe/SI475/current/particleFilter/particleFilter.pdf
 https://local.strw.leidenuniv.nl/events/phdtheses/haas/05.pdf

.. math::
    p(m)=\frac{dN}{dm}=A*m^{-\alpha }
    :label: salpeter

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

Spheroid/Bulge - Chabrier (2003)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: InitialConditions::bulgeIMF

.. plot:: pyplots/initialConditionsMassBulge.py

Per rejection sampling the following function, given by :cite:`2003PASP..115..763C`, the mass of stars, which belong to the bulge, is determined.

For :math:`m<0.7` the log-normal distribution equation :eq:`lognormal` is used. Parameters are :math:`A=3.6*10^{-4}`, :math:`m_{c}=0.22` and :math:`\sigma=0.33`.
For :math:`m>0.7` a Salpeter slope :eq:`salpeter` with parameters :math:`A=7.1*10^{-5}` and :math:`x=1.3` is chosen.


Present Day Mass Function (PDMF)
--------------------------------

.. Bulge: (m>1) http://adsabs.harvard.edu/full/1999A%26A...348..457M (m<1) https://hubblesite.org/uploads/science_paper/file_attachment/200/pdf.pdf

Disk Stellar Mass Function
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: InitialConditions::diskIMF

.. plot:: pyplots/initialConditionsMassDisk.py

Stars belonging to the disk are given a mass by rejection sampling the PDMF as given by :cite:`2003PASP..115..763C`.

For :math:`m<1` the PDMF reads

.. math::
    \xi\left(\mathrm{log}(m)\right) = \frac{\mathrm{dN}}{\mathrm{dlog}(m))} = A*\mathrm{exp}[ \frac{-( \mathrm{log}(m) -\mathrm{log}( m_{c} ))^{2} }{2 \sigma^{2}}]
    :label: lognormal

or equivalently (this version is sampled)

.. math::
    \xi\left(m\right) = \frac{\mathrm{dN}}{\mathrm{dm}} = \frac{A}{m\mathrm{ln}(10)}*\mathrm{exp}[ \frac{-( \mathrm{log}(m) -\mathrm{log}( m_{c} ))^{2} }{2 \sigma^{2}}]

For :math:`m>1` the PDMF has the form

.. math::
    \xi\left(\mathrm{log}(m)\right) = \frac{\mathrm{dN}}{\mathrm{dlog}(m))} = A m^{-x}

or depending on :math:`m` rather than :math:`\mathrm{log}(m)`

.. math::
    \xi\left(m\right) = \frac{\mathrm{dN}}{\mathrm{dm}} = \frac{A}{m\mathrm{ln}(10)} m^{-x}


Positions
---------

Individual Star positions within some volume are sampled directly from the density via rejection sampling.

.. doxygenfunction:: sampleDiskPositions

.. doxygenfunction:: sampleBulgePositions

.. plot:: pyplots/potentialPositions.py


.. bibliography:: bibtex.bib
