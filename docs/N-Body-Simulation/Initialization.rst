==============
Initialization
==============

Initial Mass Function (IMF)
---------------------------

Salpeter (1955)
^^^^^^^^^^^^^^^

.. doxygenfunction:: initialMassSalpeter

Background

In order to sample an IMF one must integrate it over the mass, yielding a cumulative probability function (cpf), and normalize it (ncpf).
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
