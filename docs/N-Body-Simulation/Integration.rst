===========
Integration
===========

Numerical integration is needed during initialization and simulation.
Various options for time integration have been implemented.
The GNU Scientific Library (:cite:`galassi_2018`, GSL) is used for integrations performed during initialization.

Quadrature, Adaptive, General-purpose (QAG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ported to C in GLS from the Fortran library QUADPACK :cite:`1983Q:as`, this algorithm makes use of adaptive Gauss-Kronrod Quadrature to
estimate the definite integral of a given function.

Quadrature sums are defined as

.. math::
   Q_n[a,b] \equiv \sum_i^nw_if(x_i) \cong \int_a^bw(x)f(x)\textup{d}x




.. bibliography:: bibtex.bib
