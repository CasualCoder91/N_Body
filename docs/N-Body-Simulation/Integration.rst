===========
Integration
===========

Numerical integration is needed during initialization and simulation.
Various options for time integration have been implemented.
The GNU Scientific Library (:cite:`galassi_2018`, GSL) is used for integrations performed during initialization.

Quadrature, Adaptive, General-purpose (QAG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ported to C in GLS from the Fortran library QUADPACK :cite:`1983Q:as`, this algorithm makes use of adaptive Gauss-Kronrod quadrature to
estimate the definite integral of a given function.

Quadrature sums are defined as

.. math::
   Q_n[a,b] \equiv \sum_i^nw_if(x_i) \cong \int_a^bw(x)f(x)\textup{d}x

where :math:`w_i` are weights, :math:`x_i` nodes, :math:`w(x)` a weight function.
The highest possible degree of precision is 2n-1. With this maximum precision :math:`Q_n` is exact for polynomials of degree smaller or equal to 2n-1.

Using classical Gaussian quadrature formulae, error estimation, by increasing n to n+1, requires n+1 evaluations of :math:`f(x)` in addition to n evaluations from
calculation of the original sum, since the respective nodes have no common points. By doing so, the degree of precision is only increased from 2n-1 to 2n+1.
Therefor, the error estimation obtained by subtracting the two sums could be unreliable.

Adding n+1 points to the Gauss-Legendre formula - here :math:`w(x)=1` and the nodes are zeros of the Legendre polynomial -
Kronrod introduced the option of increasing the precision to 3n+1, again requiring n+1 additional evaluations of :math:`f(x)` (:cite:`Notaris_2016`).

.. math::
   Q_n^K[a,b] \equiv \sum_i^nw_if(x_i) + \sum_j^{n+1}w_j^*f(x_j^*) \cong \int_{-1}^1f(x)\textup{d}x

QAG makes use of this option, subdividing the interval with the largest local absolute error estimate in each step.
This division is repeated until either the absolute or relative global error estimate are smaller than required by the caller.



.. bibliography:: bibtex.bib
