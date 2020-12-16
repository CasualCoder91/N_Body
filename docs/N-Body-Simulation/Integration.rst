===========
Integration
===========

Numerical integration is needed during initialization and simulation.
Various options for time integration have been implemented.
The GNU Scientific Library (:cite:`galassi_2018`, GSL) is used for integrations performed during initialization.
The relevant functions are all based on QAG or QAGI, which have ben ported from the Fortran library QUADPACK :cite:`1983Q:as` to C in GLS.
The decision trees given on page 79-80 in :cite:`1983Q:as` help with the decision on when and how to use the respective methods.

Quadrature, Adaptive, General-purpose (QAG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This algorithm makes use of adaptive Gauss-Kronrod quadrature to estimate the definite integral of a given function.

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

QAG makes use of this option, bisecting the interval with the largest local absolute error estimate in each step.
This division is repeated until either the absolute or relative global error estimate are smaller than required by the caller.

Quadrature, Adaptive, General-purpose, Infinite interval (QAGI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In case of a semi infinite interval :math:`(a,\pm \infty)`, the integration variable is transformed

.. math::
   x = a\pm\frac{1-t}{t}

leading to

.. math::
   \int_a^{\pm \infty}f(x)\textup{d}x = \pm \int_0^1 f \left(a\pm\frac{1-t}{t}\right)t^{-2}\textup{d}t

For a infinite interval

.. math::
   \int_{-\infty}^{\infty}f(x)\textup{d}x =
   \int_0^\infty f(x)+f(-x)\textup{d}x =
   \int_0^1  \left ( f \left(\frac{1-t}{t}\right) + f \left(\frac{t-1}{t}\right) \right )  t^{-2}\textup{d}t

After the transformation QAGS with the 15-point Kronrod rule is used.
QAGS, in addition to the adaptive bisection (see QAG), makes use of the Wynn :math:`\epsilon`-algorithm to accelerate the convergence.

Leapfrog
^^^^^^^^

For cluster members the acceleration is a combination of the force resulting from the presence of all other cluster stars (see Barens Huts Algorithm)
and from the milky way potential. The acceleration of field stars solely comes from the milky way potential.
In each time step both velocity and acceleration of each star is evaluated.




.. bibliography:: bibtex.bib
