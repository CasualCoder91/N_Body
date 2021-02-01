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

Velocity Verlet Algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^

For cluster members the acceleration is a combination of the force resulting from the presence of all other cluster stars (see Barens Huts Algorithm)
and from the milky way potential. The acceleration of field stars solely comes from the milky way potential.
In each time step both velocity and acceleration of each star is evaluated.

Since the velocity, :math:`v(t)` changes over time, it's value at the midpoint between the current (:math:`t_n`) and the next timestep :math:`t_{n+1}` is intuitively a better approximation than
:math:`v(t_n)` or :math:`v(t_{n+1})`. The same holds true for the acceleration. This leads to :cite:`feynman_1965`

.. math::
   x_{n+1} = x_{n} + hv_{n+0.5}\\
   v_{n+1.5} = v_{n+0.5} + \frac{h}{m}F(x_{n+1})

which is the Leapfrog algorithm. If one needs :math:`x` and :math:`v` at the same time, (?) can be split into two halve steps.

.. math::
   v_{n+0.5} = v_{n} + \frac{h}{2m}F(x_{n})
   x_{n+1} = x_{n} + hv_{n+0.5}\\
   v_{n+1} = v_{n+0.5} + \frac{h}{2m}F(x_{n+1})

:math:`F(x)` does not have to be calculated twice, because :math:`F(x_{n+1})` can be used as :math:`F(x_{n})` in the next timestep.

symplectic

Barnes-Hut Algorithm (BH)
^^^^^^^^^^^^^^^^^^^^^^^^^

When determining the gravitational force acting on a star which belonging to a cluster, the positions of all other stars in that cluster
have to be taken into account.

.. math::
   m_i\vec{x}_i = G\sum_{j=1,j\not\equiv i}^{N}\frac{m_im_j\left ( \vec{x}_j-\vec{x}_i \right )}{\left | \vec{x}_j-\vec{x}_i \right |^3}

Calculating this force for all stars requires :math:`O(n^2)` operations.
However, the simulated clusters consist of :math:`\sim 10^4 - 10^5` stars.
While the resulting amount of calculations is possible, it is not feasible for a typical desktop computer.
Therefor, the Barnes-Hut algorithm (BH) has been implemented which is of order :math:`O(n\log(n))`.

The gist of the BH is to approximate a set of stars by their total mass and center of mass (com) if the distance between them
and the star, for which the force is to be calculated, is large enough.

The total mass and com of a set of :math:`m` stars is

.. math::
   m_{com} = \sum_{i=1}^mm_i \\
   \vec{x}_{com} = \frac{1}{m_{com}}\sum_{i=1}^mm_i\vec{x}_i

All cluster stars are stored in an octree.
An octree is a data structure where each node in the tree has up to eight child nodes.
These nodes split the space represented by their parent node into eight cubes.
External nodes are nodes without any children. Each external node contains at most one star.
Internal nodes have at least one child. They represent stars stored in their child nodes by storing their total mass and com.
The root node contains the whole space occupied by the cluster. Each node stores the following information: total mass, amount and center of mass of stars
contained within the cube, two points defining the volume of the cube, one point at the center of the cube, since cpu time is more valuable than ram,
as well as links (pointers) to each child node and to the parent node.
If a child pointer is null, it does not exist jet.

Stars are added recursively starting at the root node. If the current node is already a internal node, the star is passed
to the appropriate child. Mass and com of the internal node are updated.
The appropriate child is determined by comparing the position of the star with the center of the node.
If the considered node is a external node but already contains a star,
both the newly added and already present star are passed down to the appropriate child or children.
Consequently the current node becomes a internal node.
Since both stars can lie in the same octant, this can lead to additional recursions until the stars are assigned to different child nodes.
If the current node is external and does not yet contain a star, the star is added to the node and the recursion ends.

When calculating the gravitational force on a star, the octree is travelled through recursively beginning with the root node.
In case the distance :math:`d` between the star and a node is sufficiently large, the stars within that node are approximated by the mass and com of that node,
otherwise all child nodes within the current node are considered. Whenever the distance criterion is met,
the acceleration vector stemming from the force is calculated, added to the overall acceleration of the star and the recursion for the considered branch stops.

Whether or not :math:`d` is big enough, is determined by the quotient :math:`\theta`.

.. math::
   \theta = \frac{s}{d} < \theta_{max}

with :math:`s` the side length of the cube and :math:`\theta_{max}` a set threshold value.
In the special case :math:`\theta_{max}=0`, BH becomes a direct-sum algorithm. :math:`\theta_{max}=0.5` is a commonly chosen value.

doto: explain smoothing

.. bibliography:: bibtex.bib
