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

Broken Power Law (BPL)
^^^^^^^^^^^^^^^^^^^^^^

The following is a generalization of the equation given by :cite:`Kroupa:2001jy` for :math:`n-1` intervals.

.. math::
    \xi(m)=A
    \begin{cases}
        \ k_{1}m^{-\alpha_{1}} &\quad\text{if }m_{1}\leqslant m< m_{2}\\
        \ k_{2}m^{-\alpha_{2}} &\quad\text{if }m_{2}\leqslant m< m_{3}\\
        \ ...  \\
        \ k_{n-1}m^{-\alpha_{n-1}} &\quad\text{if }m_{n-1}\leqslant m< m_{n}\\
    \end{cases}

where :math:`A` is a normalization constant and :math:`k_{i}` is defined as

.. math::
    k_{1} = m_{2}^{\alpha_{1}} \\
    k_{2} = m_{2}^{\alpha_{2}} \\
    k_{i} = k_{i-1}m_{i}^{\alpha_{i}-\alpha_{i-1}}

This IMF has the benefit of being highly adaptable. It can be utilized to approximate any other IMF.
The BPL has been used to describe both globular cluster IMFs :cite:`Baumgardt_2017` as well as young star cluster IMFs :cite:`Porras_2003` (multi part power law), :cite:`Kroupa_2003` (single power law)

As in the case of Salpeter, random samples are drawn with inverse transform sampling.

The normalized cumulative distribution function (NCFD) can be calculated by integrating :math:`\xi(m)` over the mass interval.

.. math::
    F(m)=\int_{m_{1}}^{m}\xi(m)\,dm

Where A is defined by the normalization constraint:

.. math::
    A * \sum_{i=1}^{n-1} \left ( k_{i}\int_{m_{i}}^{m_{i+1}}m^{-\alpha_{i}}\,dm \right ) = 1

Inverting the NCFD leads to

.. math::
    F^{-1}(y)=
     \begin{cases}
       \ \left [ \frac{(1-\alpha_{1})y}{A*k_{1}} +m_{1}^{1-\alpha_{1}}\right ]^{\frac{1}{1-\alpha_{1}}} &\quad\text{if }0\leqslant y< \frac{A*k_{1}}{1-\alpha_{1}}\left ( m_{2}^{1-\alpha_{1}}-m_{1}^{1-\alpha_{1}} \right )\\
       \ \left \{ \left [y-\frac{A*k_{1}}{1-\alpha_{1}}\left ( m_{2}^{1-\alpha_{1}}-m_{1}^{1-\alpha_{1}} \right )\right ]  \frac{1-\alpha_{2}}{A*k_{2}} + m_{2}^{1-\alpha_{2}} \right \}^{\frac{1}{1-\alpha_{2}}} &\quad\text{if }\frac{A*k_{1}}{1-\alpha_{1}}\left ( m_{2}^{1-\alpha_{1}}-m_{1}^{1-\alpha_{1}} \right )\leqslant y< \sum_{i=1}^{2}\frac{A*k_{i}}{1-\alpha_{i}}\left ( m_{i+1}^{1-\alpha_{i}}-m_{i}^{1-\alpha_{i}} \right )\\
       \ ...  \\
       \ \left \{ \left [ \sum_{i=1}^{n-2} y- \frac{A*k_{i}}{1- \alpha_{i}}\left ( m_{i+1}^{1-\alpha_{i}}-m_{i}^{1-\alpha_{i}} \right )\right ]  \frac{1-\alpha_{n-1}}{A*k_{i}} + m_{n-1}^{1-\alpha_{n-1}} \right \}^{\frac{1}{1-\alpha_{n-1}}} &\quad\text{if }\sum_{i=1}^{n-2}\frac{A*k_{i}}{1-\alpha_{i}}\left ( m_{i+1}^{1-\alpha_{i}}-m_{i}^{1-\alpha_{i}} \right )\leqslant y< \sum_{i=1}^{n-1}\frac{A*k_{i}}{1-\alpha_{i}}\left ( m_{i+1}^{1-\alpha_{i}}-m_{i}^{1-\alpha_{i}} \right )=1\\
     \end{cases}

where y is a random number in range [0,1]

.. doxygenfunction:: InitialConditions::brokenPowerLaw

Spheroid/Bulge - Chabrier (2003)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

The positions of the field stars within the cone of vision are generated in two steps of rejection sampling followed by a transformation.
The cone of vision is defined by the angle of view :math:`\alpha`, the view distance :math:`h` (height of the cone), the view point :math:`vP` (location of the observer) and the focus :math:`F` (a point along the line of sight).

In the first step trial positions are drawn from a uniform distribution within a cuboid containing the cone.
The boundaries of the cuboid are given by

.. math::
    |x|\leq R \\
    |y|\leq R \\
    0\leq z\leq h

where :math:`R=h*\textup{tan}\left ( \frac{\alpha}{2} \right )` is the base radius of the cone.

Those trial positions are rejected in case they are outside the boundaries of the cone.
The conditions for acceptance are:

.. math::
    \sqrt{x^{2}+y^{2}}\leq R \\
    z\geq h*\frac{\sqrt{x^{2}+y^{2}}}{R} \\

This method ensures that the positions are indeed homogeneously distributed which is essential for the second step.

The second step consists of rejection sampling the density distribution.
The test variable is drawn from a uniform distribution ranging from the smallest to the largest possible density within the cone volume.
If this test variable is smaller than the density at the trial position generated in step two, the trial position is accepted and rejected otherwise.

Then the accepted position is transformed via a transformation matrix.
Per this transformation the tip of the cone is displaced from the origin to the view point :math:`vP` and its axis is rotated to align with the line of sight :math:`l`.
Consequently, the transformation consists of both translation and rotation illustrated in the following figure.

.. figure:: Images/cone/cone.svg
    :align: center

    Transformation of the cone of vision

Rotation and translation are both isometric transformations meaning area and volume of the cone are preserved :cite:`Gentle_2007` (p.175).

A unit quaterion :math:`\textbf{q}` is used in order to construct the rotation matrix. With rotation axis :math:`\vec{b}` and angle :math:`\beta` the quaternion is given by

.. math::
    \textbf{q} = \left ( \textup{cos}\left (\frac{\beta}{2}\right ), \vec{b}\textup{ sin}\left ( \frac{\beta}{2} \right )\right )

The rotation axis :math:`\vec{b}` is the normalized cross product of the original (:math:`\vec{z}`) and target (:math:`l`) cone axis

.. math::
    \vec{b}=\frac{\vec{z}\times\vec{l}}{\left \| \vec{z}\times\vec{l} \right \|}

The angle :math:`\beta` between the vectors of interest can be calculated as follows

.. math::
    \beta
    =\textup{atan2}\left ( \textup{tan}\left ( \beta \right ) \right )
    =\textup{atan2}\left ( \frac{\textup{sin}\left ( \beta \right )}{\textup{cos}\left ( \beta \right )} \right )
    =\textup{atan2}\left ( \frac{\left \| \vec{z}\times\vec{l} \right \|}{\vec{z}\cdot \vec{l}} \right )

Next, quarterion is converted to the rotation matrix :cite:`Lee_1991`. Using the homogeneous notation :cite:`Vince_2006` (p. 57) the matrix becomes:

.. math::
    \mathbf{R}=\begin{bmatrix}
    q_{1}^{2}+q_{2}^{2}-q_{3}^{2}-q_{4}^{2} & -2q_{1}q_{4}+2q_{2}q_{3} & 2q_{1}q_{3}+2q_{2}q_{4} & 0\\
    2q_{1}q_{4}+2q_{2}q_{3} & q_{1}^{2}-q_{2}^{2}+q_{3}^{2}-q_{4}^{2} & -2q_{1}q_{2}+2q_{3}q_{4} & 0\\
    -2q_{1}q_{3}+2q_{2}q_{4} & 2q_{1}q_{2}+2q_{3}q_{4} & q_{1}^{2}-q_{2}^{2}-q_{3}^{2}+q_{4}^{2} & 0\\
    0  & 0 & 0 & 1
    \end{bmatrix}

The translation matrix for the translation vector :math:`\vec{t}` reads :cite:`Vince_2006` (p. 66):

.. math::
    \mathbf{T_{translation}}=\begin{bmatrix}
    1 & 0 & 0 & t_{x}\\
    0 & 1 & 0 & t_{y}\\
    0 & 0 & 1 & t_{z}\\
    0  & 0 & 0 & 1
    \end{bmatrix}

The transformation matrix :math:`\mathbf{T}` is the product of :math:`\mathbf{R}` and :math:`\mathbf{T_{translation}}`

.. math::
    \mathbf{T}=\begin{bmatrix}
    q_{1}^{2}+q_{2}^{2}-q_{3}^{2}-q_{4}^{2} & -2q_{1}q_{4}+2q_{2}q_{3} & 2q_{1}q_{3}+2q_{2}q_{4} & t_{x}\\
    2q_{1}q_{4}+2q_{2}q_{3} & q_{1}^{2}-q_{2}^{2}+q_{3}^{2}-q_{4}^{2} & -2q_{1}q_{2}+2q_{3}q_{4} & t_{y}\\
    -2q_{1}q_{3}+2q_{2}q_{4} & 2q_{1}q_{2}+2q_{3}q_{4} & q_{1}^{2}-q_{2}^{2}-q_{3}^{2}+q_{4}^{2} & t_{z}\\
    0  & 0 & 0 & 1
    \end{bmatrix}

.. doxygenfunction:: InitialConditions::sampleDiskPositions(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double coneR, double distance, Matrix *transformationMatrix)

.. doxygenfunction:: InitialConditions::sampleBulgePositions(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double coneR, double distance, Matrix *transformationMatrix)

.. plot:: pyplots/potentialPositions.py

Velocities
----------

Hamiltonian with axisymmetric potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The momenta in cylindrical coordinates are given by

.. math::
    p_{R} = m\dot{R} \\
    p_{\varphi} = mR^2\dot{\varphi} \\
    p_{z} = m\dot{z} \\

therefor the Hamiltonian with an axisymmetric potential reads (:cite:`Binney_2011` p. 278)

.. math::
    H = \frac{1}{2m}\left( p_{R}^2+\frac{p_\varphi^2}{R^2}+p_z^2 \right)+\Phi\left(R,z\right)

using Hamilton's equations gives

.. math::
    \dot{p}_{R} = -\frac{\partial H}{\partial R} = \frac{p_\varphi^2}{R^3}-\frac{\partial \Phi }{\partial R} \\
    \dot{p}_{\varphi} = -\frac{\partial H}{\partial \varphi} = -\frac{\partial \Phi }{\partial \varphi} = 0 \\
    \dot{p}_{z} = -\frac{\partial H}{\partial z} = -\frac{\partial \Phi }{\partial z}

Since :math:`\vec{L} = \vec{p} \times \vec{\dot{p}}` and thus :math:`L_z = R^2\dot{\varphi}`, the second equation above signifies that, in the case of an axisymmetric potential,
the z component of the angular momentum is conserved.

Jeans equations
^^^^^^^^^^^^^^^

Using Einstein notation for :math:`i=1,2,3` the collisionless Boltzmann Equation is given by (:cite:`Binney_2011` p. 277):

.. math::
    \frac{\partial f}{\partial t} + \frac{\partial f}{\partial q_i}\frac{\partial H}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial H}{\partial q_i} = 0

Since the galactic potential (?) is axisymmetric, it is convenient to express this equation in cylindrical coordinates.

.. math::
    \frac{\partial f}{\partial t}
    + p_R\frac{\partial f}{\partial R}
    + \frac{p_\varphi}{R^2}\frac{\partial f}{\partial \varphi}
    + p_z\frac{\partial f}{\partial z}
    - \left(\frac{\partial \Phi}{\partial R}-\frac{p_\varphi^2}{R^3}\right)\frac{\partial f}{\partial p_R}
    - \frac{\partial \Phi}{\partial \varphi}\frac{\partial f}{\partial p_\varphi}
    - \frac{\partial \Phi}{\partial z}\frac{\partial f}{\partial p_z} = 0

It is assumed that the galaxy is statistically in a steady state (:cite:`Binney_2013`) ie :math:`\frac{\partial f}{\partial t}=0`.
Due to this assumption and taking (?) into account (?) simplifies to

.. math::
    p_R\frac{\partial f}{\partial R}
    + \frac{p_\varphi}{R^2}\frac{\partial f}{\partial \varphi}
    + p_z\frac{\partial f}{\partial z}
    - \left(\frac{\partial \Phi}{\partial R}-\frac{p_\varphi^2}{R^3}\right)\frac{\partial f}{\partial p_R}
    - \frac{\partial \Phi}{\partial z}\frac{\partial f}{\partial p_z} = 0

Multiplying equation (?) by :math:`p_R` and integrating over all momenta leads to (todo: derive?)

.. math::
    \frac{\partial \nu \overline{v_R^2}}{\partial R}+\frac{\partial \nu \overline{v_Rv_z}}{\partial z} +
    \nu \left ( \frac{\overline{v_R^2}-\overline{v_\varphi^2}}{R} + \frac{\partial\Phi}{\partial R}\right ) = 0


The Epicyclic Approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Individual stars in the disk are on nearly circular orbits.
Such orbits can be approximated by circular orbits with additional retrograde elliptical orbits around the guiding center.

The derivation of this approximation starts with Hamilton's equations for an axisymmetric potential.

Rearranging equation (?) and using the constant :math:`L_z` gives

.. math::
    H = \frac{1}{2m}\left( p_{R}^2+p_z^2 \right)+\frac{mL_z^2}{2R^2}+\Phi\left(R,z\right)

With the effective potential given by

.. math::
    \Phi_{\textup{eff}}\left(R,z\right) = \frac{mL_z^2}{2R^2}+\Phi\left(R,z\right)

leads to

.. math::
    H_{\textup{eff}} = \frac{1}{2m}\left( p_{R}^2+p_z^2 \right)+\Phi_{\textup{eff}}\left(R,z\right)

Here :math:`\frac{1}{2m}\left( p_{R}^2+p_z^2 \right)` is the kinetic energy in the :math:`\left( R,z \right)` plane or meridional plane.
The angular momentum term in the effective potential is not a real potential energy even though sometimes called centrifugal potential.
It really is the angular kinetic energy. The given definition of :math:`\Phi_{\textup{eff}}` is only valid because :math:`L_z` is conserved.

with this (?) and (?) can be written as

.. math::
    \dot{p}_{R} = -\frac{\partial \Phi_{\textup{eff}} }{\partial R} \\
    \dot{p}_{z} = -\frac{\partial \Phi_{\textup{eff}} }{\partial z}

These equations describe harmonic oscillations in the effective potential.
The minimum of effective potential is the minimum of the real potential energy, together with a contribution from the angular kinetic energy.

.. math::
    \frac{\partial \Phi_{\textup{eff}} }{\partial R} = \frac{\partial \Phi }{\partial R} - \frac{mL_z^2}{2R^3} =0 \\
    \frac{\partial \Phi_{\textup{eff}} }{\partial z} = 0

The first condition states that the attractive force (:math:`-\frac{\partial \Phi_{\textup{eff}} }{\partial R}`) has to balance the “centrifugal force”.
This is the condition for circular orbits with angular momentum :math:`L_z`.
The second condition is clearly satisfied in the equatorial plane (:math:`z=0`).
The coordinates of this guiding center are defined as :math:`(R_g,\varphi_g,z_g)`.

In preparation for a Taylor series expansion about the guiding center :math:`x` is defined as

.. math::
    x \equiv R - R_g

If :math:`R = R_g` then :math:`x=0` and therefor the guiding center is at :math:`(x,z) = (0,0)`.

.. math::
    \Phi_{\textup{eff}} = \Phi_{\textup{eff}}(R_g,0) +
    \frac{\partial \Phi_{\textup{eff}} }{\partial R}\bigg|_{(R_g,0)}x +
    \frac{\partial \Phi_{\textup{eff}} }{\partial z}\bigg|_{(R_g,0)}z +
    \frac{1}{2}\frac{\partial^2 \Phi_{\textup{eff}} }{\partial R^2}\bigg|_{(R_g,0)}x^2 +
    \frac{1}{2}\frac{\partial^2 \Phi_{\textup{eff}} }{\partial z^2}\bigg|_{(R_g,0)}z^2 +
    \frac{1}{2}\frac{\partial^2 \Phi_{\textup{eff}} }{\partial x\partial z}\bigg|_{(R_g,0)}xz +
    \mathcal{O}(xz^2)

The first order terms are zero (since :math:`\Phi_{\textup{eff}}` is minimized at the guiding center) and so is the :math:`xz`, the later due to symmetric about :math:`z=0`.

In the epicyclic approximation :math:`\mathcal{O}(xz^2)` and higher order terms are neglected.

With this approximation (?) and (?) become

.. math::
    \dot{p}_{x} = -\frac{\partial \Phi_{\textup{eff}} }{\partial x} \approx
    -\frac{\partial^2 \Phi_{\textup{eff}} }{\partial R^2}\bigg|_{(R_g,0)}x \equiv
    -\kappa^2x   \\
    \dot{p}_{z} = -\frac{\partial \Phi_{\textup{eff}} }{\partial z} \approx
    -\frac{\partial^2 \Phi_{\textup{eff}} }{\partial z^2}\bigg|_{(R_g,0)}z \equiv
    -\nu^2z   \\

where the epicyclic frequency :math:`\kappa` is the frequency of small radial and the vertical frequency :math:`\nu` the frequency of small vertical oscillations.

with potential :math:`\Phi\left(R,z\right)` (?) can be written as

.. math::
    \kappa^2 = \frac{\partial^2\Phi}{\partial R^2}\bigg|_{(R_g,0)} + \frac{3L_z}{R_g^4}

The circular angular frequency (see eq. (?) with :math:`\Omega = \dot{\varphi}`) is given by

.. math::
    \Omega^2 = \frac{1}{R} \frac{\partial \Phi }{\partial R}\bigg|_{(R_g,0)} = \frac{L_z^2}{R^4}

The derivative of (?) leads to

.. math::
    \frac{\partial^2\Phi}{\partial R^2} = \Omega^2 + R \frac{d\Omega^2}{dR}

Inserting equation (?) and (?) into (?) yields

.. math::
    \kappa^2(R_g) = \left ( R\frac{d\Omega^2}{dR} + 4\Omega^2 \right )\bigg|_{R=R_g}


Disk
^^^^

The velocity distribution of stars in the milky way disk is approximated with the help of Jeans equations as well as relations and constraints based on observations.

For a flat rotation curve the radial velocity dispersion exponentially decreases with increasing radius :cite:`Kruit_1981` (p. 114)

.. math::
    \sigma_{v_{R}} \propto e^{-\frac{R}{h}}

where :math:`h` in the case of the Miyamoto Nagai potential is the radial scale length :math:`a`.

Relation (...) still requires a constant factor :math:`k`, which can be determined by means of the Toomre parameter :math:`Q` at some distance :math:`R_{ref}`

:math:`Q` is the ratio between the actual and minimum velocity dispersion :math:`\sigma_{v_{R,min}}` :cite:`Toomre_1964` (p. 1234)

.. math::
    \sigma_{v_{R,min}} = \frac{3.36G\Sigma}{\kappa} \\
    Q \equiv \frac{\sigma_{v_{R}}}{\sigma_{v_{R,min}}} = \frac{\kappa \sigma_{v_{R}}}{3.36G\Sigma }

where :math:`\kappa` denotes the epicyclic frequency (eq. (?)).

In the solar neighborhood :math:`Q_{\ast} = 2.7 \pm 0.4` and :math:`\sigma_{v_{R}} = (38 \pm 2) \left [ \frac{km}{s} \right ]` :cite:`Binney_2011` (p. 497)

The constant :math:`k` can therefor be approximated via

.. math::
    k \cong Q \sigma_{v_{R,min}}e^{\frac{R}{h}}

Under the approximation of isothermal sheets (introduced in :cite:`Kruit_1981`), the vertical velocity dispersion only depends on the surface density :cite:`Kruit_1988`

.. math::
    \sigma_{v_{z}} = \pi G \Sigma \left ( R \right )z_{0}

with :math:`z_{0}` being the vertical scale length :math:`b` when using the Miyamoto Nagai potential.

The first moments of the collisionless Boltzmann equation (CBE) for cylindrically symmetric systems are given by






Bulge
^^^^^


.. bibliography:: bibtex.bib
