==============
Galactic Potential
==============

https://iopscience.iop.org/article/10.1088/0004-637X/794/1/4
https://academic.oup.com/mnras/article/487/3/4025/5511782
(halo irrelevant(?)): https://arxiv.org/abs/1907.13132)

The potential consists of three parts: black hole, buldge and disc. The dark matter halo is not being considered since it is irrelevant for this simulation.

The black hole is represented by a Keplerian potential.

..  math::
    \Phi_{bh}\left ( r \right ) = -\frac{G*M_{bh}}{r}

the disk can be modeled via a Miyamoto Nagai potential

.. math::
    \Phi_{disk}\left ( r \right ) = -\frac{G*M_{disk}}{\sqrt{R^{2+\left ( a+\sqrt{z^{2}+b^{2}} \right )^{2}}}}

and for the buldge the Hernquist potential is used

.. math::
    \Phi_{buldge}\left ( r \right ) = -\frac{G*M_{buldge}}{\left ( r+a \right )}
