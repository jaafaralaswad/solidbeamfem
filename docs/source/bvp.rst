Boundary Value Problem
=========

In order to construct a complete BVP, we need to bring together kinematics, balance relations and material law, in addition to the boundary conditions. We formulate the problem using convective coordinates for better structural mechanical interpretation and later implementation flexibility.

Kinematics
=============

We distinguish between the reference configuration :math:`\mathcal{B}_0`
with material coordinates :math:`\mathbf{X}`, and the current configuration
:math:`\mathcal{B}` with spatial coordinates :math:`\mathbf{x}`.
The displacement is given by

.. math::

   \mathbf{u} = \mathbf{x} - \mathbf{X}.

The deformation gradient is computed via the chain rule:

.. math::

   \mathbf{F}
   = \frac{\partial \mathbf{x}}{\partial \mathbf{X}}
   = \frac{\partial \mathbf{x}}{\partial \xi^i}
     \otimes
     \frac{\partial \xi^i}{\partial \mathbf{X}}
   = \mathbf{g}_i \otimes \mathbf{G}^i,

where:

- :math:`\xi^i \in \{\xi, \eta, \zeta\}` are the convective coordinates
- :math:`\mathbf{g}_i` are the covariant base vectors in the current configuration
- :math:`\mathbf{G}^i` are the contravariant base vectors in the reference configuration

The Green–Lagrange strain tensor is given by

.. math::

   \mathbf{E}
   = \frac{1}{2} \left( \mathbf{C} - \mathbf{I} \right),
   \qquad
   \mathbf{C} = \mathbf{F}^\mathsf{T} \mathbf{F}.

The contravariant components of the strain tensor are written as

.. math::

   \mathbf{E}
   = E^{ij} \, \mathbf{G}_i \otimes \mathbf{G}_j,
   \qquad
   E^{ij} = \frac{1}{2} \left( g^{ij} - G^{ij} \right).

The metric tensors

.. math::

   g^{ij} = \mathbf{g}_i \cdot \mathbf{g}_j,
   \qquad
   G^{ij} = \mathbf{G}_i \cdot \mathbf{G}_j,

allow transformation between covariant and contravariant components.

In Voigt notation, we exploit the symmetry of the strain tensor to write

.. math::

   \hat{\mathbf{E}}
   =
   \begin{bmatrix}
   E_{\xi\xi} &
   E_{\eta\eta} &
   E_{\zeta\zeta} &
   2E_{\xi\eta} &
   2E_{\eta\zeta} &
   2E_{\xi\zeta}
   \end{bmatrix}^{\mathsf{T}} .
