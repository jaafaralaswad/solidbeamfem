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

Balance Relations
====================

The strong form of the equilibrium equations is given by

.. math::

   \operatorname{Div} \, \mathbf{P}
   + \rho_0 \left( \mathbf{b} - \ddot{\mathbf{x}} \right)
   = \mathbf{0}
   \qquad \text{in } \mathcal{B}_0,

where:

- :math:`\mathbf{P}` is the first Piola–Kirchhoff stress
- :math:`\mathbf{b}` is the body force per unit reference volume
- :math:`\rho_0` is the reference density

For static problems, inertial effects are neglected, and

.. math::

   \ddot{\mathbf{x}} = \mathbf{0}.


Material Law: Saint-Venant–Kirchhoff Material
================================================

The Helmholtz free energy is defined in the reference configuration as

.. math::

   \psi_0(\mathbf{C})
   =
   \frac{1}{8}\lambda (I_C - 3)^2
   +
   \frac{1}{4}\mu
   \left(
      I_C^2 - 2 I_C - 2 I\!I_C + 3
   \right),

where

.. math::

   I_C = \operatorname{tr}(\mathbf{C}),
   \qquad
   I\!I_C = \operatorname{tr}(\mathbf{C}^2).

The second Piola–Kirchhoff stress is given by

.. math::

   \mathbf{S}
   =
   \lambda \, \operatorname{tr}(\mathbf{E}) \, \mathbf{I}
   +
   2\mu \, \mathbf{E}.

The contravariant componenets of the fourth-order elasticity tensor are give by

.. math::

   \mathbb{C}^{ijkl}
   =
   \lambda \, G^{ij} G^{kl}
   +
   \mu
   \left(
      G^{ik} G^{jl}
      +
      G^{il} G^{jk}
   \right).


Boundary Conditions
======================

The problem is completed with the following boundary conditions.

Displacement boundary conditions on :math:`\Gamma_u`:

.. math::

   \mathbf{u} = \bar{\mathbf{u}}
   \qquad \text{on } \Gamma_u.

Traction boundary conditions on :math:`\Gamma_t`:

.. math::

   \mathbf{P} \, \mathbf{N}
   =
   \bar{\mathbf{t}}
   \qquad \text{on } \Gamma_t.
