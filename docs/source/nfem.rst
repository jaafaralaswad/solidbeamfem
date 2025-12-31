Nonlinear Finite Element Procedure
=========

We work in convective curvilinear coordinates, and while we follow the
general framework presented in Nonlinear Finite Element Methods by
Wriggers (2008), several steps are modified to suit this coordinate system.
Only the key results and necessary adjustments are summarized here.
We assume the reader is familiar with the general procedure and provide
only brief reminders.

The weak form of the problem, neglecting inertial effects and body
forces, is given by

.. math::

   \int_{\Omega_0} \mathbf{S} : \delta \mathbf{E} \, \mathrm{d}V
   -
   \int_{\partial \Omega_0} \mathbf{t}_0 \cdot \delta \mathbf{x} \, \mathrm{d}A
   = 0,

where:

- :math:`\mathbf{S}` is the second Piola–Kirchhoff stress tensor
- :math:`\mathbf{t}_0` is the surface traction in the reference configuration

Equivalently, one may write

.. math::

   g(\mathbf{x}, \delta \mathbf{x})
   =
   g^{\mathrm{int}}(\mathbf{x}, \delta \mathbf{x})
   -
   g^{\mathrm{ext}}(\mathbf{x}, \delta \mathbf{x})
   = 0,

with

.. math::

   g^{\mathrm{int}}(\mathbf{x}, \delta \mathbf{x})
   :=
   \int_{\Omega_0} \mathbf{S} : \delta \mathbf{E} \, \mathrm{d}V,

.. math::

   g^{\mathrm{ext}}(\mathbf{x}, \delta \mathbf{x})
   :=
   \int_{\partial \Omega_0} \mathbf{t}_0 \cdot \delta \mathbf{x} \, \mathrm{d}A.

The linearization of the internal virtual work is given by

.. math::

   \Delta g^{\mathrm{int}}
   =
   \int_{\Omega_0}
   \Delta \mathbf{E} : \mathbb{C}(\mathbf{E}) : \delta \mathbf{E} \, \mathrm{d}V
   +
   \int_{\Omega_0}
   \mathbf{S} : \Delta \delta \mathbf{E} \, \mathrm{d}V,

where:

- :math:`\int_{\Omega_0} \Delta \mathbf{E} : \mathbb{C}(\mathbf{E}) : \delta \mathbf{E} \, \mathrm{d}V`
  is the material contribution
- :math:`\int_{\Omega_0} \mathbf{S} : \Delta \delta \mathbf{E} \, \mathrm{d}V`
  is the geometric contribution

The linearization of the external virtual work vanishes when the applied
load is conservative. However, for non-conservative loads (e.g. follower
loads), the linearization is generally non-zero and must be taken into account.
This leads to the following consequences:

- The load must be updated during every Newton–Raphson iteration
  to remain consistent.
- An additional (third) contribution to the tangent matrix must be added,
  but only for elements that are directly loaded. As a result, the tangent
  matrix is no longer symmetric.

For details on the linearization of the external virtual work, see
Wriggers (2008), Section 4.2.5. The procedure followed here is exactly as
outlined therein.

In this work, a bending moment acting at the tip of the beam is modeled
as a follower load, defined as

.. math::

   \hat{p}
   =
   -
   \left(
      \frac{12 M}{I}
   \right)
   \zeta,
   \qquad
   \zeta \in \left[ -\frac{h}{2}, \frac{h}{2} \right].

The load is applied as a first Piola–Kirchhoff stress tensor, meaning that
it is defined over the undeformed area. The magnitude of the load remains
constant, while its direction changes with deformation.

.. figure:: images/bending_traction.jpg
   :width: 65%
   :align: center
   :alt: Bending applied as follower load


Within the isoparametric concept, both the geometry and the displacement
field are approximated using the same shape functions
:math:`N_I(\xi,\eta,\zeta)`:

.. math::

   \mathbf{u}(\xi,\eta,\zeta)
   \approx
   \sum_{I=1}^{n_e} N_I(\xi,\eta,\zeta)\, \mathbf{u}_I.

The mapping between the parametric coordinates
:math:`(\xi,\eta,\zeta)` and the physical coordinates
:math:`(X_1,X_2,X_3)` is defined through the Jacobian matrix:

.. math::

   \begin{bmatrix}
   \dfrac{\partial (\cdot)}{\partial \xi} \\
   \dfrac{\partial (\cdot)}{\partial \eta} \\
   \dfrac{\partial (\cdot)}{\partial \zeta}
   \end{bmatrix}
   =
   \mathbf{J}
   \begin{bmatrix}
   \dfrac{\partial (\cdot)}{\partial X_1} \\
   \dfrac{\partial (\cdot)}{\partial X_2} \\
   \dfrac{\partial (\cdot)}{\partial X_3}
   \end{bmatrix},

with

.. math::

   \mathbf{J}
   =
   \begin{bmatrix}
   \dfrac{\partial X_1}{\partial \xi} &
   \dfrac{\partial X_2}{\partial \xi} &
   \dfrac{\partial X_3}{\partial \xi} \\
   \dfrac{\partial X_1}{\partial \eta} &
   \dfrac{\partial X_2}{\partial \eta} &
   \dfrac{\partial X_3}{\partial \eta} \\
   \dfrac{\partial X_1}{\partial \zeta} &
   \dfrac{\partial X_2}{\partial \zeta} &
   \dfrac{\partial X_3}{\partial \zeta}
   \end{bmatrix}.

Inverted, the gradient of the shape functions in physical coordinates is

.. math::

   \begin{bmatrix}
   \dfrac{\partial N_I}{\partial X_1} \\
   \dfrac{\partial N_I}{\partial X_2} \\
   \dfrac{\partial N_I}{\partial X_3}
   \end{bmatrix}
   =
   \mathbf{J}^{-1}
   \begin{bmatrix}
   \dfrac{\partial N_I}{\partial \xi} \\
   \dfrac{\partial N_I}{\partial \eta} \\
   \dfrac{\partial N_I}{\partial \zeta}
   \end{bmatrix}.

The variation of the Green–Lagrange strain tensor in Voigt notation
can be written as

.. math::

   \delta \hat{\mathbf{E}}
   =
   \begin{bmatrix}
   \delta E_{\xi\xi} &
   \delta E_{\eta\eta} &
   \delta E_{\zeta\zeta} &
   2 \delta E_{\xi\eta} &
   2 \delta E_{\eta\zeta} &
   2 \delta E_{\xi\zeta}
   \end{bmatrix}^{\mathsf{T}}
   \approx
   \sum_{I=1}^{n_e} \mathbf{B}_I \, \delta \mathbf{u}_I,

where the strain–displacement matrix :math:`\mathbf{B}_I` is given by

.. math::

   \mathbf{B}_I^e
   =
   \begin{bmatrix}
   N_{I,\xi}\,\mathbf{g}_\xi^{\mathsf{T}} &
   N_{I,\eta}\,\mathbf{g}_\eta^{\mathsf{T}} &
   N_{I,\zeta}\,\mathbf{g}_\zeta^{\mathsf{T}} &
   N_{I,\xi}\,\mathbf{g}_\eta^{\mathsf{T}} + N_{I,\eta}\,\mathbf{g}_\xi^{\mathsf{T}} &
   N_{I,\eta}\,\mathbf{g}_\zeta^{\mathsf{T}} + N_{I,\zeta}\,\mathbf{g}_\eta^{\mathsf{T}} &
   N_{I,\xi}\,\mathbf{g}_\zeta^{\mathsf{T}} + N_{I,\zeta}\,\mathbf{g}_\xi^{\mathsf{T}}
   \end{bmatrix}^{\mathsf{T}}.

Here, :math:`\mathbf{g}_\xi`, :math:`\mathbf{g}_\eta`, and :math:`\mathbf{g}_\zeta`
are the covariant base vectors in the current configuration.

The internal virtual work discretized at the element level is

.. math::

   g_{\mathrm{int}}
   \approx
   \sum_{e=1}^{n_{\mathrm{el}}}
   \delta \mathbf{u}_e^{\mathsf{T}} \, \mathbf{f}_{\mathrm{int}}^{e},

with the element internal force vector

.. math::

   \mathbf{f}_{\mathrm{int}}^{e}
   =
   \int_{\Omega^{\square}}
   \mathbf{B}_e^{\mathsf{T}} \, \hat{\mathbf{S}} \, J_e
   \, \mathrm{d}\xi \, \mathrm{d}\eta \, \mathrm{d}\zeta,

where:

- :math:`\hat{\mathbf{S}}` is the second Piola–Kirchhoff stress in Voigt notation
- :math:`J_e` is the determinant of the Jacobian

The residual vector at the element level is defined as

.. math::

   \mathbf{f}_I
   =
   \mathbf{f}_{\mathrm{int},I}
   -
   \mathbf{f}_{\mathrm{ext},I}.

The linearization of the internal virtual work yields two contributions.

Material contribution:

.. math::

   \Delta g_{\mathrm{int}}^{\mathrm{mat}}
   \approx
   \delta \mathbf{u}^{\mathsf{T}}
   \sum_{e=1}^{n_{\mathrm{el}}}
   \mathbf{K}_m^{e}
   \delta \mathbf{u},

with

.. math::

   \mathbf{K}_m^{e}
   =
   \int_{\Omega^{\square}}
   \mathbf{B}_e^{\mathsf{T}} \, \hat{\mathbb{C}} \, \mathbf{B}_e \, J_e
   \, \mathrm{d}\xi \, \mathrm{d}\eta \, \mathrm{d}\zeta.

Geometric contribution:

.. math::

   \Delta g_{\mathrm{int}}^{\mathrm{geo}}
   \approx
   \sum_{e=1}^{n_{\mathrm{el}}}
   \sum_{I,J=1}^{n_e}
   \delta \mathbf{u}_I^{\mathsf{T}}
   \mathbf{K}_g^{IJ}
   \delta \mathbf{u}_J,

with

.. math::

   \mathbf{K}_g^{IJ}
   =
   \int_{\Omega^{\square}}
   S_{IJ} \, \mathbf{I}_3 \, J_e
   \, \mathrm{d}\xi \, \mathrm{d}\eta \, \mathrm{d}\zeta.

Here, :math:`S_{IJ}` is constructed from the second Piola–Kirchhoff stress
and derivatives of the shape functions.

When follower loads are present, the linearization introduces an additional
load stiffness matrix :math:`\mathbf{K}_l^{e}`, which is derived analogously
and subtracted from the tangent matrix. This contribution affects only elements
that are directly loaded.

The total tangent stiffness matrix is therefore given by

.. math::

   \mathbf{K}
   =
   \sum_{e=1}^{n_{\mathrm{el}}}
   \left(
      \mathbf{K}_m^{e}
      +
      \mathbf{K}_g^{e}
      -
      \mathbf{K}_l^{e}
   \right).