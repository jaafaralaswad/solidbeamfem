Assumed Natural Strain (ANS) Method
==================================

The Assumed Natural Strain (ANS) method is employed to alleviate
geometry-induced locking effects in solid-beam finite elements.
The formulation follows the procedure proposed by Caseiro *et al.* (2014),
originally developed for solid-shell elements, and is adapted here to beam
problems. The method is further extended to the
geometrically nonlinear regime.

Only the key ideas and final expressions are summarized here. Detailed
derivations are omitted, as they can be found in the original reference.

The ANS method targets the following locking mechanisms:

- Membrane locking
- Transverse shear locking
- Curvature–thickness locking

The formulation operates entirely at the element level and does not introduce
additional degrees of freedom.



Core Idea
---------

In the standard finite element formulation, strain components are evaluated
directly at Gauss integration points using the compatible strain–displacement
matrix. This can lead to locking because certain strain components become
over-constrained.

The ANS method alleviates locking by:

1. Evaluating selected strain components at predefined tying points
   in the parametric domain.
2. Interpolating these strains back to the Gauss points.
3. Replacing the corresponding rows of the compatible
   strain–displacement matrix with ANS-modified counterparts.



ANS Strain Definition
---------------------

Let :math:`E^c` denote a compatible strain component evaluated at a tying point
:math:`(\hat{\xi}_I,\eta,\zeta)`. The ANS strain evaluated at a Gauss point
:math:`(\xi_{\mathrm{GP}},\eta_{\mathrm{GP},\zeta_{\mathrm{GP})` is defined as

.. math::

   E^{\mathrm{ANS}}(\xi_{\mathrm{GP}},\eta_{\mathrm{GP},\zeta_{\mathrm{GP})
   =
   \sum_{I=1}^{n_t}
   \bar{N}_I(\xi_{\mathrm{GP}})
   \, E^c(\hat{\xi}_I,\eta_{\mathrm{GP},\zeta_{\mathrm{GP}),

where:

- :math:`\bar{N}_I` are interpolation functions associated with the tying points
- :math:`n_t` is the number of tying points

This leads to an ANS-modified strain–displacement matrix
:math:`\mathbf{B}^{\mathrm{ANS}}` such that

.. math::

   \mathbf{E}^{\mathrm{ANS}} = \mathbf{B}^{\mathrm{ANS}} \mathbf{u}.



Locking-Specific Modifications
------------------------------

Different locking modes are alleviated by modifying specific strain components:

- Membrane and transverse shear locking* 
  The axial strain components
  :math:`E_{\xi\xi}`, :math:`E_{\xi\eta}`, and :math:`E_{\xi\zeta}`
  are evaluated using reduced-order tying points along the beam axis.

- Curvature–thickness locking  
  The transverse normal strain components
  :math:`E_{\eta\eta}` and :math:`E_{\zeta\zeta}`
  are evaluated at nodal tying points, exploiting the vanishing of parasitic
  strain terms at nodes.

All remaining strain components are left unchanged.



Tangent Stiffness Consistency
-----------------------------

Because the strain–displacement matrix is modified, the material and geometric
stiffness matrices must be computed consistently using
:math:`\mathbf{B}^{\mathrm{ANS}}`. In the geometrically nonlinear regime, the
geometric stiffness contribution is reconstructed component-wise to preserve
variational consistency.



Summary
-------

The ANS method implemented in this work:

- Alleviates membrane, shear, and curvature–thickness locking
- Is applicable to geometrically nonlinear solid-beam formulations
- Preserves the standard finite element structure
- Requires only element-level modifications
- Allows selective activation of individual locking alleviation modes


More details
-------

.. [Caseiro2014]
   J. F. Caseiro, R. F. Valente, A. Reali, J. Kiendl,
   F. Auricchio, and R. Alves de Sousa,
   *On the Assumed Natural Strain method to alleviate locking
   in solid-shell NURBS-based finite elements*,
   Computational Mechanics, **53**, 1341–1353 (2014).
   https://doi.org/10.1007/s00466-014-0978-4
