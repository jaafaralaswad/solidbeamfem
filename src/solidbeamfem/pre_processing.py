"""Pre-processing utilities: input validation, mesh generation, supports, and material parameters."""
# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np

def validate_inputs(
    width, height, length,
    E, nu,
    load_type, k_max,
    numel, ne_L,
    ngp_c, ngp_l,
    n_load_steps, max_iterations,
    tolerance,
):
    """
    Validate input parameters for the finite element simulation.

    Raises
    ------
    ValueError
        If any input parameter is outside its admissible range.
    """
    # Geometry
    if width <= 0 or height <= 0 or length <= 0:
        raise ValueError(
            "All geometry dimensions (width, height, length) must be positive."
        )

    # Material properties
    if E <= 0:
        raise ValueError("Young's modulus 'E' must be positive.")
    if not (-1 < nu < 0.5):
        raise ValueError("Poisson's ratio 'nu' must lie in (-1, 0.5).")

    # Loading
    if load_type not in ("moment", "shear"):
        raise ValueError("Load type must be either 'moment' or 'shear'.")

    if load_type == "moment" and not (-1.0 <= k_max <= 1.0):
        raise ValueError(
            "For moment loading, |k_max| must not exceed 1. "
            "Larger values lead to self-penetration when the beam closes into a circle."
        )

    # Discretization
    if numel < 1:
        raise ValueError("'numel' (number of elements) must be at least 1.")
    if ne_L < 2:
        raise ValueError(
            "'ne_L' (nodes per element in axial direction) must be at least 2."
        )
    if not (1 <= ngp_c <= 10):
        raise ValueError(
            "'ngp_c' (Gauss points per cross-sectional direction) must be in [1, 10]."
        )
    if not (1 <= ngp_l <= 10):
        raise ValueError(
            "'ngp_l' (Gauss points in axial direction) must be in [1, 10]."
        )

    # Solver parameters
    if n_load_steps < 1:
        raise ValueError("'n_load_steps' must be at least 1.")
    if max_iterations < 1:
        raise ValueError("'max_iterations' must be at least 1.")
    if tolerance <= 0:
        raise ValueError("'tolerance' must be positive.")


def create_mesh(width, height, length, numel, ne_L):
    """
    Create a structured 3D solid-beam mesh (2x2 nodes in cross-section, ne_L nodes per element in length).

    Parameters
    ----------
    width : float
        Beam width (cross-sectional dimension, local y-direction).
    height : float
        Beam height (cross-sectional dimension, local z-direction).
    length : float
        Beam length (axial direction, global x-direction).
    numel : int
        Number of elements along the beam length.
    ne_L : int
        Number of nodes per element along the beam length (>= 2).

    Returns
    -------
    coords : (nnode, 3) ndarray
        Nodal coordinates ordered along length first, then height, then width.
        Columns are [x, y, z] = [length, width, height].
    nnode_W, nnode_H, nnode_L : int
        Number of nodes in width, height, and length directions, respectively.
    ndof : int
        Total number of degrees of freedom (3 * nnode).
    connect : (numel, ne) ndarray of int
        Element connectivity (global node indices, 0-based).
    ne : int
        Total number of nodes per element (4 * ne_L).
    """
    ne = 4 * ne_L  # nodes per element (2x2 cross-section times ne_L along length)

    # Nodes per direction (fixed 2x2 in cross-section)
    nnode_W = 2
    nnode_H = 2
    nnode_L = numel + 1 + (ne_L - 2) * numel  # = 1 + numel*(ne_L - 1)

    nnode = nnode_W * nnode_H * nnode_L
    ndof = 3 * nnode

    # Element sizes and nodal spacing
    elen_W = width
    elen_H = height
    elen_L = length / numel
    nelen_L = elen_L / (ne_L - 1)  # node spacing along length

    # Nodal coordinates
    coords = np.zeros((nnode, 3))
    node_id = 0

    for i in range(nnode_L):                 # length direction
        x_L = nelen_L * i
        for j in range(nnode_H):             # height direction
            x_H = elen_H * j
            for k in range(nnode_W):         # width direction
                x_W = elen_W * k

                coords[node_id, 0] = x_L     # x (length)
                coords[node_id, 1] = x_W     # y (width)
                coords[node_id, 2] = x_H     # z (height)
                node_id += 1

    # Connectivity (0-based global node indices)
    connect = np.zeros((numel, ne), dtype=int)
    stride_cs = nnode_W * nnode_H  # nodes per cross-sectional slice

    for i in range(numel):
        base = i * stride_cs + i * (ne_L - 2) * stride_cs
        for j in range(ne):
            connect[i, j] = j + base

    return coords, nnode_W, nnode_H, nnode_L, ndof, connect, ne

def impose_supports(prescribed_nodes):
    """
    Construct the list of constrained degrees of freedom (DOFs)
    corresponding to prescribed nodes.
    """
    num_nodes = len(prescribed_nodes)
    prescribed_dofs = np.zeros(3 * num_nodes, dtype=int)

    # Each node has 3 translational DOFs: (x, y, z)
    for i, node in enumerate(prescribed_nodes):
        prescribed_dofs[3*i    ] = 3*node      # x-DOF
        prescribed_dofs[3*i + 1] = 3*node + 1  # y-DOF
        prescribed_dofs[3*i + 2] = 3*node + 2  # z-DOF

    return prescribed_dofs

def compute_lame_parameters(E, nu):
    """
    Compute the Lamé parameters for isotropic linear elasticity.

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.

    Returns
    -------
    lambda_ : float
        First Lamé parameter.
    mu : float
        Second Lamé parameter (shear modulus).
    """
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))  # first Lamé parameter
    mu = E / (2 * (1 + nu))                       # shear modulus

    return lambda_, mu