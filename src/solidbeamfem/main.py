# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np
import matplotlib.pyplot as plt

from pre_processing import (
    validate_inputs,
    create_mesh,
    impose_supports,
    compute_lame_parameters,
)
from processing import (
    newton_raphson_solver,
    shape_functions_3D,
)
from post_processing import compute_tip_displacement, plot_tip_displacement
import visualize as viz

# =============================================================================
# Problem setup
# =============================================================================

# Geometry (rectangular cross-section)
width = 1.0      # beam width
height = 0.5     # beam height
length = 20.0    # beam length

# Material (Saint Venant–Kirchhoff)
E = 1.2e7        # Young's modulus
nu = 0.0        # Poisson's ratio (nu=0 used to avoid Poisson/volumetric locking)

# Loading
load_type = "moment"     # "moment" or "shear"
k_max = 0.5             # target nondimensional load

# Discretization
numel = 10               # elements along beam length
ne_L = 2                 # nodes per element in axial direction (>= 2)

# Numerical integration
ngp_c = 2                # Gauss points per cross-sectional direction
ngp_l = ne_L             # Gauss points along beam length (<= 10)

# Locking alleviation (Assumed Natural Strain)
ANS_membrane = True      # membrane locking
ANS_shear = True         # transverse shear locking
ANS_curvature = True     # curvature–thickness locking

# Incremental-iterative solver
n_load_steps = 10        # load increments
max_iterations = 20      # Newton–Raphson iterations per step
tolerance = 1e-15        # convergence tolerance on energy norm

# Visualization
visualize_reference = True      # Visualize reference configuration
visualize_final = True          # Visualize final configuration
plot_displacement = True        # Plot normalized tip displacement 

# =============================================================================
# Pre-processing
# =============================================================================

# Validate inputs
validate_inputs(
    width, height, length,
    E, nu,
    load_type, k_max,
    numel, ne_L,
    ngp_c, ngp_l,
    n_load_steps, max_iterations, tolerance,
)

# Mesh: coordinates, connectivity, DOFs
coords, nnode_W, nnode_H, nnode_L, ndof, connect, ne = create_mesh(
    width, height, length, numel, ne_L
)

# Boundary conditions (clamped end)
prescribed_nodes = [0, 1, 2, 3]                  # clamped nodes
prescribed_dofs = impose_supports(prescribed_nodes)  # constrained DOFs

# Material parameters (Lamé)
lambda_, mu = compute_lame_parameters(E, nu)

# =============================================================================
# Visualize reference (undeformed) configuration
# =============================================================================

if visualize_reference:
    viz.plot_mesh_3D(coords, connect, prescribed_nodes)

# =============================================================================
# Solving the system
# =============================================================================

u, displacements_all = newton_raphson_solver(
    E, width, height, length, load_type, k_max, lambda_, mu,
    coords, connect, prescribed_dofs,
    numel, ne, ne_L, ngp_c, ngp_l,
    n_load_steps, max_iterations,
    tolerance, ANS_membrane,
    ANS_shear, ANS_curvature
)

# =============================================================================
# Post-processing
# =============================================================================

# Interpolate the displacement at the center of the tip cross-section
if plot_displacement:
    # Tip displacement for each load step
    tip_displacements = [
        compute_tip_displacement(u_step, connect, shape_functions_3D, ne_L)
        for u_step in displacements_all
    ]

    # Plot the load–displacement curve
    plot_tip_displacement(tip_displacements, length, load_type, k_max)

# =============================================================================
# Visualize final (deformed) configuration
# =============================================================================

if visualize_final:
    viz.plot_deformed_mesh_3D(
        coords, connect, u,
        scale=1.0,
        prescribed_nodes=prescribed_nodes,
    )