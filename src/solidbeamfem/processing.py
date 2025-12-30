"""
Processing routines for the nonlinear 3D solid-beam solver.

This module contains the Newton–Raphson driver, global assembly, element routines,
loading (including follower-load stiffness), and small numerical utilities.
"""

# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np


def assemble_global_system(
    coords, connect, u,
    numel, ne, ne_L,
    lambda_, mu,
    ngp_c, ngp_l,
    ke_l,
    ANS_membrane, ANS_shear, ANS_curvature,
):
    """
    Assemble the global stiffness matrix and internal force vector at the current state.

    Parameters
    ----------
    coords : (nnode, 3) ndarray
        Reference nodal coordinates.
    connect : (numel, ne) ndarray of int
        Element connectivity (0-based global node indices).
    u : (ndof, 1) ndarray
        Global displacement vector.
    numel, ne, ne_L, ngp_c, ngp_l : int
        Discretization and quadrature parameters.
    lambda_, mu : float
        Lamé parameters.
    ke_l : (ndof_e, ndof_e) ndarray
        Load stiffness contribution (added to the last element only).
    ANS_membrane, ANS_shear, ANS_curvature : bool
        ANS toggles.

    Returns
    -------
    K : (ndof, ndof) ndarray
        Global stiffness matrix.
    F_int : (ndof, 1) ndarray
        Global internal force vector.
    """
    ndof = coords.shape[0] * 3
    ndof_e = 3 * ne

    K = np.zeros((ndof, ndof))
    F_int = np.zeros((ndof, 1))

    for ie in range(numel):
        x0_e = np.zeros((ndof_e, 3))  # reference nodal positions (element)
        u_e = np.zeros((ndof_e, 3))   # nodal displacements (element)

        # Gather element nodal coordinates and displacements
        for a in range(ne):
            node = connect[ie, a]
            x0_e[a, :] = coords[node, :]
            u_e[a, :] = u[3*node:3*node+3, 0]

        ke, fe = element_3D_nonlinear(
            ne, ne_L, x0_e, u_e, lambda_, mu,
            ngp_c, ngp_l,
            ANS_membrane, ANS_shear, ANS_curvature,
        )

        # Add follower-load stiffness only to the loaded (last) element
        if ie == numel - 1:
            ke += ke_l

        # Assemble global stiffness and internal force
        for a in range(ne):
            ia = connect[ie, a]
            F_int[3*ia:3*ia+3] += fe[3*a:3*a+3]

            for b in range(ne):
                ib = connect[ie, b]
                K[3*ia:3*ia+3, 3*ib:3*ib+3] += ke[3*a:3*a+3, 3*b:3*b+3]

    return K, F_int


def compute_convergence_metrics(R_r, R, du):
    """
    Compute convergence metrics used by the Newton–Raphson loop.

    Returns
    -------
    res_norm : float
        Euclidean norm of the reduced residual.
    energy_norm : float
        Energy norm: |R^T du|.
    """
    res_norm = np.linalg.norm(R_r)
    energy_norm = abs(R.T @ du)[0, 0]
    return res_norm, energy_norm


def newton_raphson_solver(
    E, width, height, length, load_type, k_max, lambda_, mu,
    coords, connect, prescribed_dofs,
    numel, ne, ne_L, ngp_c, ngp_l,
    n_load_steps, max_iterations, tolerance,
    ANS_membrane, ANS_shear, ANS_curvature,
):
    """
    Newton–Raphson solver for incremental loading of a nonlinear 3D cantilever beam.

    Parameters
    ----------
    E, width, height, length : float
        Material and geometric properties used by the load definition.
    load_type : {"moment", "shear"}
        Loading type.
    k_max : float
        Target nondimensional load parameter.
    lambda_, mu : float
        Lamé parameters.
    coords : (nnode, 3) ndarray
        Reference nodal coordinates.
    connect : (numel, ne) ndarray of int
        Element connectivity (0-based global node indices).
    prescribed_dofs : array-like of int
        Constrained DOF indices.
    numel : int
        Number of elements.
    ne : int
        Nodes per element.
    ne_L : int
        Nodes per element along the axial direction.
    ngp_c, ngp_l : int
        Gauss points in cross-sectional and axial directions.
    n_load_steps : int
        Number of load increments.
    max_iterations : int
        Maximum Newton–Raphson iterations per increment.
    tolerance : float
        Convergence tolerance on the energy norm.
    ANS_membrane, ANS_shear, ANS_curvature : bool
        Assumed Natural Strain (ANS) toggles.

    Returns
    -------
    u : (ndof, 1) ndarray
        Final displacement vector.
    displacements_all : list of (ndof, 1) ndarray
        Converged displacement vector at each load step.
    """
    ndof = coords.shape[0] * 3  # 3 DOFs per node

    u = np.zeros((ndof, 1))

    # Convergence history
    ResNorm = np.zeros((n_load_steps, max_iterations))   # residual norm
    EnerNorm = np.zeros((n_load_steps, max_iterations))  # energy norm

    displacements_all = []

    # Incremental loading loop
    for inc in range(1, n_load_steps + 1):

        # Newton–Raphson iterations
        for itn in range(1, max_iterations + 1):

            # External force (follower load updated every iteration for consistency)
            F_ext, ke_l = compute_external_force(
                E, length, height, width, load_type, k_max,
                n_load_steps, inc,
                coords, connect, ne, ne_L, numel,
                ngp_c, ndof, u,
            )

            # Assemble stiffness and internal force
            K, F_int = assemble_global_system(
                coords, connect, u,
                numel, ne, ne_L,
                lambda_, mu,
                ngp_c, ngp_l,
                ke_l,
                ANS_membrane, ANS_shear, ANS_curvature,
            )

            # Residual (internal - external)
            R = F_int - F_ext

            # Linear solve with displacement constraints
            du, R_r = solve_linear_system(K, R, prescribed_dofs)

            # Update solution
            u += du

            # Convergence metrics
            res_norm, energy_norm = compute_convergence_metrics(R_r, R, du)
            ResNorm[inc - 1, itn - 1] = res_norm
            EnerNorm[inc - 1, itn - 1] = energy_norm

            print(
                f"Load Step {inc}, Iteration {itn}: "
                f"Residual Norm = {res_norm:10.6e}, "
                f"Energy Norm = {energy_norm:10.6e}"
            )

            if energy_norm < tolerance:
                print("=====")
                break

        displacements_all.append(u.copy())

    return u, displacements_all


def compute_external_force(
    E, length, height, width,
    load_type, k_max,
    n_load_steps, inc,
    coords, connect,
    ne, ne_L, numel,
    ngp_c, ndof, u,
):
    """
    Compute the external force vector and the load-stiffness contribution.

    Supported loading types
    -----------------------
    "moment"
        Follower bending moment applied via surface traction on the end face.
        The follower-load traction and the associated load stiffness follow
        Wriggers (2008), Section 4.2.5.
    "shear"
        Conservative transverse shear force on the end face (ke_l = 0).

    Returns
    -------
    F_ext : (ndof, 1) ndarray
        External force vector.
    ke_l : (ndof_e, ndof_e) ndarray
        Load stiffness matrix contribution for the loaded element.
    """
    ndof_e = 3 * ne
    F_ext = np.zeros((ndof, 1))
    ke_l = np.zeros((ndof_e, ndof_e))

    # Loaded region: only the last element (end face)
    loaded_elements = [numel - 1]

    # Current nodal coordinates (reference + displacement)
    x_current = coords + u.reshape(-1, 3)

    # -------------------------------------------------------------------------
    # MOMENT LOADING (follower surface traction on end face)
    # -------------------------------------------------------------------------
    if load_type == "moment":

        # Gauss quadrature on the end face (cross-section)
        gauss_pts, gauss_wts = gauss(ngp_c)
        xi_end = 1.0  # end face: xi = 1
        eta = gauss_pts
        zeta = gauss_pts

        # Rectangular section second moment of area
        Iyy = height**3 * width / 12.0

        # Incremental bending moment (nondimensional scaling used in benchmark)
        k = k_max * (inc / n_load_steps)
        M = (2.0 * np.pi * E * Iyy / length) * k

        for ie in loaded_elements:
            fe_local = np.zeros((ndof_e, 1))

            for j_eta in range(ngp_c):
                for j_zeta in range(ngp_c):

                    N, dN = shape_functions_3D(ne_L, xi_end, eta[j_eta], zeta[j_zeta])

                    # Build the interpolation matrix for traction integration
                    Nmat = np.zeros((3, ndof_e))

                    # Interpolate the reference z-coordinate at the Gauss point
                    z_gp = 0.0
                    for a in range(ne):
                        node = connect[ie, a]
                        Nmat[0, 3*a]     = N[a]
                        Nmat[1, 3*a + 1] = N[a]
                        Nmat[2, 3*a + 2] = N[a]
                        z_gp += N[a] * coords[node, 2]

                    # Linear traction distribution over height to represent bending moment
                    lever_arm = z_gp - height / 2.0
                    traction = (-12.0 * M / height**3) * lever_arm

                    # Tangents on the end face (used for normal and follower-load stiffness)
                    t_eta = np.zeros((3, 1))
                    t_zeta = np.zeros((3, 1))
                    for a in range(ne):
                        node = connect[ie, a]
                        t_eta  += (dN[1, a] * x_current[node, :]).reshape(3, 1)
                        t_zeta += (dN[2, a] * x_current[node, :]).reshape(3, 1)

                    # Surface normal = t_eta x t_zeta
                    normal = np.zeros((3, 1))
                    normal[0] = t_eta[1]*t_zeta[2] - t_eta[2]*t_zeta[1]
                    normal[1] = t_eta[2]*t_zeta[0] - t_eta[0]*t_zeta[2]
                    normal[2] = t_eta[0]*t_zeta[1] - t_eta[1]*t_zeta[0]

                    # Auxiliary matrices for follower-load stiffness (Wriggers 2008, Sec. 4.2.5)
                    Nbar_eta = np.zeros((3, 3))
                    Nbar_eta[0, 1] = t_eta[2]
                    Nbar_eta[1, 0] = -t_eta[2]
                    Nbar_eta[0, 2] = -t_eta[1]
                    Nbar_eta[2, 0] = t_eta[1]
                    Nbar_eta[1, 2] = t_eta[0]
                    Nbar_eta[2, 1] = -t_eta[0]

                    Nbar_zeta = np.zeros((3, 3))
                    Nbar_zeta[0, 1] = t_zeta[2]
                    Nbar_zeta[1, 0] = -t_zeta[2]
                    Nbar_zeta[0, 2] = -t_zeta[1]
                    Nbar_zeta[2, 0] = t_zeta[1]
                    Nbar_zeta[1, 2] = t_zeta[0]
                    Nbar_zeta[2, 1] = -t_zeta[0]

                    w = gauss_wts[j_eta] * gauss_wts[j_zeta]

                    # Load stiffness contribution (only for the loaded element)
                    for a in range(ne):
                        for b in range(ne):
                            ke_l[3*a:3*a+3, 3*b:3*b+3] -= (
                                traction * N[a] *
                                (dN[1, b] * Nbar_zeta - dN[2, b] * Nbar_eta) * w
                            )

                    # External force contribution
                    fe_local += (Nmat.T @ (traction * normal)) * w

            # Scatter local external force to global vector
            for a in range(ne):
                node = connect[ie, a]
                F_ext[3*node:3*node+3, 0] += fe_local[3*a:3*a+3, 0]

    # -------------------------------------------------------------------------
    # SHEAR LOADING (conservative nodal force on end face)
    # -------------------------------------------------------------------------
    elif load_type == "shear":

        shear_force = (1.0 / length**2) * (inc / n_load_steps)

        for ie in loaded_elements:
            fe_local = np.zeros((ndof_e, 1))

            # Apply to the top 4 nodes of the loaded end slice (z-direction)
            for k in range(4):
                a = ne - 1 - k
                fe_local[3*a + 2, 0] = shear_force

            for a in range(ne):
                node = connect[ie, a]
                F_ext[3*node:3*node+3, 0] += fe_local[3*a:3*a+3, 0]

        # ke_l remains zero for conservative shear load

    return F_ext, ke_l


def element_3D_nonlinear(
    ne, ne_L,
    x0_e, u_e,
    lambda_, mu,
    ngp_c, ngp_l,
    ANS_membrane, ANS_shear, ANS_curvature,
):
    """
    Compute element tangent stiffness matrix and internal force vector, using ANS to alleviate
    membrane, transverse shear, and curvature–thickness locking.
    """
    ndof_e = 3 * ne  # Degrees of freedom per element

    # ------------------------------------- helpers ---------------------------------
    def _covariant_basis_reference(dN, x0):
        """Reference covariant basis (Jacobian) and metric tensors."""
        H = np.zeros((3, 3))
        for a in range(ne):
            H[:, 0] += dN[:, a] * x0[a, 0]
            H[:, 1] += dN[:, a] * x0[a, 1]
            H[:, 2] += dN[:, a] * x0[a, 2]

        detJ = np.linalg.det(H)

        G_1 = H[:, 0]
        G_2 = H[:, 1]
        G_3 = H[:, 2]

        Gco = np.array([
            [np.dot(G_1, G_1), np.dot(G_1, G_2), np.dot(G_1, G_3)],
            [np.dot(G_2, G_1), np.dot(G_2, G_2), np.dot(G_2, G_3)],
            [np.dot(G_3, G_1), np.dot(G_3, G_2), np.dot(G_3, G_3)],
        ])
        Gcontra = np.linalg.inv(Gco)
        return detJ, Gco, Gcontra

    def _covariant_basis_current(dN, xcur):
        """Current covariant basis vectors and covariant metric tensor."""
        H = np.zeros((3, 3))
        for a in range(ne):
            H[0, :] += dN[:, a] * xcur[a, 0]
            H[1, :] += dN[:, a] * xcur[a, 1]
            H[2, :] += dN[:, a] * xcur[a, 2]

        g_1 = H[:, 0]
        g_2 = H[:, 1]
        g_3 = H[:, 2]

        gco = np.array([
            [np.dot(g_1, g_1), np.dot(g_1, g_2), np.dot(g_1, g_3)],
            [np.dot(g_2, g_1), np.dot(g_2, g_2), np.dot(g_2, g_3)],
            [np.dot(g_3, g_1), np.dot(g_3, g_2), np.dot(g_3, g_3)],
        ])
        return g_1, g_2, g_3, gco

    def _stress_svk(Gcontra, E_GL):
        """2nd Piola–Kirchhoff stress for SVK material: S = C : E."""
        S = np.zeros((3, 3))
        for I in range(3):
            for J in range(3):
                acc = 0.0
                for Kidx in range(3):
                    for Lidx in range(3):
                        acc += compute_C_SVK_component(Gcontra, lambda_, mu, I, J, Kidx, Lidx) * E_GL[Kidx, Lidx]
                S[I, J] = acc
        return S

    def _stress_to_voigt(S):
        """Cast symmetric stress tensor to Voigt notation (6x1)."""
        return np.array([S[0, 0], S[1, 1], S[2, 2], S[0, 1], S[0, 2], S[1, 2]]).reshape(6, 1)

    def _build_B_standard(dN, g_1, g_2, g_3):
        """Standard strain–displacement matrix Be (6 x 3*ne)."""
        Be = np.zeros((6, ndof_e))
        for a in range(ne):
            B = np.zeros((6, 3))
            B[0, :] = g_1 * dN[0, a]
            B[1, :] = g_2 * dN[1, a]
            B[2, :] = g_3 * dN[2, a]
            B[3, :] = g_2 * dN[0, a] + g_1 * dN[1, a]
            B[4, :] = g_3 * dN[0, a] + g_1 * dN[2, a]
            B[5, :] = g_3 * dN[1, a] + g_2 * dN[2, a]
            Be[:, 3*a:3*a+3] = B
        return Be

    def _material_tangent_voigt(Gcontra):
        """Material tangent in Voigt notation (6x6, symmetric)."""
        matC = np.zeros((6, 6))
        for I in range(3):
            for J in range(I, 3):
                for Kidx in range(3):
                    for Lidx in range(Kidx, 3):
                        row = voigt_index(I, J)
                        col = voigt_index(Kidx, Lidx)
                        matC[row, col] = compute_C_SVK_component(Gcontra, lambda_, mu, I, J, Kidx, Lidx)
                        if row != col:
                            matC[col, row] = matC[row, col]
        return matC

    # Initialize element material and geometric stiffnesses and internal force vector
    ke_m = np.zeros((ndof_e, ndof_e))
    ke_g = np.zeros((ndof_e, ndof_e))
    fe = np.zeros((ndof_e, 1))

    # Gauss points
    eta, wgp_c = gauss(ngp_c)
    zeta = eta
    xi, wgp_l = gauss(ngp_l)

    # Current configuration
    x_current = x0_e + u_e

    # Loop over Gauss points
    for k1 in range(ngp_l):
        for k2 in range(ngp_c):
            for k3 in range(ngp_c):

                # Evaluate shape function derivatives at Gauss point
                _, dN = shape_functions_3D(ne_L, xi[k1], eta[k2], zeta[k3])

                # Reference metric tensors and Jacobian determinant
                detJ, Gco, Gcontra = _covariant_basis_reference(dN, x0_e)

                # Current metric tensors
                g_1, g_2, g_3, gco = _covariant_basis_current(dN, x_current)

                # Green-Lagrange strain tensor
                E_GL = 0.5 * (gco - Gco)

                # 2nd Piola-Kirchhoff stress (SVK- material)
                S = _stress_svk(Gcontra, E_GL)

                # Cast stress in Voigt notation
                Svgt = _stress_to_voigt(S)

                # Calculate standard strain-displacement matrix Be
                Be = _build_B_standard(dN, g_1, g_2, g_3)

                # ANS - refer to accompanying documentation for details
                # Modify standard Be using ANS for membrane and transverse shear locking
                if ANS_membrane or ANS_shear:

                    ntying1 = ne_L - 1  # Number of tying points based on reduced integration
                    xibar1, _ = gauss(ntying1)  # Tying points in xi
                    Nmat1 = lagrange_shape_1D(ntying1, xi[k1])  # Evaluate 1D Lagrange shape functions for tying points

                    M1 = np.zeros((ntying1, ntying1))  # Initialize Matrix M
                    for i in range(ntying1):  # Calculate Matrix M
                        M1[i, :] = lagrange_shape_1D(ntying1, xibar1[i])
                    invM1 = np.linalg.inv(M1)

                    if ANS_membrane:  # Initialize modified 1st row for membrane
                        Bebar_1 = np.zeros((ntying1, 3 * ne))

                    if ANS_shear:  # Initialize modified 4th and 5th rows for transverse shear
                        Bebar_4 = np.zeros((ntying1, 3 * ne))
                        Bebar_5 = np.zeros((ntying1, 3 * ne))

                    for k in range(ntying1):
                        _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zeta[k3])  # Evaluate shape functions at tying points

                        # Covariant basis in current configuration at tying points
                        H_tying1 = np.zeros((3, 3))
                        for m in range(ne):
                            H_tying1[0, :] += dNbar1[0:3, m] * x_current[m, 0]
                            H_tying1[1, :] += dNbar1[0:3, m] * x_current[m, 1]
                            H_tying1[2, :] += dNbar1[0:3, m] * x_current[m, 2]

                        # Construct covariant basis vectors at tying points
                        g_tying_1 = H_tying1[:, 0]
                        g_tying_2 = H_tying1[:, 1]
                        g_tying_3 = H_tying1[:, 2]

                        # Construct the ANS modified B matrix for membrane locking
                        if ANS_membrane:
                            for m in range(ne):
                                Bbar1 = np.zeros(3)
                                Bbar1[0] = g_tying_1[0] * dNbar1[0, m]
                                Bbar1[1] = g_tying_1[1] * dNbar1[0, m]
                                Bbar1[2] = g_tying_1[2] * dNbar1[0, m]
                                Bebar_1[k, 3*m:3*m+3] = Bbar1

                        # Construct the ANS modified B matrix for transverse shear locking
                        if ANS_shear:
                            for m in range(ne):
                                Bbar4 = np.zeros(3)
                                Bbar4[0] = g_tying_2[0] * dNbar1[0, m] + g_tying_1[0] * dNbar1[1, m]
                                Bbar4[1] = g_tying_2[1] * dNbar1[0, m] + g_tying_1[1] * dNbar1[1, m]
                                Bbar4[2] = g_tying_2[2] * dNbar1[0, m] + g_tying_1[2] * dNbar1[1, m]
                                Bebar_4[k, 3*m:3*m+3] = Bbar4

                                Bbar5 = np.zeros(3)
                                Bbar5[0] = g_tying_3[0] * dNbar1[0, m] + g_tying_1[0] * dNbar1[2, m]
                                Bbar5[1] = g_tying_3[1] * dNbar1[0, m] + g_tying_1[1] * dNbar1[2, m]
                                Bbar5[2] = g_tying_3[2] * dNbar1[0, m] + g_tying_1[2] * dNbar1[2, m]
                                Bebar_5[k, 3*m:3*m+3] = Bbar5

                    # Replace 1st row for membrane
                    if ANS_membrane:
                        Row1 = Nmat1 @ (invM1 @ Bebar_1)
                        Be[0, :] = Row1

                    # Replace 4th and 5th rows for transverse shear
                    if ANS_shear:
                        Row4 = Nmat1 @ (invM1 @ Bebar_4)
                        Be[3, :] = Row4
                        Row5 = Nmat1 @ (invM1 @ Bebar_5)
                        Be[4, :] = Row5

                # Modify standard Be using ANS for curvature-thickness locking
                if ANS_curvature:

                    ntying2 = ne_L  # Number of tying points based on equidistant nodes
                    xibar2 = np.linspace(-1.0, 1.0, ntying2)  # Use equally spaced tying points

                    Nmat2 = lagrange_shape_1D(ntying2, xi[k1])  # Tying points in xi
                    M2 = np.eye(ntying2)  # Matrix M collapses to identity - equidistant nodes
                    invM2 = M2

                    # Initialize modified 2nd and 3rd rows of B matrix
                    Bebar_2 = np.zeros((ntying2, 3 * ne))
                    Bebar_3 = np.zeros((ntying2, 3 * ne))

                    for k in range(ntying2):
                        _, dNbar2 = shape_functions_3D(ne_L, xibar2[k], eta[k2], zeta[k3])  # Evaluate shape functions derivatives at tying points

                        # Covariant basis in current configuration at tying points
                        H_tying2 = np.zeros((3, 3))
                        for m in range(ne):
                            H_tying2[0, :] += dNbar2[0:3, m] * x_current[m, 0]
                            H_tying2[1, :] += dNbar2[0:3, m] * x_current[m, 1]
                            H_tying2[2, :] += dNbar2[0:3, m] * x_current[m, 2]

                        # Construct covariant vectors at tying points in current configuration
                        g_tying_1 = H_tying2[:, 0]
                        g_tying_2 = H_tying2[:, 1]
                        g_tying_3 = H_tying2[:, 2]

                        # Calculate 2nd row of B matrix
                        for m in range(ne):
                            Bbar2 = np.zeros(3)
                            Bbar2[0] = g_tying_2[0] * dNbar2[1, m]
                            Bbar2[1] = g_tying_2[1] * dNbar2[1, m]
                            Bbar2[2] = g_tying_2[2] * dNbar2[1, m]
                            Bebar_2[k, 3*m:3*m+3] = Bbar2

                        # Calculate 3rd row of B matrix
                        for m in range(ne):
                            Bbar3 = np.zeros(3)
                            Bbar3[0] = g_tying_3[0] * dNbar2[2, m]
                            Bbar3[1] = g_tying_3[1] * dNbar2[2, m]
                            Bbar3[2] = g_tying_3[2] * dNbar2[2, m]
                            Bebar_3[k, 3*m:3*m+3] = Bbar3

                    # Calculate and replace 2nd row of B matrix
                    Row2 = Nmat2 @ (invM2 @ Bebar_2)
                    Be[1, :] = Row2

                    # Calculate and replace 3rd row of B matrix
                    Row3 = Nmat2 @ (invM2 @ Bebar_3)
                    Be[2, :] = Row3

                weight = detJ * wgp_l[k1] * wgp_c[k2] * wgp_c[k3]

                # Internal force vector
                fe += Be.T @ Svgt * weight

                # Material stiffness matrix
                matC = _material_tangent_voigt(Gcontra)

                # Material stiffness contribution
                ke_m += Be.T @ matC @ Be * weight

                # Geometric stiffness contribution

                # Initialize sub-matrices
                ke_g11 = np.zeros((ndof_e, ndof_e))
                ke_g12 = np.zeros((ndof_e, ndof_e))
                ke_g13 = np.zeros((ndof_e, ndof_e))
                ke_g21 = np.zeros((ndof_e, ndof_e))
                ke_g22 = np.zeros((ndof_e, ndof_e))
                ke_g23 = np.zeros((ndof_e, ndof_e))
                ke_g31 = np.zeros((ndof_e, ndof_e))
                ke_g32 = np.zeros((ndof_e, ndof_e))
                ke_g33 = np.zeros((ndof_e, ndof_e))

                # Initialize tying arrays for membrane
                if ANS_membrane:
                    tying_dNxi_dNxi = np.zeros((ntying1, 1))

                # Initialize tying arrays for transverse shear
                if ANS_shear:
                    tying_dNxi_dNeta = np.zeros((ntying1, 1))
                    tying_dNeta_dNxi = np.zeros((ntying1, 1))
                    tying_dNxi_dNzeta = np.zeros((ntying1, 1))
                    tying_dNzeta_dNxi = np.zeros((ntying1, 1))

                # Initialize tying arrays for curvature-thickness
                if ANS_curvature:
                    tying_dNeta_dNeta = np.zeros((ntying2, 1))
                    tying_dNzeta_dNzeta = np.zeros((ntying2, 1))

                I3 = np.eye(3)

                # Loop over element nodes
                for K in range(ne):
                    for L in range(ne):

                        # Lump constants and define identity matrix
                        factor = weight

                        # In-plane shear contributions always calculated the standard way (never modified by ANS)
                        ke_g23[3*K:3*K+3, 3*L:3*L+3] += dN[1, K] * S[1, 2] * I3 * dN[2, L] * factor
                        ke_g32[3*K:3*K+3, 3*L:3*L+3] += dN[2, K] * S[2, 1] * I3 * dN[1, L] * factor

                        # --- Membrane contribution to stress ---
                        if ANS_membrane:  # calculated stresses in a modified way if alleviated
                            for k in range(ntying1):
                                _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zeta[k3])
                                tying_dNxi_dNxi[k, 0] = dNbar1[0, K] * dNbar1[0, L]
                            proj_dNxi_dNxi = Nmat1 @ (invM1 @ tying_dNxi_dNxi)
                            Sxi_xi = S[0, 0] * proj_dNxi_dNxi
                        else:  # Calculate stresses the standard way if not alleviated
                            Sxi_xi = S[0, 0] * dN[0, K] * dN[0, L]
                        ke_g11[3*K:3*K+3, 3*L:3*L+3] += Sxi_xi * I3 * factor  # Calculate membrane geometric contribution

                        # --- Transverse shear contribution to stress ---
                        if ANS_shear:  # calculated stresses in a modified way if alleviated
                            for k in range(ntying1):
                                _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zeta[k3])
                                tying_dNxi_dNeta[k, 0] = dNbar1[0, K] * dNbar1[1, L]
                                tying_dNeta_dNxi[k, 0] = dNbar1[1, K] * dNbar1[0, L]
                                tying_dNxi_dNzeta[k, 0] = dNbar1[0, K] * dNbar1[2, L]
                                tying_dNzeta_dNxi[k, 0] = dNbar1[2, K] * dNbar1[0, L]

                            proj_dNxi_dNeta = Nmat1 @ (invM1 @ tying_dNxi_dNeta)
                            proj_dNeta_dNxi = Nmat1 @ (invM1 @ tying_dNeta_dNxi)
                            proj_dNxi_dNzeta = Nmat1 @ (invM1 @ tying_dNxi_dNzeta)
                            proj_dNzeta_dNxi = Nmat1 @ (invM1 @ tying_dNzeta_dNxi)

                            Sxi_eta = S[0, 1] * proj_dNxi_dNeta
                            Seta_xi = S[1, 0] * proj_dNeta_dNxi
                            Sxi_zeta = S[0, 2] * proj_dNxi_dNzeta
                            Szeta_xi = S[2, 0] * proj_dNzeta_dNxi
                        else:  # Calculate stresses the standard way if not alleviated
                            Sxi_eta = S[0, 1] * dN[0, K] * dN[1, L]
                            Seta_xi = S[1, 0] * dN[1, K] * dN[0, L]
                            Sxi_zeta = S[0, 2] * dN[0, K] * dN[2, L]
                            Szeta_xi = S[2, 0] * dN[2, K] * dN[0, L]

                        # Calculate transvere shear geometric contributions
                        ke_g12[3*K:3*K+3, 3*L:3*L+3] += Sxi_eta * I3 * factor
                        ke_g21[3*K:3*K+3, 3*L:3*L+3] += Seta_xi * I3 * factor
                        ke_g13[3*K:3*K+3, 3*L:3*L+3] += Sxi_zeta * I3 * factor
                        ke_g31[3*K:3*K+3, 3*L:3*L+3] += Szeta_xi * I3 * factor

                        # --- Curvature-thickness contribution to stress ---
                        if ANS_curvature:
                            for k in range(ntying2):
                                _, dNbar2 = shape_functions_3D(ne_L, xibar2[k], eta[k2], zeta[k3])
                                tying_dNeta_dNeta[k, 0] = dNbar2[1, K] * dNbar2[1, L]
                                tying_dNzeta_dNzeta[k, 0] = dNbar2[2, K] * dNbar2[2, L]

                            proj_dNeta_dNeta = Nmat2 @ (invM2 @ tying_dNeta_dNeta)
                            proj_dNzeta_dNzeta = Nmat2 @ (invM2 @ tying_dNzeta_dNzeta)

                            Seta_eta = S[1, 1] * proj_dNeta_dNeta
                            Szeta_zeta = S[2, 2] * proj_dNzeta_dNzeta
                        else:  # Calculate stresses the standard way if not alleviated
                            Seta_eta = S[1, 1] * dN[1, K] * dN[1, L]
                            Szeta_zeta = S[2, 2] * dN[2, K] * dN[2, L]

                        # Calculate curvature-thickness geometric contributions
                        ke_g22[3*K:3*K+3, 3*L:3*L+3] += Seta_eta * I3 * factor
                        ke_g33[3*K:3*K+3, 3*L:3*L+3] += Szeta_zeta * I3 * factor

                # Sum all parts to form final geometric stiffness matrix
                ke_g += (
                    ke_g11 + ke_g12 + ke_g13 +
                    ke_g21 + ke_g22 + ke_g23 +
                    ke_g31 + ke_g32 + ke_g33
                )

    # Add material and geometric parts of element stiffness matrix
    ke = ke_m + ke_g
    return ke, fe

def gauss(ngp):
    """
    Return Gauss–Legendre quadrature points and weights on [-1, 1].

    Parameters
    ----------
    ngp : int
        Number of Gauss points (1 to 10 supported).

    Returns
    -------
    xi : (ngp,) ndarray
        Gauss points.
    w : (ngp,) ndarray
        Corresponding weights.
    """
    if not (1 <= ngp <= 10):
        raise ValueError("'ngp' must be an integer in [1, 10].")

    data = {
        1: ([0.0], [2.0]),
        2: ([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)], [1.0, 1.0]),
        3: (
            [-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)],
            [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0],
        ),
        4: (
            [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116],
            [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451],
        ),
        5: (
            [-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459],
            [0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851],
        ),
        6: (
            [-0.9324695142, -0.6612093865, -0.2386191861,
              0.2386191861,  0.6612093865,  0.9324695142],
            [0.1713244924, 0.3607615730, 0.4679139346,
             0.4679139346, 0.3607615730, 0.1713244924],
        ),
        7: (
            [-0.9491079123, -0.7415311856, -0.4058451514,
              0.0,
              0.4058451514, 0.7415311856, 0.9491079123],
            [0.1294849662, 0.2797053915, 0.3818300505,
             0.4179591837,
             0.3818300505, 0.2797053915, 0.1294849662],
        ),
        8: (
            [-0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425,
              0.1834346425,  0.5255324099,  0.7966664774,  0.9602898565],
            [0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
             0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363],
        ),
        9: (
            [-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0.0,
              0.3242534234, 0.6133714327, 0.8360311073, 0.9681602395],
            [0.0812743884, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550,
             0.3123470770, 0.2606106964, 0.1806481607, 0.0812743884],
        ),
        10: (
            [-0.9739065285, -0.8650633667, -0.6794095683, -0.4333953941, -0.1488743390,
              0.1488743390,  0.4333953941,  0.6794095683,  0.8650633667,  0.9739065285],
            [0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193, 0.2955242247,
             0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443],
        ),
    }

    xi, w = data[ngp]
    return np.asarray(xi, dtype=float), np.asarray(w, dtype=float)

def solve_linear_system(K, R, prescribed_dofs):
    """
    Solve the linear system K d = -R with Dirichlet boundary conditions.

    Parameters
    ----------
    K : (ndof, ndof) ndarray
        Global stiffness matrix.
    R : (ndof, 1) ndarray
        Residual force vector.
    prescribed_dofs : array-like of int
        Constrained degrees of freedom (Dirichlet BCs).

    Returns
    -------
    d : (ndof, 1) ndarray
        Incremental displacement vector.
    R_r : (nfree, 1) ndarray
        Reduced residual vector associated with free DOFs.
    """
    ndof = K.shape[0]

    # Identify free degrees of freedom
    all_dofs = np.arange(ndof)
    free_dofs = np.setdiff1d(all_dofs, prescribed_dofs)

    # Initialize solution vector
    d = np.zeros((ndof, 1))

    # Reduced system
    R_r = R[free_dofs]
    d_free = -np.linalg.solve(K[np.ix_(free_dofs, free_dofs)], R_r)

    # Scatter solution back to full vector
    d[free_dofs, 0] = d_free[:, 0]

    return d, R_r

def compute_C_SVK_component(Gcontra, lambda_, mu, I, J, K, L):
    """
    Compute a single component of the Saint Venant–Kirchhoff elasticity tensor.

    Parameters
    ----------
    Gcontra : (3, 3) ndarray
        Contravariant metric tensor in the reference configuration.
    lambda_ : float
        First Lamé parameter.
    mu : float
        Second Lamé parameter (shear modulus).
    I, J, K, L : int
        Tensor indices (0-based).

    Returns
    -------
    comp : float
        Value of the elasticity tensor component C_{IJKL}.
    """
    comp = (
        lambda_ * Gcontra[I, J] * Gcontra[K, L]
        + mu * (
            Gcontra[I, K] * Gcontra[J, L]
            + Gcontra[I, L] * Gcontra[J, K]
        )
    )
    return comp

def shape_functions_3D(ne_L, xi, eta, zeta):
    """
    Compute shape functions and their derivatives for a 3D brick element
    with 2x2 nodes per cross-section and ne_L nodes along the xi direction.

    Natural coordinates (xi, eta, zeta) are in [-1, 1]^3.

    Parameters
    ----------
    ne_L : int
        Number of nodes along xi (axial direction).
    xi, eta, zeta : float
        Natural coordinates.

    Returns
    -------
    N : (4*ne_L,) ndarray
        Shape functions evaluated at (xi, eta, zeta).
    dN : (3, 4*ne_L) ndarray
        Derivatives with respect to (xi, eta, zeta):
        rows correspond to [d/dxi, d/deta, d/dzeta].
    """
    # -------------------------------------------------------------------------
    # 1D Lagrange shape functions along xi (ne_L nodes)
    # -------------------------------------------------------------------------
    xi_nodes = np.linspace(-1.0, 1.0, ne_L)

    N_xi = np.array([
        np.prod([
            (xi - xi_nodes[m]) / (xi_nodes[a] - xi_nodes[m])
            for m in range(ne_L) if m != a
        ])
        for a in range(ne_L)
    ])

    dN_xi = np.zeros(ne_L)
    for a in range(ne_L):
        for b in range(ne_L):
            if b == a:
                continue
            term = 1.0 / (xi_nodes[a] - xi_nodes[b])
            for m in range(ne_L):
                if m != a and m != b:
                    term *= (xi - xi_nodes[m]) / (xi_nodes[a] - xi_nodes[m])
            dN_xi[a] += term

    # -------------------------------------------------------------------------
    # Bilinear shape functions in (eta, zeta) for a 2x2 cross-section
    # -------------------------------------------------------------------------
    N_eta = np.array([0.5 * (1.0 - eta), 0.5 * (1.0 + eta)])
    dN_eta = np.array([-0.5, 0.5])

    N_zeta = np.array([0.5 * (1.0 - zeta), 0.5 * (1.0 + zeta)])
    dN_zeta = np.array([-0.5, 0.5])

    # -------------------------------------------------------------------------
    # Tensor-product assembly
    # Ordering: xi, eta, zeta
    # -------------------------------------------------------------------------
    ne = 4 * ne_L
    N = np.zeros(ne)
    dN = np.zeros((3, ne))  # rows: d/dxi, d/deta, d/dzeta

    idx = 0
    for a_xi in range(ne_L):
        for a_z in range(2):
            for a_eta in range(2):
                N[idx] = N_eta[a_eta] * N_zeta[a_z] * N_xi[a_xi]

                dN[0, idx] = N_eta[a_eta] * N_zeta[a_z] * dN_xi[a_xi]
                dN[1, idx] = dN_eta[a_eta] * N_zeta[a_z] * N_xi[a_xi]
                dN[2, idx] = N_eta[a_eta] * dN_zeta[a_z] * N_xi[a_xi]

                idx += 1

    return N, dN

def voigt_index(i, j):
    """
    Map a symmetric tensor index (i, j) to its Voigt notation index.

    Voigt mapping for a 3×3 symmetric tensor:
        (0,0) → 0   (xx)
        (1,1) → 1   (yy)
        (2,2) → 2   (zz)
        (0,1),(1,0) → 3   (xy)
        (0,2),(2,0) → 4   (xz)
        (1,2),(2,1) → 5   (yz)

    Parameters
    ----------
    i, j : int
        Tensor indices (0-based).

    Returns
    -------
    int
        Voigt index in {0, …, 5}.
    """
    if i == j:
        return i
    elif (i, j) in ((0, 1), (1, 0)):
        return 3
    elif (i, j) in ((0, 2), (2, 0)):
        return 4
    elif (i, j) in ((1, 2), (2, 1)):
        return 5
    else:
        raise ValueError(f"Invalid tensor indices ({i}, {j}); expected values in {{0,1,2}}.")

def lagrange_shape_1D(ne_L, xi):
    """
    Compute 1D Lagrange shape functions at a given point xi.

    The shape functions are defined on the reference interval [-1, 1]
    using equispaced interpolation nodes.

    Parameters
    ----------
    ne_L : int
        Number of nodes in the xi direction (polynomial order + 1).
    xi : float
        Evaluation point in the reference coordinate [-1, 1].

    Returns
    -------
    N_xi : (ne_L,) ndarray
        Values of the Lagrange shape functions at xi.
    """
    xi_nodes = np.linspace(-1.0, 1.0, ne_L)
    N_xi = np.zeros(ne_L)

    for a in range(ne_L):
        value = 1.0
        for b in range(ne_L):
            if b != a:
                value *= (xi - xi_nodes[b]) / (xi_nodes[a] - xi_nodes[b])
        N_xi[a] = value

    return N_xi