"""Post-processing utilities: tip displacement interpolation, load–displacement plotting, and analytical elastica benchmarks (moment, shear)."""

# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np
import matplotlib.pyplot as plt

def compute_tip_displacement(u, connect, shape_functions_3D, ne_L):
    """
    Interpolate the displacement at the center of the beam tip cross-section.

    Parameters
    ----------
    u : (ndof,) or (ndof, 1) ndarray
        Global displacement vector.
    connect : (numel, ne) ndarray of int
        Element connectivity (0-based global node indices).
    shape_functions_3D : callable
        Shape function evaluator: shape_functions_3D(ne_L, xi, eta, zeta) -> (N, dN).
    ne_L : int
        Number of nodes along the axial direction per element.

    Returns
    -------
    u_tip_center : (3,) ndarray
        Displacement vector interpolated at (xi, eta, zeta) = (1, 0, 0) on the last element.
    """
    # Reshape global displacement vector to nodal form (nnode, 3)
    u_nodes = np.asarray(u).reshape(-1, 3)

    # Last element is assumed to be the tip element
    tip_nodes = connect[-1]
    u_tip_nodes = u_nodes[tip_nodes, :]

    # Tip face center in natural coordinates: xi=1, eta=0, zeta=0
    N, _ = shape_functions_3D(ne_L, 1.0, 0.0, 0.0)

    # Interpolate displacement at the tip center
    u_tip_center = N @ u_tip_nodes
    return u_tip_center

def plot_tip_displacement(tip_displacements, length, load_type, k_max):
    """
    Plot the normalized tip displacement versus the normalized load parameter.
    """
    if load_type not in ("moment", "shear"):
        raise ValueError("load_type must be either 'moment' or 'shear'.")

    # Extract the displacement component relevant to the loading type
    disp_vals = []
    for d in tip_displacements:
        a = np.asarray(d)

        if a.ndim == 0:
            val = float(a)
        else:
            flat = a.reshape(-1)
            val = float(flat[0] if load_type == "moment" else flat[2])

        disp_vals.append(abs(val))

    # Convert to array and prepend undeformed state
    disp = np.asarray(disp_vals, dtype=float)
    disp = np.insert(disp, 0, 0.0)

    # Normalize displacement by beam length
    disp_norm = disp / length
    n_steps = disp_norm.size

    plt.figure(figsize=(6, 4))

    # --------------------------------------------------
    # MOMENT LOADING
    # --------------------------------------------------
    if load_type == "moment":
        k = np.linspace(0.0, k_max, n_steps)
        analytical = analytical_u_over_L_for_moment(k)

        plt.plot(k, disp_norm, marker="o", linestyle="-",
                 color="blue", label="FEM solution")
        plt.plot(k, analytical, linestyle="-",
                 color="black", label="Analytical solution (elastica)")

        plt.xlabel(r"$k = \dfrac{ML}{2\pi EI}$", fontsize=12)
        plt.ylabel(r"$\dfrac{u}{L}$", fontsize=12)
        plt.xlim(0.0, k_max)
        plt.ylim(0.0, 1.2)

    # --------------------------------------------------
    # SHEAR LOADING
    # --------------------------------------------------
    else:
        k = np.linspace(0.0, 4.0, n_steps)
        analytical = analytical_w_over_L_for_shear(k)

        plt.plot(k, disp_norm, marker="o", linestyle="-",
                 color="blue", label="FEM solution")
        plt.plot(k, analytical, linestyle="-",
                 color="black", label="Analytical solution (elastica)")

        plt.xlabel(r"$k = \dfrac{PL^2}{EI}$", fontsize=12)
        plt.ylabel(r"$\dfrac{w}{L}$", fontsize=12)
        plt.xlim(0.0, 4.0)
        plt.ylim(0.0, 1.2)

    # Final plot formatting
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def analytical_u_over_L_for_moment(k: np.ndarray) -> np.ndarray:
    """
    Analytical elastica solution for a cantilever beam under an end moment.

    The normalized tip displacement is
        u / L = 1 - sin(2πk) / (2πk),
    where k = ML / (2πEI).

    """
    k = np.asarray(k, dtype=float)

    return 1.0 - np.sinc(2.0 * k)

def analytical_w_over_L_for_shear(k_vals, n_int=4000, tol=1e-10, max_iter=80):
    """
    Cantilever elastica under an end shear load.

    The governing elastica boundary-value problem is solved numerically using a
    shooting method. The unknown initial slope is determined by bisection such
    that the free-end condition is satisfied, while the resulting initial-value
    problem is integrated with a fourth-order Runge–Kutta (RK4) scheme in the
    normalized arc-length s̄ = s/L.

    Returns
    -------
    w_over_L : ndarray
        Normalized vertical tip deflection w/L = |y(1)| for each k, where
        y'(s̄) = sin(theta) and k = P L^2 / (E I).
    """
    k_vals = np.asarray(k_vals, dtype=float)
    w_out = np.zeros_like(k_vals)

    # Normalized arc-length step size on s̄ ∈ [0, 1]
    h = 1.0 / n_int

    def integrate(k, alpha):
        """
        Integrate the elastica ODE system on s̄ ∈ [0, 1] using RK4.

        The first-order system is
            theta' = p
            p'     = k cos(theta)
            y'     = sin(theta)
        with initial conditions theta(0)=0, p(0)=alpha, y(0)=0.

        Returns
        -------
        p_end : float
            Value p(1) used by the shooting residual.
        y_end : float
            Value y(1) giving the normalized tip deflection.
        """
        theta = 0.0
        p = alpha
        y = 0.0

        def rhs(state):
            th, pp, yy = state
            return np.array([pp, k * np.cos(th), np.sin(th)], dtype=float)

        for _ in range(n_int):
            s0 = np.array([theta, p, y], dtype=float)

            k1 = rhs(s0)
            k2 = rhs(s0 + 0.5 * h * k1)
            k3 = rhs(s0 + 0.5 * h * k2)
            k4 = rhs(s0 + h * k3)

            s1 = s0 + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            theta, p, y = s1

        return p, y

    for i, k in enumerate(k_vals):
        if k == 0.0:
            w_out[i] = 0.0
            continue

        # Shooting unknown: alpha = p(0). Target condition: p(1) = 0.
        a_lo, a_hi = -50.0, 0.0
        f_lo = integrate(k, a_lo)[0]
        f_hi = integrate(k, a_hi)[0]

        # Bracket the root by expanding the lower bound until a sign change occurs.
        expand = 0
        while f_lo * f_hi > 0.0 and expand < 20:
            a_lo *= 2.0
            f_lo = integrate(k, a_lo)[0]
            expand += 1

        if f_lo * f_hi > 0.0:
            raise RuntimeError(
                f"Could not bracket the shooting parameter for k={k}. "
                "Increase the bracket range or check the sign convention."
            )

        # Bisection on alpha for the scalar shooting residual f(alpha) = p(1).
        for _ in range(max_iter):
            a_mid = 0.5 * (a_lo + a_hi)
            f_mid, y_mid = integrate(k, a_mid)

            if abs(f_mid) < tol:
                w_out[i] = abs(y_mid)  # w/L = |y(1)|
                break

            if f_lo * f_mid <= 0.0:
                a_hi = a_mid
                f_hi = f_mid
            else:
                a_lo = a_mid
                f_lo = f_mid
        else:
            w_out[i] = abs(y_mid)

    return w_out