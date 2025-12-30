"""Visualization utilities: Visualize undeformed and final configurations"""

import numpy as np
import matplotlib.pyplot as plt

def plot_mesh_3D(
    coords: np.ndarray,
    connect: np.ndarray,
    prescribed_nodes: np.ndarray | None = None,
    show_node_numbers: bool = False,
    elev: float = 20.0,
    azim: float = -60.0,
):
    """
    Plot the 3D reference mesh.

    Plotting convention
    -------------------
    - Axial lines: connect slice s → s+1 for all slices in each element.
    - Transverse lines: draw only on the first and last axial slice of each element.

    Notes
    -----
    This routine assumes each axial slice contains 4 nodes and that each element
    connectivity is ordered by axial slices:
        [slice0(4 nodes), slice1(4 nodes), ..., sliceEnd(4 nodes)].
    """
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    connect = np.asarray(connect, dtype=int)

    fig = plt.figure(figsize=(8, 6), dpi=150)
    ax = fig.add_subplot(111, projection="3d")

    xs, ys, zs = coords[:, 0], coords[:, 1], coords[:, 2]

    # Use unique coordinate values as ticks
    ax.set_xticks(np.unique(np.round(xs, decimals=6)))
    ax.set_yticks(np.unique(np.round(ys, decimals=6)))
    ax.set_zticks(np.unique(np.round(zs, decimals=6)))

    if prescribed_nodes is None:
        prescribed_nodes = np.array([], dtype=int)
    else:
        prescribed_nodes = np.asarray(prescribed_nodes, dtype=int)

    all_nodes = np.arange(coords.shape[0], dtype=int)
    free_nodes = np.setdiff1d(all_nodes, prescribed_nodes)

    # Nodes
    ax.scatter(xs[free_nodes], ys[free_nodes], zs[free_nodes],
               c="blue", s=10, depthshade=True, label="Free nodes")
    if prescribed_nodes.size > 0:
        ax.scatter(xs[prescribed_nodes], ys[prescribed_nodes], zs[prescribed_nodes],
                   c="red", s=20, depthshade=True, label="Clamped nodes")

    n_per_slice = 4

    for element in connect:
        element = np.asarray(element, dtype=int)

        if element.size % n_per_slice != 0:
            raise ValueError(
                f"Connectivity length {element.size} is not a multiple of {n_per_slice}. "
                "This plotting routine assumes 4 nodes per axial slice."
            )

        n_slices = element.size // n_per_slice
        if n_slices < 2:
            continue

        # Transverse connections (only at first and last slice)
        for s in (0, n_slices - 1):
            i = s * n_per_slice
            n0, n1, n2, n3 = element[i:i + 4]

            ax.plot(*zip(coords[n0], coords[n1]), color="gray", linewidth=1)
            ax.plot(*zip(coords[n1], coords[n3]), color="gray", linewidth=1)
            ax.plot(*zip(coords[n3], coords[n2]), color="gray", linewidth=1)
            ax.plot(*zip(coords[n2], coords[n0]), color="gray", linewidth=1)

        # Axial connections (between every consecutive slice)
        for s in range(n_slices - 1):
            i = s * n_per_slice
            j = (s + 1) * n_per_slice

            n0, n1, n2, n3 = element[i:i + 4]
            m0, m1, m2, m3 = element[j:j + 4]

            ax.plot(*zip(coords[n0], coords[m0]), color="gray", linewidth=1)
            ax.plot(*zip(coords[n1], coords[m1]), color="gray", linewidth=1)
            ax.plot(*zip(coords[n2], coords[m2]), color="gray", linewidth=1)
            ax.plot(*zip(coords[n3], coords[m3]), color="gray", linewidth=1)

    if show_node_numbers:
        for node_id, (x, y, z) in enumerate(coords):
            ax.text(x, y, z, str(node_id), fontsize=8, color="black",
                    verticalalignment="bottom", horizontalalignment="right")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=elev, azim=azim)

    # Preserve aspect ratio in data units
    xrange = np.ptp(xs) or 1.0
    yrange = np.ptp(ys) or 1.0
    zrange = np.ptp(zs) or 1.0
    ax.set_box_aspect([xrange, yrange, zrange])

    ax.grid(True)
    ax.set_title("Reference configuration", fontsize=14, fontweight="bold")
    ax.legend()

    plt.tight_layout()
    plt.show()
    plt.close(fig)
    
def plot_deformed_mesh_3D(
    coords: np.ndarray,
    connect: np.ndarray,
    u: np.ndarray,
    scale: float = 1.0,
    prescribed_nodes: np.ndarray | None = None,
    show_node_numbers: bool = False,
    elev: float = 20.0,
    azim: float = -60.0,
):
    """
    Plot the deformed 3D mesh.

    Plotting convention
    -------------------
    - Transverse lines: only on the first and last axial slice of each element.
    - Axial lines: between every consecutive slice.
    """
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    connect = np.asarray(connect, dtype=int)
    u = np.asarray(u, dtype=float).reshape(-1, 3)

    coords_def = coords + scale * u

    fig = plt.figure(figsize=(8, 6), dpi=150)
    ax = fig.add_subplot(111, projection="3d")

    xs, ys, zs = coords_def[:, 0], coords_def[:, 1], coords_def[:, 2]

    # Use unique coordinate values as ticks (useful for structured meshes)
    ax.set_xticks(np.unique(np.round(xs, decimals=6)))
    ax.set_yticks(np.unique(np.round(ys, decimals=6)))
    ax.set_zticks(np.unique(np.round(zs, decimals=6)))

    if prescribed_nodes is None:
        prescribed_nodes = np.array([], dtype=int)
    else:
        prescribed_nodes = np.asarray(prescribed_nodes, dtype=int)

    all_nodes = np.arange(coords_def.shape[0], dtype=int)
    free_nodes = np.setdiff1d(all_nodes, prescribed_nodes)

    # Nodes
    ax.scatter(xs[free_nodes], ys[free_nodes], zs[free_nodes],
               c="blue", s=10, depthshade=True, label="Free nodes")
    if prescribed_nodes.size > 0:
        ax.scatter(xs[prescribed_nodes], ys[prescribed_nodes], zs[prescribed_nodes],
                   c="red", s=20, depthshade=True, label="Clamped nodes")

    n_per_slice = 4

    for element in connect:
        element = np.asarray(element, dtype=int)

        if element.size % n_per_slice != 0:
            raise ValueError(
                f"Connectivity length {element.size} is not a multiple of {n_per_slice}. "
                "This plotting routine assumes 4 nodes per axial slice."
            )

        n_slices = element.size // n_per_slice
        if n_slices < 2:
            continue

        # Transverse connections (only at first and last slice)
        for s in (0, n_slices - 1):
            i = s * n_per_slice
            n0, n1, n2, n3 = element[i:i + 4]

            ax.plot(*zip(coords_def[n0], coords_def[n1]), color="gray", linewidth=1)
            ax.plot(*zip(coords_def[n1], coords_def[n3]), color="gray", linewidth=1)
            ax.plot(*zip(coords_def[n3], coords_def[n2]), color="gray", linewidth=1)
            ax.plot(*zip(coords_def[n2], coords_def[n0]), color="gray", linewidth=1)

        # Axial connections (between every consecutive slice)
        for s in range(n_slices - 1):
            i = s * n_per_slice
            j = (s + 1) * n_per_slice

            n0, n1, n2, n3 = element[i:i + 4]
            m0, m1, m2, m3 = element[j:j + 4]

            ax.plot(*zip(coords_def[n0], coords_def[m0]), color="gray", linewidth=1)
            ax.plot(*zip(coords_def[n1], coords_def[m1]), color="gray", linewidth=1)
            ax.plot(*zip(coords_def[n2], coords_def[m2]), color="gray", linewidth=1)
            ax.plot(*zip(coords_def[n3], coords_def[m3]), color="gray", linewidth=1)

    if show_node_numbers:
        for node_id, (x, y, z) in enumerate(coords_def):
            ax.text(x, y, z, str(node_id), fontsize=8, color="black",
                    verticalalignment="bottom", horizontalalignment="right")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=elev, azim=azim)

    # Preserve aspect ratio in data units
    xrange = np.ptp(xs) or 1.0
    yrange = np.ptp(ys) or 1.0
    zrange = np.ptp(zs) or 1.0
    ax.set_box_aspect([xrange, yrange, zrange])

    ax.grid(True)
    ax.set_title("Deformed configuration", fontsize=14, fontweight="bold")
    ax.legend()

    plt.tight_layout()
    plt.show()
    plt.close(fig)