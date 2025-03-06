# %%
if __name__ == "__main__":
    from simsopt.mhd.vmec import Vmec
    from rabe.boozer_modes import read_modes_bc
    from rabe.coordinate_orientations import get_bc_theta_orientation

    par_bc_file = "python/output/quasi_symmetric.bc"
    get_mode_idx_par = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max + n) + (n_max + 1)
    )
    rmnc, zmns, vmns, bmnc = read_modes_bc(par_bc_file, get_mode_idx=get_mode_idx_par)
    print("For Landreman bc file")
    get_bc_theta_orientation(rmnc, zmns)

    from rabe.boozer_modes import (
        get_xyz_surface,
        get_theta_phi_surface,
        get_axis_projection,
    )
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    x, y, z, B = get_xyz_surface(rmnc, zmns, vmns, bmnc, 4)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    norm = plt.Normalize(B.min(), B.max())
    colors = plt.cm.viridis(norm(B))
    surface = ax.plot_surface(x, y, z, facecolors=colors, edgecolor="none")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.axis("equal")
    mappable = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    mappable.set_array(B)
    fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10)

    theta, phi_boozer, B = get_theta_phi_surface(bmnc, 4)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    norm = plt.Normalize(B.min(), B.max())
    colors = plt.cm.viridis(norm(B))
    surface = ax.contour(phi_boozer, theta, B)
    ax.set_ylabel(r"$\theta_\mathrm{B}$")
    ax.set_xlabel(r"$\varphi_\mathrm{B}$")
    mappable = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    mappable.set_array(B)
    fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10)

    phi_boozer, R, z = get_axis_projection(rmnc, zmns, nfp=4)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.plot(R, z, "r.", label="approx. magnetic axis")
    ax.set_xlabel("R")
    ax.set_ylabel("z")
    plt.axis("equal")
    plt.legend()
    plt.show()
