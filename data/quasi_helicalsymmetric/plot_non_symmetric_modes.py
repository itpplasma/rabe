# %%
if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    from rabe.boozer_modes import read_modes_bc, split_off_symmetric_modes

    nfp = 4
    par_bc_file = "output/all_modes.bc"
    get_mode_idx_par = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max + n) + (n_max + 1)
    )
    _, _, _, bmnc = read_modes_bc(par_bc_file, get_mode_idx=get_mode_idx_par)

    symmetric_modes, non_symmetric_modes = split_off_symmetric_modes(bmnc, -1)

    stor = non_symmetric_modes.rho_tor**2
    max_B_nonsym = np.max(np.abs(non_symmetric_modes.coefs), axis=1)
    B_00_sym = symmetric_modes.coefs[:, 0]

    plt.figure()
    plt.title("maximal symmetry-breaking mode")
    plt.plot(stor, max_B_nonsym / B_00_sym)
    plt.ylabel(r"$\max{(B^\mathrm{nonsym}_{mn})}/B^\mathrm{sym}_{00}$")
    plt.xlabel(r"$s_\mathrm{tor}$")
    plt.show()
