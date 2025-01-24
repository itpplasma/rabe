# %%
import matplotlib.pyplot as plt
from rabe.boozer_modes import read_modes_bc


def compare_coefs(modes_par, modes_ql, ms, ns):
    _, ax = plt.subplots(nrows=len(ns), ncols=len(ms), figsize=(20, 12))
    for idx_m in range(len(ms)):
        for idx_n in range(len(ns)):
            m = ms[idx_m]
            n = ns[idx_n]
            mode_par = modes_par.get_mode_idx(m, -n)
            mode_ql = modes_ql.get_mode_idx(m, n)
            if m == 0 and n < 0:
                continue
            elif m == 0 and n >= 0:
                mode_par = mode_ql
            modenumber_par = (
                " ("
                + str(modes_par.m[0][mode_par])
                + ","
                + str(modes_par.n[0][mode_par])
                + ")"
            )
            modenumber_ql = (
                " ("
                + str(modes_ql.m[0][mode_ql])
                + ","
                + str(modes_ql.n[0][mode_ql])
                + ")"
            )
            modenumbers = "par:" + modenumber_par + " / ql:" + modenumber_ql
            ax[idx_n, idx_m].plot(
                modes_par.rho_tor, modes_par.coefs[:, mode_par], "b-", label=modenumbers
            )
            ax[idx_n, idx_m].plot(modes_ql.rho_tor, modes_ql.coefs[:, mode_ql], "b--")
            ax[idx_n, idx_m].set_xlabel(r"$\rho_\mathrm{tor}$")
            ax[idx_n, idx_m].set_ylabel(r"coef")
            ax[idx_n, idx_m].legend()


if __name__ == "__main__":
    import sys
    import os

    output_dir = sys.argv[1]

    ms = [0, 1, 3]
    ns = [-3, -1, 0, 1, 3]

    bc_file_par = os.path.join(output_dir, "booz_xform_field.bc")
    get_mode_idx_par = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max + n) + (n_max + 1)
    )

    bc_file_ql = os.path.join(output_dir, "field.bc")
    get_mode_idx_ql = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max + n) + (n_max + 1)
    )

    rmnc_par, zmns_par, vmns_par, bmnc_par = read_modes_bc(
        bc_file_par, get_mode_idx_par
    )
    rmnc_ql, zmns_ql, vmns_ql, bmnc_ql = read_modes_bc(bc_file_ql, get_mode_idx_ql)
    compare_coefs(rmnc_par, rmnc_ql, ms, ns)
    compare_coefs(zmns_par, zmns_ql, ms, ns)
    compare_coefs(vmns_par, vmns_ql, ms, ns)
    compare_coefs(bmnc_par, bmnc_ql, ms, ns)
    plt.show()
