# %%
import matplotlib.pyplot as plt
import numpy as np
import copy


class Modes:
    def init(self):
        self.rho_tor = np.array([])
        self.coefs = np.array([[]])
        self.m = np.array([[]])
        self.n = np.array([[]])
        self.m_max = 0
        self.n_max = 0
        self.get_mode_idx = lambda m, n: 0


def read_modes_bc(bc_filename: str, get_mode_idx):
    from libneo import BoozerFile

    rmnc = Modes()
    bc_file = BoozerFile(bc_filename)
    rmnc.rho_tor = np.sqrt(np.array(bc_file.s))
    rmnc.m = np.array(bc_file.m)
    rmnc.n = np.array(bc_file.n)
    rmnc.m_max = np.max(rmnc.m[0])
    rmnc.n_max = np.max(rmnc.n[0])
    rmnc.get_mode_idx = lambda m, n: get_mode_idx(m, n, n_max=rmnc.n_max)
    zmns = copy.deepcopy(rmnc)
    vmns = copy.deepcopy(rmnc)
    bmnc = copy.deepcopy(rmnc)
    rmnc.coefs = np.array(bc_file.rmnc)
    zmns.coefs = np.array(bc_file.zmns)
    vmns.coefs = np.array(bc_file.vmns)
    bmnc.coefs = np.array(bc_file.bmnc)
    return rmnc, zmns, vmns, bmnc


def compare_coefs(modes_1, modes_2, ms, ns):
    _, ax = plt.subplots(nrows=len(ns), ncols=len(ms), figsize=(20, 12))
    for idx_m in range(len(ms)):
        for idx_n in range(len(ns)):
            m = ms[idx_m]
            n = ns[idx_n]
            mode_idx1 = modes_1.get_mode_idx(m, n)
            mode_idx2 = modes_2.get_mode_idx(m, n)
            if m == 0:
                mode_idx2 = mode_idx1
            modenumber_1 = (
                " ("
                + str(modes_1.m[0][mode_idx1])
                + ","
                + str(modes_1.n[0][mode_idx1])
                + ")"
            )
            modenumber_2 = (
                " ("
                + str(modes_2.m[0][mode_idx2])
                + ","
                + str(modes_2.n[0][mode_idx2])
                + ")"
            )
            modenumbers = modenumber_1 + " / " + modenumber_2
            ax[idx_n, idx_m].plot(
                modes_1.rho_tor, modes_1.coefs[:, mode_idx1], "b-", label=modenumbers
            )
            ax[idx_n, idx_m].plot(modes_2.rho_tor, modes_2.coefs[:, mode_idx2], "b--")
            ax[idx_n, idx_m].set_xlabel(r"$\rho_\mathrm{tor}$")
            ax[idx_n, idx_m].set_ylabel(r"coef")
            ax[idx_n, idx_m].legend()


if __name__ == "__main__":
    ms = [0, 1, 3]
    ns = [0, 1, 3]
    bc_file_1 = "field.bc"
    get_mode_idx_1 = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max + n) + (n_max + 1)
    )
    bc_file_2 = "booz_xform_field.bc"
    get_mode_idx_2 = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max - n) + (n_max + 1)
    )
    rmnc_1, zmns_1, vmns_1, bmnc_1 = read_modes_bc(bc_file_1, get_mode_idx_1)
    rmnc_2, zmns_2, vmns_2, bmnc_2 = read_modes_bc(bc_file_2, get_mode_idx_2)
    compare_coefs(rmnc_1, rmnc_2, ms, ns)
    compare_coefs(zmns_1, zmns_2, ms, ns)
    compare_coefs(vmns_1, vmns_2, ms, ns)
    compare_coefs(bmnc_1, bmnc_2, ms, ns)
    plt.show()
