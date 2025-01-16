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


def split_off_symmetric_modes(modes, helicity_n):
    symmetric_modes = copy.deepcopy(modes)
    non_symmetric_modes = copy.deepcopy(modes)
    for m in range(modes.m_max + 1):
        for n in range(-modes.n_max, modes.n_max + 1):
            mode_idx = modes.get_mode_idx(m, n)
            if mode_idx < 0:
                continue
            if m * helicity_n == n:
                non_symmetric_modes.coefs[:, mode_idx] = 0.0
            else:
                symmetric_modes.coefs[:, mode_idx] = 0.0
    return symmetric_modes, non_symmetric_modes
