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


def get_xyz_surface(
    rmnc: Modes,
    zmns: Modes,
    vmns: Modes,
    bmnc: Modes,
    nfp: int,
    n_theta: int = 100,
    n_phi: int = 100,
    idx_surface: int = -1,
):
    from .fourier_series import FourierSeries, evaluate

    fourier_r = FourierSeries(
        list(rmnc.m[idx_surface]),
        list(nfp * rmnc.n[idx_surface]),
        list(rmnc.coefs[idx_surface]),
        list(np.zeros(len(rmnc.m[idx_surface]))),
    )
    fourier_z = FourierSeries(
        list(zmns.m[idx_surface].tolist()),
        list(nfp * zmns.n[idx_surface]),
        list(np.zeros(len(zmns.m[idx_surface]))),
        list(zmns.coefs[idx_surface]),
    )
    fourier_v = FourierSeries(
        list(vmns.m[idx_surface]),
        list(nfp * vmns.n[idx_surface]),
        list(np.zeros(len(vmns.m[idx_surface]))),
        list(vmns.coefs[idx_surface]),
    )

    theta, phi_boozer, surface_B = get_theta_phi_surface(
        bmnc=bmnc, nfp=nfp, n_theta=n_theta, n_phi=n_phi, idx_surface=idx_surface
    )

    surface_R = evaluate(fourier_r, theta, phi_boozer)
    surface_z = evaluate(fourier_z, theta, phi_boozer)
    phi = phi_boozer + (2 * np.pi) / nfp * evaluate(fourier_v, theta, phi_boozer)
    # The coefficients are not given in terms of phi, but phi_boozer. However,
    # to get the cartesian coordinates from the simple sin/cos relation with R,
    # phi has to be used.
    #
    # The calculation of R, z is done using phi_boozer and then we take these phi_boozer
    # values and transform them into phi values using the (back)transformation
    #
    # phi_boozer + p(s,theta_boozer,phi_boozer) = phi
    #
    # where
    #
    # p = 2*pi/nfp * v
    #
    # with v, the normalized transformation function

    phi_fullperiod = phi.copy()
    for i in range(nfp - 1):
        phi_fullperiod = np.concatenate(
            (phi_fullperiod, phi + 2 * np.pi * (i + 1) / nfp), axis=0
        )

    surface_R = np.tile(surface_R, (nfp, 1))
    surface_z = np.tile(surface_z, (nfp, 1))
    surface_B = np.tile(surface_B, (nfp, 1))

    surface_x = surface_R * np.cos(phi_fullperiod)
    surface_y = surface_R * np.sin(phi_fullperiod)

    return surface_x, surface_y, surface_z, surface_B


def get_theta_phi_surface(
    bmnc: Modes,
    nfp: int,
    n_theta: int = 100,
    n_phi: int = 100,
    idx_surface: int = -1,
):
    from .fourier_series import FourierSeries, evaluate

    fourier_b = FourierSeries(
        list(bmnc.m[idx_surface]),
        list(nfp * bmnc.n[idx_surface]),
        list(bmnc.coefs[idx_surface]),
        list(np.zeros(len(bmnc.m[idx_surface]))),
    )

    theta = np.linspace(0.0, 2 * np.pi, n_theta)
    phi_boozer = np.linspace(0.0, 2 * np.pi / nfp, n_phi)
    theta, phi_boozer = np.meshgrid(theta, phi_boozer)

    surface_B = evaluate(fourier_b, theta, phi_boozer)
    # The coefficients are not given in terms of phi, but phi_boozer on the
    # interval [0, 2pi/nfp].

    return theta, phi_boozer, surface_B


def get_axis_projection(rmnc: Modes, zmns: Modes, nfp: int, n_phi: int = 100):
    from .fourier_series import FourierSeries, evaluate

    idx_surface = 0
    fourier_r = FourierSeries(
        list(rmnc.m[idx_surface]),
        list(nfp * rmnc.n[idx_surface]),
        list(rmnc.coefs[idx_surface]),
        list(np.zeros(len(rmnc.m[idx_surface]))),
    )
    fourier_z = FourierSeries(
        list(zmns.m[idx_surface]),
        list(nfp * zmns.n[idx_surface]),
        list(np.zeros(len(zmns.m[idx_surface]))),
        list(zmns.coefs[idx_surface]),
    )

    phi_boozer = np.linspace(0, 2 * np.pi / nfp / 2, n_phi)
    theta = np.zeros_like(phi_boozer)
    R = evaluate(fourier_r, theta, phi_boozer)
    Z = evaluate(fourier_z, theta, phi_boozer)

    return phi_boozer, R, Z
