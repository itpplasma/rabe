import numpy as np
from rabe.fourier_series import FourierSeries, evaluate


def shaing_callen_bootstrap(
    bc_filename: str,
    stor: np.double,
    iota: np.double,
    b_0: np.double,
    avnabpsi: np.double,
):
    """
    all input quantities in SI units
    """
    from scipy.integrate import quad

    bmnc_spl, bmns_spl, m, n = get_bmn_splines(bc_filename)
    Bphi_spl = get_Bphi_spline(bc_filename)
    bmnc = bmnc_spl(stor).tolist()
    bmns = bmns_spl(stor).tolist()
    fourier_b = FourierSeries(m, n, bmnc, bmns)

    def integrand(x):
        def f(theta):
            b = evaluate(fourier_b, theta, phi=0)
            return np.sqrt(1 - x * b / b_0) / (b / b_0) ** 2

        theta_average, _ = quad(f, 0, 2 * np.pi)
        theta_average = theta_average / (2 * np.pi)
        return x / theta_average

    Bphi = Bphi_spl(stor)
    integral, _ = quad(integrand, 0, 1)
    avnabpsi_cgs = avnabpsi * 1e8
    lambda_bB_analytic = -(3 / 4 * integral - 1) * Bphi / (np.abs(avnabpsi_cgs) * iota)

    return lambda_bB_analytic


def get_bmn_splines(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    m = bc_file.m[0]
    n = bc_file.n[0]

    bmnc_spl = interp1d(stor, np.array(bc_file.bmnc), axis=0)
    bmns_spl = interp1d(stor, np.array(bc_file.bmns), axis=0)
    return bmnc_spl, bmns_spl, m, n


def get_Bphi_spline(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    c_cgs = 3e10
    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    Jpol_SI = np.array(bc_file.Jpol_divided_by_nper) * np.array(bc_file.nper)
    Jpol_cgs = Jpol_SI * 3e9
    Bphi_cgs = -2 * Jpol_cgs / c_cgs

    Bphi_spl = interp1d(stor, Bphi_cgs, axis=0)
    return Bphi_spl
