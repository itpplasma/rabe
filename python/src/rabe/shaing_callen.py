import numpy as np
from rabe.fourier_series import FourierSeries, evaluate


def shaing_callen_bootstrap(
    bc_filename: str,
    stor: np.double,
    iota: np.double,
    b_0: np.double,
    avnabAphi: np.double,
    N: int = 0,
):
    """
    all input quantities in SI units

    avnabAphi ... the average <|nabla Aphi|> = <|nabla psi|> with sign(Aphi)
    """
    from scipy.integrate import quad

    bmnc_spl, bmns_spl, m, n = get_bmn_splines(bc_filename)
    Bphi_spl = get_Bphi_spline(bc_filename)
    Btheta_spl = get_Btheta_spline(bc_filename)
    bmnc = bmnc_spl(stor).tolist()
    bmns = bmns_spl(stor).tolist()
    fourier_b = FourierSeries(m, n, bmnc, bmns)

    def integrand(x):
        def f(alpha):
            # for alpha = M*theta - N*phi for M=1
            # -> m = m_alpha * M = m_alpha
            # -> n = m_alpha * N
            # -> cos(m*theta - n*phi) = cos(m*alpha)
            # can therefore evaluate as if alpha = theta and phi=0
            b = evaluate(fourier_b, theta=alpha, phi=0)
            return np.sqrt(1 - x * b / b_0) / (b / b_0) ** 2

        alpha_average, _ = quad(f, 0, 2 * np.pi)
        alpha_average = alpha_average / (2 * np.pi)
        return x / alpha_average

    Bphi = Bphi_spl(stor)
    Btheta = Btheta_spl(stor)
    integral, _ = quad(integrand, 0, 1)
    avnabAphi_cgs = avnabAphi * 1e8
    lambda_bB_analytic = (
        (3 / 4 * integral - 1) * (N * Btheta + Bphi) / (avnabAphi_cgs * (N - iota))
    )

    return lambda_bB_analytic


def get_bmn_splines(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    m = bc_file.m[0]
    n = (np.array(bc_file.n[0]) * np.array(bc_file.nper)).tolist()

    bmnc_spl = interp1d(stor, np.array(bc_file.bmnc), axis=0)
    bmns_spl = interp1d(stor, np.array(bc_file.bmns), axis=0)
    return bmnc_spl, bmns_spl, m, n


def get_Bphi_spline(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    c_cgs = 3e10
    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    Jpol_SI_lhs = np.array(bc_file.Jpol_divided_by_nper) * np.array(bc_file.nper)
    Jpol_cgs_lhs = Jpol_SI_lhs * 3e9
    Bphi_cgs = -2 * Jpol_cgs_lhs / c_cgs

    Bphi_spl = interp1d(stor, Bphi_cgs, axis=0)
    return Bphi_spl


def get_Btheta_spline(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    c_cgs = 3e10
    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    Itor_SI_lhs = np.array(bc_file.Itor)
    Itor_cgs_lhs = Itor_SI_lhs * 3e9
    Btheta_cgs = -2 * Itor_cgs_lhs / c_cgs

    Btheta_spl = interp1d(stor, Btheta_cgs, axis=0)
    return Btheta_spl
