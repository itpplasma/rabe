import numpy as np
import matplotlib.pyplot as plt

from rabe.fourier_field import FourierField
from rabe.fieldline_mod import FlockOfFieldlines


# ---------------------------------------------------------------------------
# Analytical fit functions
# ---------------------------------------------------------------------------


def S_A(iota_p):
    return 0.3 * (iota_p - 0.5 * np.pi * np.sign(iota_p))


def S_B(iota_p):
    return np.sign(iota_p) * 2.0


def calc_lambda_a_analytic(
    ds_datheta, R, B_0, eps_1, eps_2, N_TOR, M_POL, iota_p, iota
):
    """Small-aspect-ratio analytical estimate for lambda_A."""
    result = -ds_datheta * np.sqrt(2.0 * np.pi) * R * B_0 * np.abs(eps_2) ** 2
    result *= N_TOR * S_A(iota_p) / (2.0 * np.abs(eps_1)) ** 0.75 * 0.5
    result /= (N_TOR - iota * M_POL) ** 1.5
    result *= 1.0 + 16.0 / 3.0 * np.abs(eps_1)
    return result


def calc_lambda_b_analytic(ds_datheta, R, eps_2, delta_B_1, iota_p):
    """Small-aspect-ratio analytical estimate for lambda_B."""
    return -abs(eps_2) * ds_datheta * R * np.pi * delta_B_1 * S_B(iota_p) * 0.25


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def plot_field(field, nfp):
    theta = np.linspace(0, 2.0 * np.pi, 100)
    phi = np.linspace(0, 2.0 * np.pi / nfp, 100)
    Theta, Phi = np.meshgrid(theta, phi)
    B_array = np.zeros_like(Theta)
    for k in range(len(theta)):
        for l in range(len(phi)):
            B_array[k, l] = field.compute_b_mod(Theta[k, l], Phi[k, l])
    plt.contourf(Phi, Theta, B_array, levels=50)
    plt.xlabel(r"$\varphi$")
    plt.ylabel(r"$\vartheta$")
    plt.title(r"$B(\vartheta,\varphi)$")
    plt.colorbar(label="B")
    plt.show()


def plot_offset_vs_nu_star(lam_a, lam_b, lam_a_ana=None, lam_b_ana=None):
    """Log-log plot of offset coefficients vs collisionality."""
    nu_star = 10.0 ** np.linspace(0.0, -8.0, 200)

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(
        nu_star,
        np.abs(lam_a) / np.sqrt(nu_star),
        "r-",
        label=f"$\\lambda_A$ (factor={lam_a:.3e})",
    )
    if lam_a_ana is not None:
        rel = abs(lam_a_ana / lam_a - 1.0)
        ax.loglog(
            nu_star,
            np.abs(lam_a_ana) / np.sqrt(nu_star),
            "ro",
            markersize=4,
            markevery=10,
            label=f"$\\lambda_A$ analytic (rel.diff.={rel:.3f})",
        )

    ax.loglog(
        nu_star,
        np.abs(lam_b) / nu_star,
        "b-",
        label=f"$\\lambda_B$ (factor={lam_b:.3e})",
    )
    if lam_b_ana is not None and lam_b != 0:
        rel_b = abs(lam_b_ana / lam_b - 1.0)
        ax.loglog(
            nu_star,
            np.abs(lam_b_ana) / nu_star,
            "bo",
            markersize=4,
            markevery=10,
            label=f"$\\lambda_B$ analytic (rel.diff.={rel_b:.3f})",
        )

    total = np.abs(lam_a / np.sqrt(nu_star) + lam_b / nu_star)
    ax.loglog(nu_star, total, "c--", label="total")

    ax.set_xlabel(r"$\nu_*$")
    ax.set_ylabel(r"$|\lambda_\mathrm{off}|$")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

PLOT_FIELD = False
EPS_1 = 0.005
EPS_2 = 0.00009375
IOTA = 0.47
NFP = 10
B_0 = 1.0  # T
B_PHI_COV = 8.0  # Tm
B_THETA_COV = 0.0  # Tm
R = 8.0  # m
PSI_TOR_EDGE = 2.812e-07  # Tm^2

m = np.array([0, 1, 1, 2, 3], dtype=np.int32)
n = np.array([0, 0, 1, 1, 1], dtype=np.int32)  # normalized to NFP
B_mn = B_0 * np.array(
    [1.0, -EPS_2, -EPS_2 * 0.5, -EPS_1, -EPS_2 * 0.5], dtype=np.float64
)
delta_B_1 = B_0 * EPS_2 * 0.01
B_mn[1] += delta_B_1

field = FourierField()
field.fourier_field_init(m, n, B_mn, B_THETA_COV, B_PHI_COV, nfp=NFP)

# M_POL is the mode of the maximal B_mn mode where n=NFP (1 if normalized to NFP)
M_POL = m[np.argmax(np.abs(B_mn * (n == 1)))]
N_TOR = float(NFP)
MAX_N_FIELDLINES = 200
flock = FlockOfFieldlines(MAX_N_FIELDLINES, IOTA, field, M_POL, N_TOR, NFP)

ds_datheta = 1.0 / PSI_TOR_EDGE  # as A_theta = PSI = PSI_TOR_EDGE*s_tor

lambda_a, lambda_b = flock.calc_offset_coefficients(R, ds_datheta)
nu_star_crit = flock.calc_nu_star_crit(R)

lam_a_ana = calc_lambda_a_analytic(
    ds_datheta, R, B_0, EPS_1, EPS_2, N_TOR, M_POL, flock.iota_p, IOTA
)
lam_b_ana = calc_lambda_b_analytic(ds_datheta, R, EPS_2, delta_B_1, flock.iota_p)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

print(f"Field: B = 1 - {EPS_1}*cos(2t-Np) - {EPS_2}*cos(t)*(1+cos(2t-Np))")
print(f"  iota={IOTA},  nfp={NFP},  R={R},  B_phi_cov={B_PHI_COV:.4f}")
print(f"Lambda_A     = {lambda_a:.6e}")
print(f"Lambda_A_ana = {lam_a_ana:.6e}  (rel.diff = {abs(lam_a_ana/lambda_a-1):.3f})")
print(f"Lambda_B     = {lambda_b:.6e}")
print(f"Lambda_B_ana = {lam_b_ana:.6e}")
print(f"nu_star_crit = {nu_star_crit:.6e}")

plot_offset_vs_nu_star(lambda_a, lambda_b, lam_a_ana, lam_b_ana)
if PLOT_FIELD:
    plot_field(field, NFP)
