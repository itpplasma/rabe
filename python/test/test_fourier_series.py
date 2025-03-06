# %%
from rabe.fourier_series import FourierSeries, evaluate
import numpy as np


def test_fourier_series_against_analytic():
    m = [0, 0, 1, 1, 1]
    n = [0, 1, -1, 0, 1]
    fmnc = [1.0, 0.0, 0.0, 0.0, 1.0]
    fmns = [0.0, 1.0, 0.0, 1.0, 0.0]
    angles = [[0, 0], [np.pi, 0], [0, np.pi / 2], [np.pi / 2, np.pi / 2]]
    f_check = [2.0, 0.0, 0.0, 2.0]

    fourier = FourierSeries(m, n, fmnc, fmns)
    f = []
    for i in range(len(angles)):
        f.append(evaluate(fourier, angles[i][0], angles[i][1]))

    np.testing.assert_allclose(f, f_check, atol=1e-12)


def test_fourier_series_run():
    m = np.array([0, 0, 1, 1, 1])
    n = np.array([0, 1, -1, 0, 1])
    nfp = 3
    n = n * nfp
    fmnc = np.array([1.0, 0.0, 0.0, 0.0, 1.0])
    fmns = np.array([0.0, 1.0, 0.0, 1.0, 0.0])
    fourier = FourierSeries(list(m), (nfp * n).tolist(), list(fmnc), fmns.tolist())
    _ = evaluate(fourier, np.pi / 4, -np.pi / 4)


if __name__ == "__main__":
    test_fourier_series_run()
    test_fourier_series_against_analytic()
