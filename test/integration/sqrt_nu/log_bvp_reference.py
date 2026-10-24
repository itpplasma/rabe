"""Independent BVP and Kelvin-function oracles; Python standard library only.

No RABE implementation is imported. The log model solves gamma'' =
i*A**2*log(1/x)*gamma, gamma(0)=1, gamma'(1/2)=0 on a graded mesh.
The Kelvin oracle solves (z*gamma')'=i*A**2*gamma analytically as a
normalized I0(2*sqrt(i*A**2*z)), regular at z=0. Both test the WKB
approximation within its model regime, not agreement with experiment.
"""
import cmath
import math
import subprocess
import sys


def log_bvp(a, n):
    x = [0.5 * (j / n) ** 3 for j in range(n + 1)]
    diagonal = [0j] * (n + 1)
    upper = [0j] * n
    lower = [0j] * n
    rhs = [0j] * (n + 1)
    diagonal[0] = rhs[0] = 1.0
    for j in range(1, n):
        hm, hp = x[j] - x[j - 1], x[j + 1] - x[j]
        w = 2 / (hm + hp)
        lower[j - 1] = w / hm
        upper[j] = w / hp
        diagonal[j] = -w / hm - w / hp - 1j * a * a * math.log(1 / x[j])
    h = x[n] - x[n - 1]
    lower[n - 1] = 2 / h**2
    diagonal[n] = -2 / h**2 - 1j * a * a * math.log(1 / x[n])
    for j in range(1, n + 1):
        multiplier = lower[j - 1] / diagonal[j - 1]
        diagonal[j] -= multiplier * upper[j - 1]
        rhs[j] -= multiplier * rhs[j - 1]
    gamma = [0j] * (n + 1)
    gamma[n] = rhs[n] / diagonal[n]
    for j in range(n - 1, -1, -1):
        gamma[j] = (rhs[j] - upper[j] * gamma[j + 1]) / diagonal[j]
    return -sum((x[j + 1] - x[j]) * (gamma[j].real + gamma[j + 1].real) / 2
                for j in range(n))


def kelvin_integral(a):
    # I0 series and its exact integral over z in [0,1]. The additional
    # 1/(k+1) comes from integrating z**k, so the oracle uses no WKB kernel.
    term = 1 + 0j
    value = integral = term
    for k in range(1, 400):
        term *= 1j * a * a / (k * k)
        value += term
        integral += term / (k + 1)
        if abs(term) < 1e-17 * abs(value):
            return (integral / value).real
    raise AssertionError('Kelvin series did not converge')


seen = set()
output = subprocess.check_output([sys.argv[1], '--reference-output'], text=True)
for line in output.splitlines():
    parts = line.split()
    if not parts or parts[0] not in ('LOG', 'KELVIN'):
        continue
    kind, a, wkb = parts[0], float(parts[1]), float(parts[2])
    seen.add((kind, a))
    if kind == 'LOG':
        low, reference = log_bvp(a, 8000), log_bvp(a, 16000)
        assert abs(low - reference) < 2e-6 * abs(reference), (a, low, reference)
        error = abs(wkb - reference) / abs(reference)
        assert error < 0.065, ('log-model WKB/BVP error', a, error)
    else:
        reference = kelvin_integral(a)
        error = abs(wkb - reference) / abs(reference)
        assert error < 0.04, ('Kelvin WKB/exact error', a, error)
    print(f'{kind} A={a:g}: WKB={wkb:.10e}, independent={reference:.10e}, '
          f'relative error={error:.4%}')

assert seen == {('LOG', 10.0), ('LOG', 100.0), ('LOG', 1000.0),
                ('KELVIN', 10.0), ('KELVIN', 30.0)}, seen
