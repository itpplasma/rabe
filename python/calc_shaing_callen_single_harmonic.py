# %%
from scipy.integrate import quad
import numpy as np
from rabe.shaing_callen import find_maximum

if __name__ == "__main__":
    b_harmonic = lambda alpha: 1.0 - 0.0001 * np.cos(alpha)
    _, b_0 = find_maximum(b_harmonic, 0.0, 2.0 * np.pi)

    def integrand(x):
        def f(alpha):
            b = b_harmonic(alpha)
            return np.sqrt(1 - x * b / b_0) / (b / b_0) ** 2.0

        alpha_average, _ = quad(f, 0, 2 * np.pi)
        alpha_average = alpha_average / (2 * np.pi)
        return x / alpha_average

    integral, _ = quad(integrand, 0, 1)
    print(1.0 - 3.0 / 4.0 * integral)
