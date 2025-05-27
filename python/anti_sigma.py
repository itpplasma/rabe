# %%
import matplotlib.pyplot as plt
import numpy as np

eps_0 = 0.5
eps_1 = 0.1
deltaB = 0.1


def anti_sigma(theta, phi):
    B = (
        1
        + eps_0 * np.cos(phi)
        + eps_1 * np.cos(theta) * (1 - np.cos(phi))
        + deltaB * np.cos(theta)
    )
    return B


phi = np.linspace(-np.pi, np.pi, 100)
theta = np.linspace(0, 2 * np.pi, 100)
Phi, Theta = np.meshgrid(phi, theta)

plt.figure()
plt.contour(Phi, Theta, anti_sigma(Theta, Phi), levels=20)
plt.grid(True)
plt.colorbar()
