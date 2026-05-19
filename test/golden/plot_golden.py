"""Plot rabe.nc output from the golden record run."""

import sys
import xarray as xr
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "rabe.nc"
ds = xr.open_dataset(path)
s = ds["s_tor"].values

has_sc = "lambda_SC_bB" in ds

if has_sc:
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    ax1, ax2, ax3, ax4 = axes.flat
else:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(s, ds["Lambda_A"].values, "o-", label=r"$\Lambda_A$ (1/$\sqrt{\nu_*}$)")
ax1.plot(s, ds["Lambda_B"].values, "s-", label=r"$\Lambda_B$ (1/$\nu_*$)")
ax1.set_xlabel(r"$s_\mathrm{tor}$")
ax1.set_ylabel("coefficient [1]")
ax1.set_title("off-set")
ax1.legend()
ax1.grid(True)

ax2.plot(s, ds["Lambda_S"].values, "^-")
ax2.set_xlabel(r"$s_\mathrm{tor}$")
ax2.set_ylabel("coefficient [1]")
ax2.set_title(r"$\Lambda_\mathrm{S}$")
ax2.grid(True)

if has_sc:
    ax3.plot(s, ds["lambda_SC_bB"].values, "o-")
    ax3.set_xlabel(r"$s_\mathrm{tor}$")
    ax3.set_ylabel("coefficient [1]")
    ax3.set_title(r"$\lambda^\mathrm{SC}_{bB}$ (omnigenous Shaing-Callen)")
    ax3.grid(True)

    ax4.plot(s, ds["remainder"].values, "s-")
    ax4.set_xlabel(r"$s_\mathrm{tor}$")
    ax4.set_ylabel("coefficient [1]")
    ax4.set_title("remainder (non-omnigenous Shaing-Callen)")
    ax4.grid(True)

fig.suptitle(f"rabe output — {path}")
fig.tight_layout()
plt.savefig("rabe_output.png", dpi=150)
print("saved rabe_output.png")
plt.show()
