# %%
if __name__ == "__main__":
    from simsopt.mhd.vmec import Vmec
    from rabe.boozer_modes import read_modes_bc
    from rabe.coordinate_orientations import (
        get_bc_theta_orientation,
        get_vmec_theta_orientation,
    )

    vmec_file = "../external/data/wout_w7x.nc"
    vmec = Vmec(vmec_file)
    print("For W7x VMEC")
    get_vmec_theta_orientation(vmec)
    print("poloidal field")
    print(vmec.wout.bsubumnc[0, -1])
    print("toroidal field")
    print(vmec.wout.bsubvmnc[0, -1])

    par_bc_file = "../external/data/w7x_sc1.bc"
    get_mode_idx_par = lambda m, n, n_max: (
        (2 * n_max + 1) * (m > 1) * (m - 1) + (n_max + 1) * (m > 0) + n_max - n
    )
    rmnc, zmns, _, _ = read_modes_bc(par_bc_file, get_mode_idx=get_mode_idx_par)
    print("For W7x bc file")
    get_bc_theta_orientation(rmnc, zmns)
