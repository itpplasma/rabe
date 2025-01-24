# %%
if __name__ == "__main__":
    from simsopt.mhd.vmec import Vmec
    from rabe.boozer_modes import read_modes_bc
    from rabe.coordinate_orientations import (
        get_bc_theta_orientation,
        get_vmec_theta_orientation,
    )

    vmec_file = "output/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    vmec = Vmec(vmec_file)
    print("For Landreman VMEC")
    get_vmec_theta_orientation(vmec)

    par_bc_file = "output/quasi_symmetric.bc"
    get_mode_idx_par = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max + n) + (n_max + 1)
    )
    rmnc, zmns, _, _ = read_modes_bc(par_bc_file, get_mode_idx=get_mode_idx_par)
    print("For Landreman bc file")
    get_bc_theta_orientation(rmnc, zmns)
