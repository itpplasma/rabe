# %%
def get_vmec_theta_orientation(vmec):
    idx_m1n0 = vmec.wout.ntor * 2 + 1
    dr_dtheta = -vmec.wout.rmnc[idx_m1n0, -1]  # approximately
    dz_dtheta = vmec.wout.zmns[idx_m1n0, -1]  # approximately

    if dr_dtheta * dz_dtheta < 0:
        print("Theta goes counter-clockwise")
    elif dr_dtheta * dz_dtheta > 0:
        print("Theta goes clockwise")
    else:
        print("Can not determine orientation")


def get_bc_theta_orientation(rmnc, zmns):
    idx_m1n0 = rmnc.get_mode_idx(m=1, n=0)
    dr_dtheta = -rmnc.coefs[-1, idx_m1n0]  # approximately
    dz_dtheta = zmns.coefs[-1, idx_m1n0]  # approximately

    if dr_dtheta * dz_dtheta < 0:
        print("Theta goes counter-clockwise")
    elif dr_dtheta * dz_dtheta > 0:
        print("Theta goes clockwise")
    else:
        print("Can not determine orientation")
