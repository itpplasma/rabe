import numpy as np
from rabe.write_bc import write_stellarator_bc_file

if __name__ == "__main__":
    import sys
    from simsopt.mhd.vmec import Vmec
    from simsopt.mhd.boozer import Boozer

    vmec_file = sys.argv[1]
    output_file = sys.argv[2]

    vmec = Vmec(vmec_file)
    boozer = Boozer(vmec)
    boozer.register(vmec.s_half_grid)
    boozer.run()

    write_stellarator_bc_file(
        filename=output_file,
        rmnc=boozer.bx.rmnc_b.T,
        zmns=boozer.bx.zmns_b.T,
        vmns=-boozer.bx.numns_b.T * boozer.bx.nfp / (2 * np.pi),
        bmnc=boozer.bx.bmnc_b.T,
        m=boozer.bx.xm_b,
        n=boozer.bx.xn_b,
        s_tor=boozer.bx.s_b,
        iota=boozer.bx.iota,
        b_covar_pol=boozer.bx.bsubumnc[0, :],
        b_covar_tor=boozer.bx.bsubvmnc[0, :],
        nfp=boozer.bx.nfp,
        edge_toroidal_flux=vmec.wout.phi[-1],
        minor_radius=vmec.wout.Aminor_p,
        major_radius=vmec.wout.Rmajor_p,
    )
