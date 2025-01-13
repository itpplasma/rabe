import numpy as np


def write_stellarator_symmetric_bc_file(
    filename: str,
    rmnc: np.array,
    zmns: np.array,
    vmns: np.array,
    bmnc: np.array,
    m: np.array,
    n: np.array,
    s_tor: np.array,
    iota: np.array,
    b_covar_pol: np.array,
    b_covar_tor: np.array,
    nfp: int,
    edge_toroidal_flux: float,
    minor_radius: float,
    major_radius: float,
):

    from libneo import write_boozer_head, append_boozer_block_head, append_boozer_block

    n_surf = len(s_tor)
    dummy = 0.0

    write_boozer_head(
        filename,
        "",
        0,
        np.max(m),
        np.max(n) // nfp,
        n_surf,
        nfp,
        edge_toroidal_flux,
        minor_radius,
        major_radius,
    )
    for i in range(n_surf):
        rmn = np.array([complex(x, 0) for x in rmnc[i]])
        zmn = np.array([complex(0, -x) for x in zmns[i]])
        vmn = np.array([complex(0, -x) for x in vmns[i]])
        bmn = np.array([complex(x, 0) for x in bmnc[i]])
        append_boozer_block_head(
            filename,
            s_tor[i],
            iota[i],
            b_covar_tor[i],
            b_covar_pol[i],
            dummy,
            dummy,
            nfp,
        )
        append_boozer_block(filename, m, n, rmn, zmn, vmn, bmn, nfp)


if __name__ == "__main__":
    import sys
    from simsopt.mhd.vmec import Vmec
    from simsopt.mhd.boozer import Boozer

    vmec_file = sys.argv[1]
    output_file = "booz_xform_field.bc"

    vmec = Vmec(vmec_file)
    boozer = Boozer(vmec)
    boozer.register(vmec.s_half_grid)
    boozer.run()

    write_stellarator_symmetric_bc_file(
        filename=output_file,
        rmnc=boozer.bx.rmnc_b.T,
        zmns=boozer.bx.zmns_b.T,
        vmns=-boozer.bx.numns_b.T * boozer.bx.nfp / (2 * np.pi),
        bmnc=boozer.bx.bmnc_b.T,
        m=boozer.bx.xm_b,
        n=-boozer.bx.xn_b,
        s_tor=boozer.bx.s_b,
        iota=boozer.bx.iota,
        b_covar_pol=boozer.bx.bsubumnc[:, 0],
        b_covar_tor=boozer.bx.bsubvmnc[:, 0],
        nfp=boozer.bx.nfp,
        edge_toroidal_flux=vmec.wout.phi[-1],
        minor_radius=vmec.wout.Aminor_p,
        major_radius=vmec.wout.Rmajor_p,
    )
