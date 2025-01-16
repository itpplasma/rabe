import numpy as np


def write_stellarator_bc_file(
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
