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
    from libneo import (
        write_boozer_head,
        append_boozer_block_head,
        append_boozer_block_stellerator_symmetry,
    )

    n_surf = len(s_tor)
    dummy = 0.0

    write_boozer_head(
        filename,
        "",
        0,
        np.max(np.abs(m)),
        np.max(n) // nfp,
        n_surf,
        nfp,
        edge_toroidal_flux,
        minor_radius,
        major_radius,
    )
    for i in range(n_surf):
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
        append_boozer_block_stellerator_symmetry(
            filename, m, n, rmnc[i], zmns[i], vmns[i], bmnc[i], nfp
        )
