# from simsopt.mhd.vmec import Vmec
# from simsopt.mhd.boozer import Boozer
# %%
def write_stellarator_symmetric_bc_file(
    filename: str,
    rmnc: list,
    zmns: list,
    vmns: list,
    bmnc: list,
    m: list,
    n: list,
    s_tor: list,
    iota: list,
    b_covar_pol: list,
    b_covar_tor: list,
    vp: list,
    nfp: int,
    edge_toroidal_flux: float,
    minor_radius: float,
    major_radius: float,
):

    from libneo import write_boozer_head, append_boozer_block_head, append_boozer_block

    n_surf = len(s_tor)

    write_boozer_head(
        filename,
        "",
        0,
        np.max(m),
        np.max(n),
        n_surf,
        nfp,
        edge_toroidal_flux,
        minor_radius,
        major_radius,
    )
    for i in range(n_surf):
        rmn = np.array([complex(x, 0) for x in rmnc[i]])
        zmn = np.array([complex(0, x) for x in zmns[i]])
        vmn = np.array([complex(0, x) for x in vmns[i]])
        bmn = np.array([complex(x, 0) for x in bmnc[i]])
        append_boozer_block_head(
            filename, s_tor[i], iota[i], b_covar_tor[i], b_covar_pol[i], 0, vp[i], nfp
        )
        append_boozer_block(filename, m, n, rmn, zmn, vmn, bmn, nfp)


if __name__ == "__main__":
    # import sys

    # vmec_file = sys.argv[1]
    # boozer = Boozer(Vmec(vmec_file))
    import numpy as np

    dummy_file = "dummy.bc"
    s_tor = np.linspace(0.1, 0.9, 3)
    m_max = 3
    dummy_int = 1
    dummy_float = 1.0
    dummy_int_list = [dummy_int] * m_max
    dummy_float_list = [dummy_float] * len(s_tor)
    dummy_list_list = [[dummy_float] * m_max] * len(s_tor)
    dummy_int_array = np.array(dummy_int_list)
    dummy_float_array = np.array(dummy_float)
    write_stellarator_symmetric_bc_file(
        filename=dummy_file,
        rmnc=dummy_list_list,
        zmns=dummy_list_list,
        vmns=dummy_list_list,
        bmnc=dummy_list_list,
        m=dummy_int_array,
        n=dummy_int_array,
        s_tor=s_tor,
        iota=dummy_float_list,
        b_covar_pol=dummy_float_list,
        b_covar_tor=dummy_float_list,
        vp=dummy_float_list,
        nfp=dummy_int,
        edge_toroidal_flux=dummy_float,
        minor_radius=dummy_float,
        major_radius=dummy_float,
    )
