"""``rabe.ad``: NetCDF-free, flang-compiled AD backend (coexistence spike).

This subpackage deliberately does NOT mirror the f90wrap OO API used by the
gfortran ``rabe`` package (``rabe.FourierField``, ``rabe.FlockOfFieldlines``,
...). Instead it exposes a single flat function over plain numpy arrays,
calling straight into a ``bind(C)`` Fortran symbol via :mod:`ctypes` - no
Fortran derived type or object handle ever crosses from this extension into
the gfortran one, so there is nothing that looks-alike-but-is-incompatible
between the two.

This is also the shape autodiff wants: Enzyme differentiates a plain
``inputs[] -> outputs[]`` function. The eventual adjoint will be exposed as
one more ``bind(C)`` symbol in the same shared library, wired up as one more
function here - no change to the calling convention.

If the flang backend was not built (``-DBUILD_AD_BACKEND=ON`` was not passed
at configure time, or no flang toolchain was available), importing this
module still succeeds, but calling into it raises a clear
:class:`NotImplementedError` rather than crashing on a missing library.
"""

from __future__ import annotations

import ctypes
import glob
import os

import numpy as np

_LIB = None
_LIB_ERROR: Exception | None = None


def _find_library_path() -> str | None:
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = sorted(
        glob.glob(os.path.join(here, "_rabe_ad*.so"))
        + glob.glob(os.path.join(here, "_rabe_ad*.dylib"))
        + glob.glob(os.path.join(here, "_rabe_ad*.dll"))
    )
    return candidates[0] if candidates else None


def _load_library():
    global _LIB, _LIB_ERROR
    if _LIB is not None or _LIB_ERROR is not None:
        return
    path = _find_library_path()
    if path is None:
        _LIB_ERROR = NotImplementedError(
            "rabe.ad backend not built: no _rabe_ad shared library found "
            "next to this module. Configure with -DBUILD_AD_BACKEND=ON and "
            "a flang toolchain to enable rabe.ad (see cmake/ad_backend.cmake)."
        )
        return
    try:
        lib = ctypes.CDLL(path)
        fn = lib.rabe_fourier_offset_coefficients
        fn.restype = None
        fn.argtypes = [
            ctypes.c_int,  # mn_max
            np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # m
            np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # n
            np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # B_mn
            ctypes.c_double,  # B_theta_cov
            ctypes.c_double,  # B_phi_cov
            ctypes.c_int,  # nfp
            ctypes.c_int,  # n_grid
            ctypes.c_double,  # iota
            ctypes.c_double,  # M_pol
            ctypes.c_double,  # N_tor
            ctypes.c_int,  # max_n_fieldlines
            ctypes.c_double,  # R
            ctypes.c_double,  # dr_dAtheta
            ctypes.POINTER(ctypes.c_double),  # lambda_a (out)
            ctypes.POINTER(ctypes.c_double),  # lambda_b (out)
            ctypes.POINTER(ctypes.c_double),  # nu_star_crit (out)
        ]
        _LIB = fn
    except OSError as exc:  # pragma: no cover - environment dependent
        _LIB_ERROR = NotImplementedError(
            f"rabe.ad backend found at {path!r} but failed to load: {exc}"
        )


def offset_coefficients(
    m,
    n,
    B_mn,
    B_theta_covariant,
    B_phi_covariant,
    nfp,
    n_grid,
    iota,
    M_pol,
    N_tor,
    max_n_fieldlines,
    R,
    dr_dAtheta,
):
    """Compute (Lambda_A, Lambda_B, nu_star_crit) for a Fourier-mode field.

    Mirrors the gfortran path ``FlockOfFieldlines.calc_offset_coefficients`` /
    ``calc_nu_star_crit`` used in ``python/example_fourier.py``, but runs
    entirely through the NetCDF-free flang-compiled ``_rabe_ad`` extension.
    No adjoints yet - this only proves the flang/gfortran coexistence and
    numerical path.
    """
    _load_library()
    if _LIB is None:
        raise _LIB_ERROR

    m = np.ascontiguousarray(m, dtype=np.int32)
    n = np.ascontiguousarray(n, dtype=np.int32)
    B_mn = np.ascontiguousarray(B_mn, dtype=np.float64)
    if not (m.shape == n.shape == B_mn.shape):
        raise ValueError("m, n, B_mn must have the same shape")

    lambda_a = ctypes.c_double()
    lambda_b = ctypes.c_double()
    nu_star_crit = ctypes.c_double()

    _LIB(
        m.size,
        m,
        n,
        B_mn,
        float(B_theta_covariant),
        float(B_phi_covariant),
        int(nfp),
        int(n_grid),
        float(iota),
        float(M_pol),
        float(N_tor),
        int(max_n_fieldlines),
        float(R),
        float(dr_dAtheta),
        ctypes.byref(lambda_a),
        ctypes.byref(lambda_b),
        ctypes.byref(nu_star_crit),
    )
    return lambda_a.value, lambda_b.value, nu_star_crit.value
