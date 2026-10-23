program test_wout_chain_modes
    !! Input mode numbers of a wout chain on libneo's default construction:
    !! the (m, n) of the wout harmonics that a chain tangent indexes must match
    !! the file's xm and xn, read here directly with netCDF. The check needs
    !! no rabe field or spline state and uses the Landreman-Paul QA wout.
    use constants, only: dp
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, &
                      nf90_inq_dimid, nf90_inquire_dimension, NF90_NOWRITE, &
                      NF90_NOERR
    use wout_chain, only: wout_chain_t, wout_chain_init, wout_chain_modes

    implicit none

    integer, parameter :: mboz = 6, nboz = 6, grid_refinement = 6
    type(wout_chain_t) :: chain
    character(len=*), parameter :: wout_file = &
        "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
    integer, allocatable :: xm(:), xn(:), xm_file(:), xn_file(:)

    call wout_chain_init(chain, wout_file, mboz, nboz, grid_refinement)
    call read_file_modes(wout_file, xm_file, xn_file)

    call wout_chain_modes(chain, xm, xn)

    if (size(xm) /= chain%ctx%nstrm .or. size(xn) /= chain%ctx%nstrm) &
        error stop "wout chain mode count differs from the tangent shape"
    if (size(xm) /= size(xm_file)) error stop "wout chain mode count differs from the file"
    if (any(xm /= xm_file) .or. any(xn /= xn_file)) &
        error stop "wout chain modes differ from the file's xm, xn"
    print "(A,I0,A)", "wout chain modes match the file (", size(xm), " harmonics)"

contains

    subroutine read_file_modes(path, xm_file, xn_file)
        character(len=*), intent(in) :: path
        integer, allocatable, intent(out) :: xm_file(:), xn_file(:)
        real(dp), allocatable :: buf(:)
        integer :: ncid, dimid, varid, mnmax

        call check(nf90_open(path, NF90_NOWRITE, ncid))
        call check(nf90_inq_dimid(ncid, "mn_mode", dimid))
        call check(nf90_inquire_dimension(ncid, dimid, len=mnmax))
        allocate (buf(mnmax))
        call check(nf90_inq_varid(ncid, "xm", varid))
        call check(nf90_get_var(ncid, varid, buf))
        xm_file = nint(buf)
        call check(nf90_inq_varid(ncid, "xn", varid))
        call check(nf90_get_var(ncid, varid, buf))
        xn_file = nint(buf)
        call check(nf90_close(ncid))
    end subroutine read_file_modes

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) error stop "netCDF error reading the wout file"
    end subroutine check

end program test_wout_chain_modes
