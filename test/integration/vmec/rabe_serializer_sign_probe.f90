!> Emit native boozmn-reader and RABE-NetCDF serializer values.
!>
!> The companion verifier independently decodes both generated NetCDF files,
!> evaluates the signed Fourier series, and checks attributes and array order.
program rabe_serializer_sign_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozmn_file, only: boozmn_data_t, write_boozmn
    use boozer_field, only: boozer_field_t
    use constants, only: pi
    use netcdf_mod, only: netcdf_t

    implicit none

    character(len=*), parameter :: boozmn_path = 'rabe_serializer_probe_boozmn.nc'
    character(len=*), parameter :: output_path = 'rabe_serializer_probe_output.nc'
    integer, parameter :: nmode = 4
    integer, parameter :: ns_full = 5
    integer, parameter :: nsurf = 3
    integer, parameter :: ncoeff = 5
    integer, parameter :: ixm(nmode) = [0, 1, 1, 2]
    integer, parameter :: ixn(nmode) = [0, 1, -1, 0]
    real(dp), parameter :: bmn(nmode) = [2.0_dp, 0.1_dp, -0.07_dp, 0.02_dp]
    real(dp), parameter :: signed_coefficients(ncoeff) = &
        [-1.25_dp, 2.5_dp, -3.75_dp, 5.0_dp, -6.25_dp]
    integer, parameter :: flags(ncoeff) = [1, 0, -1, 2, -2]
    real(dp), parameter :: stor = 0.5_dp
    real(dp), parameter :: theta = 2.0_dp*pi*7.0_dp/48.0_dp
    real(dp), parameter :: zeta = 2.0_dp*pi*11.0_dp/96.0_dp

    type(boozmn_data_t) :: data
    type(boozer_field_t) :: field
    type(netcdf_t) :: output
    real(dp) :: bmod, iota, btheta_cov, bzeta_cov

    data%ns = ns_full
    data%nfp = 1
    data%nmodes = nmode
    data%nsurf = nsurf
    data%asym = .false.
    allocate (data%jlist(nsurf), data%ixm(nmode), data%ixn(nmode))
    allocate (data%iota(ns_full), data%buco(ns_full), data%bvco(ns_full))
    allocate (data%phi(ns_full), data%bmnc(nmode, nsurf))
    allocate (data%rmnc(nmode, nsurf), data%zmns(nmode, nsurf))
    allocate (data%pmns(nmode, nsurf))

    data%jlist = [2, 3, 4]
    data%ixm = ixm
    data%ixn = ixn
    data%iota = [0.5_dp, 0.55_dp, 0.6_dp, 0.65_dp, 0.7_dp]
    data%buco = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp]
    data%bvco = [0.2_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.24_dp]
    data%phi = [0.0_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp]
    data%bmnc = spread(bmn, dim=2, ncopies=nsurf)
    data%rmnc = 0.0_dp
    data%rmnc(1, :) = [1.8_dp, 1.9_dp, 2.0_dp]
    data%zmns = 0.0_dp
    data%pmns = 0.0_dp

    call write_boozmn(boozmn_path, data)
    call field%init_from_boozmn(boozmn_path)
    if (.not. field%initialized) error stop 'boozmn field did not initialize'
    call field%fix_to_surface(stor)
    call field%compute_B_mod(theta, zeta, bmod)
    call field%get_iota(stor, iota)
    call field%get_covariant_components(btheta_cov, bzeta_cov)

    call output%create(output_path)
    call output%add_global_attribute('title', 'RABE serializer sign probe')
    call output%def_dim('coefficient', ncoeff)
    call output%add_real('psi_tor_edge')
    call output%add_attr('psi_tor_edge', 'units', 'Wb')
    call output%add_real('bmod')
    call output%add_attr('bmod', 'units', 'T')
    call output%add_real_1d('signed_coefficients', 'coefficient')
    call output%add_attr('signed_coefficients', 'long_name', &
        'signed bootstrap and Ware coefficient probe')
    call output%add_int_1d('flags', 'coefficient')
    call output%write_real('psi_tor_edge', field%psi_tor_edge)
    call output%write_real('bmod', bmod)
    call output%write_real_1d('signed_coefficients', signed_coefficients)
    call output%write_int_1d('flags', flags)
    call output%close()

    write (*, '(a,a)') 'RABE_SERIALIZER boozmn_path=', boozmn_path
    write (*, '(a,a)') 'RABE_SERIALIZER output_path=', output_path
    write (*, '(a,3(1x,es26.17))') 'RABE_SERIALIZER s_theta_zeta=', &
        stor, theta, zeta
    write (*, '(a,es26.17)') 'RABE_SERIALIZER bmod_T=', bmod
    write (*, '(a,es26.17)') 'RABE_SERIALIZER psi_tor_edge_Wb=', &
        field%psi_tor_edge
    write (*, '(a,es26.17)') 'RABE_SERIALIZER iota=', iota
    write (*, '(a,2(1x,es26.17))') 'RABE_SERIALIZER B_covariant_Tm=', &
        btheta_cov, bzeta_cov
    write (*, '(a,5(1x,es26.17))') 'RABE_SERIALIZER signed_coefficients=', &
        signed_coefficients
    write (*, '(a,5(1x,i0))') 'RABE_SERIALIZER flags=', flags
end program rabe_serializer_sign_probe
