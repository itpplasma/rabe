program test_booz_xform_signs
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_field, only: boozxform_field_t
    use netcdf
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: boozmn_file = 'test_booz_xform_signs.nc'
    integer, parameter :: nmode = 4
    integer, dimension(nmode), parameter :: ixm = [0, 1, 1, 2]
    integer, dimension(nmode), parameter :: ixn = [0, 1, -1, 0]
    real(dp), dimension(nmode), parameter :: bmn = [2.0_dp, 0.1_dp, -0.07_dp, 0.02_dp]
    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: stor = 0.5_dp
    real(dp), parameter :: theta = twopi*7.0_dp/48.0_dp
    real(dp), parameter :: zeta = twopi*11.0_dp/96.0_dp
    real(dp), parameter :: reltol = 1.0e-6_dp
    real(dp), parameter :: abstol = 1.0e-10_dp

    type(boozxform_field_t) :: field
    real(dp) :: bmod, expected

    call write_boozmn_fixture(boozmn_file)
    call field%boozer_field_init(boozmn_file)
    if (.not. field%initialized) error stop 'booz_xform field not initialized'

    call field%fix_to_surface(stor)
    call field%compute_B_mod(theta, zeta, bmod)
    expected = expected_bmod(theta, zeta)

    print *, 'booz_xform signed-mode Bmod:', bmod, ' expected:', expected
    if (not_same(bmod, expected, reltol_in=reltol, abstol_in=abstol)) then
        error stop 'booz_xform signed-mode Bmod mismatch'
    end if

contains

    function expected_bmod(theta_in, zeta_in) result(value)
        real(dp), intent(in) :: theta_in, zeta_in
        real(dp) :: value
        integer :: imode

        value = 0.0_dp
        do imode = 1, nmode
            value = value + bmn(imode)*cos(real(ixm(imode), dp)*theta_in &
                                           - real(ixn(imode), dp)*zeta_in)
        end do
    end function expected_bmod

    subroutine write_boozmn_fixture(path)
        character(len=*), intent(in) :: path
        integer, parameter :: ns_full = 5
        integer, parameter :: nsurf = 3
        integer :: ncid, dim_radius, dim_mode, dim_surf
        integer :: var_ns, var_nfp, var_mnboz, var_lasym
        integer :: var_jlist, var_ixm, var_ixn
        integer :: var_iota, var_buco, var_bvco, var_phi
        integer :: var_bmnc, var_rmnc, var_zmns, var_pmns
        integer, dimension(nsurf) :: jlist
        real(dp), dimension(ns_full) :: iota, buco, bvco, phi
        real(dp), dimension(nmode, nsurf) :: coeff

        call check_nc(nf90_create(path, nf90_clobber, ncid), 'create boozmn')
        call check_nc(nf90_def_dim(ncid, 'radius', ns_full, dim_radius), &
                      'define radius')
        call check_nc(nf90_def_dim(ncid, 'mn_mode', nmode, dim_mode), &
                      'define mn_mode')
        call check_nc(nf90_def_dim(ncid, 'comput_surfs', nsurf, dim_surf), &
                      'define comput_surfs')

        call check_nc(nf90_def_var(ncid, 'ns_b', nf90_int, varid=var_ns), &
                      'define ns_b')
        call check_nc(nf90_def_var(ncid, 'nfp_b', nf90_int, varid=var_nfp), &
                      'define nfp_b')
        call check_nc(nf90_def_var(ncid, 'mnboz_b', nf90_int, varid=var_mnboz), &
                      'define mnboz_b')
        call check_nc(nf90_def_var(ncid, 'lasym__logical__', nf90_int, &
                      varid=var_lasym), 'define lasym')
        call check_nc(nf90_def_var(ncid, 'jlist', nf90_int, [dim_surf], var_jlist), &
                      'define jlist')
        call check_nc(nf90_def_var(ncid, 'ixm_b', nf90_int, [dim_mode], var_ixm), &
                      'define ixm_b')
        call check_nc(nf90_def_var(ncid, 'ixn_b', nf90_int, [dim_mode], var_ixn), &
                      'define ixn_b')
        call check_nc(nf90_def_var(ncid, 'iota_b', nf90_double, [dim_radius], &
                      var_iota), 'define iota_b')
        call check_nc(nf90_def_var(ncid, 'buco_b', nf90_double, [dim_radius], &
                      var_buco), 'define buco_b')
        call check_nc(nf90_def_var(ncid, 'bvco_b', nf90_double, [dim_radius], &
                      var_bvco), 'define bvco_b')
        call check_nc(nf90_def_var(ncid, 'phi_b', nf90_double, [dim_radius], &
                      var_phi), 'define phi_b')
        call check_nc(nf90_def_var(ncid, 'bmnc_b', nf90_double, &
                      [dim_mode, dim_surf], var_bmnc), 'define bmnc_b')
        call check_nc(nf90_def_var(ncid, 'rmnc_b', nf90_double, &
                      [dim_mode, dim_surf], var_rmnc), 'define rmnc_b')
        call check_nc(nf90_def_var(ncid, 'zmns_b', nf90_double, &
                      [dim_mode, dim_surf], var_zmns), 'define zmns_b')
        call check_nc(nf90_def_var(ncid, 'pmns_b', nf90_double, &
                      [dim_mode, dim_surf], var_pmns), 'define pmns_b')
        call check_nc(nf90_enddef(ncid), 'end definitions')

        jlist = [2, 3, 4]
        iota = [0.5_dp, 0.55_dp, 0.6_dp, 0.65_dp, 0.7_dp]
        buco = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp]
        bvco = [0.2_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.24_dp]
        phi = [0.0_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp]

        call check_nc(nf90_put_var(ncid, var_ns, ns_full), 'write ns_b')
        call check_nc(nf90_put_var(ncid, var_nfp, 1), 'write nfp_b')
        call check_nc(nf90_put_var(ncid, var_mnboz, nmode), 'write mnboz_b')
        call check_nc(nf90_put_var(ncid, var_lasym, 0), 'write lasym')
        call check_nc(nf90_put_var(ncid, var_jlist, jlist), 'write jlist')
        call check_nc(nf90_put_var(ncid, var_ixm, ixm), 'write ixm_b')
        call check_nc(nf90_put_var(ncid, var_ixn, ixn), 'write ixn_b')
        call check_nc(nf90_put_var(ncid, var_iota, iota), 'write iota_b')
        call check_nc(nf90_put_var(ncid, var_buco, buco), 'write buco_b')
        call check_nc(nf90_put_var(ncid, var_bvco, bvco), 'write bvco_b')
        call check_nc(nf90_put_var(ncid, var_phi, phi), 'write phi_b')

        coeff = spread(bmn, dim=2, ncopies=nsurf)
        call check_nc(nf90_put_var(ncid, var_bmnc, coeff), 'write bmnc_b')
        coeff = 0.0_dp
        coeff(1, :) = [1.8_dp, 1.9_dp, 2.0_dp]
        call check_nc(nf90_put_var(ncid, var_rmnc, coeff), 'write rmnc_b')
        coeff = 0.0_dp
        call check_nc(nf90_put_var(ncid, var_zmns, coeff), 'write zmns_b')
        call check_nc(nf90_put_var(ncid, var_pmns, coeff), 'write pmns_b')
        call check_nc(nf90_close(ncid), 'close boozmn')
    end subroutine write_boozmn_fixture

    subroutine check_nc(status, context)
        integer, intent(in) :: status
        character(len=*), intent(in) :: context

        if (status /= nf90_noerr) then
            print *, trim(context), ': ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check_nc

end program test_booz_xform_signs
