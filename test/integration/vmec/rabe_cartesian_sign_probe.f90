!> Emit native RABE Boozer-field quantities for an independent Cartesian sign gate.
!>
!> This is diagnostic-only source.  The verifier in plasma-sign-conventions reads
!> the same chartmap geometry, differentiates its Cartesian map independently,
!> and reconstructs B = B^theta e_theta + B^zeta e_zeta.  Grid nodes are used
!> deliberately so the native field spline and the stored chartmap payload can
!> be compared without an independent interpolation convention entering the gate.
program rabe_cartesian_sign_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_field, only: boozer_field_t
    use fieldline_integrands, only: local_radial_drift
    use fieldline_mod, only: fieldline_t
    use fourier_field, only: fourier_field_init, fourier_field_t

    implicit none

    character(len=*), parameter :: chartmap_file = 'input/circ_chartmap.nc'
    integer, parameter :: n_points = 4
    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
    integer, parameter :: rho_index(n_points) = [5, 9, 13, 17]
    real(dp), parameter :: rho(n_points) = &
        0.001_dp + real(rho_index, dp)*(1.0_dp - 0.001_dp)/19.0_dp
    integer, parameter :: theta_index(n_points) = [3, 8, 14, 19]
    integer, parameter :: zeta_index(n_points) = [4, 10, 16, 21]

    type(boozer_field_t) :: field
    type(fourier_field_t) :: fourier
    type(fieldline_t) :: fieldline
    real(dp) :: u(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: iota, btheta_cov, bphi_cov
    real(dp) :: theta_test, zeta_test, eta_test, theta_along, drift
    real(dp) :: bmod_shifted, bder_shifted(3)
    integer, parameter :: m_modes(4) = [0, 2, -1, 3]
    integer, parameter :: n_modes(4) = [0, 1, -2, -1]
    real(dp), parameter :: b_modes(4) = [2.5_dp, 0.17_dp, -0.09_dp, 0.04_dp]
    integer :: point

    call field%init_from_chartmap(chartmap_file)
    if (.not. field%initialized) error stop 'RABE chartmap field did not initialize'

    write (*, '(a,es26.17)') 'RABE_PROBE nfp=', field%nfp
    write (*, '(a,es26.17)') 'RABE_PROBE psi_tor_edge_Wb=', field%psi_tor_edge
    write (*, '(a,es26.17)') 'RABE_PROBE major_radius_m=', field%R

    do point = 1, n_points
        u = [rho(point)**2, &
            twopi*real(theta_index(point), dp)/24.0_dp, &
            twopi*real(zeta_index(point), dp)/24.0_dp]
        call field%evaluate(u, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        call field%get_iota(u(1), iota)
        call field%fix_to_surface(u(1))
        call field%get_covariant_components(btheta_cov, bphi_cov)

        write (*, '(a,i0)') 'RABE_PROBE point=', point
        write (*, '(a,3(1x,es26.17))') 'RABE_PROBE u_s_thetaB_zetaB=', u
        write (*, '(a,es26.17)') 'RABE_PROBE bmod_T=', bmod
        write (*, '(a,es26.17)') 'RABE_PROBE sqrtg_m3=', sqrtg
        write (*, '(a,3(1x,es26.17))') 'RABE_PROBE grad_log_bmod_covariant=', bder
        write (*, '(a,3(1x,es26.17))') 'RABE_PROBE bhat_covariant_m=', hcovar
        write (*, '(a,3(1x,es26.17))') 'RABE_PROBE bhat_contravariant_per_m=', hctrvr
        write (*, '(a,3(1x,es26.17))') 'RABE_PROBE curl_bhat_contravariant_per_m2=', &
            hcurl
        write (*, '(a,es26.17)') 'RABE_PROBE iota=', iota
        write (*, '(a,2(1x,es26.17))') 'RABE_PROBE B_covariant_Tm=', &
            btheta_cov, bphi_cov
    end do

    ! Phase-sensitive manufactured spectrum.  Positive and negative n labels,
    ! NFP packing, derivatives, the field-line slope, and the local radial-drift
    ! sign are all emitted for independent analytic reconstruction.
    call fourier_field_init(fourier, m_modes, n_modes, b_modes, &
        -0.13_dp, 2.7_dp, nfp=3, n_grid=257)
    theta_test = 0.73_dp
    zeta_test = 0.41_dp
    eta_test = 0.21_dp
    call fourier%compute_B_and_dB_dx(theta_test, zeta_test, bmod, bder)
    drift = local_radial_drift(fourier, theta_test, zeta_test, eta_test)

    fieldline%theta_0 = -0.37_dp
    fieldline%phi_0 = 0.22_dp
    fieldline%iota = -0.61_dp
    theta_along = fieldline%get_theta(zeta_test)

    write (*, '(a,3(1x,i0))') 'RABE_PROBE fourier_nfp_mn_count=', &
        nint(fourier%nfp), fourier%mn_max, size(m_modes)
    write (*, '(a,4(1x,i0))') 'RABE_PROBE fourier_m=', m_modes
    write (*, '(a,4(1x,i0))') 'RABE_PROBE fourier_n_normalized=', n_modes
    write (*, '(a,4(1x,es26.17))') 'RABE_PROBE fourier_Bmn_T=', b_modes
    write (*, '(a,2(1x,es26.17))') 'RABE_PROBE fourier_thetaB_zetaB=', &
        theta_test, zeta_test
    write (*, '(a,es26.17)') 'RABE_PROBE fourier_bmod_T=', bmod
    write (*, '(a,3(1x,es26.17))') 'RABE_PROBE fourier_dB_du_T=', bder
    call fourier%compute_B_and_dB_dx(theta_test, zeta_test + twopi/fourier%nfp, &
        bmod_shifted, bder_shifted)
    write (*, '(a,es26.17)') 'RABE_PROBE fourier_period_shift_bmod_T=', bmod_shifted
    write (*, '(a,3(1x,es26.17))') 'RABE_PROBE fourier_period_shift_dB_du_T=', &
        bder_shifted
    write (*, '(a,es26.17)') 'RABE_PROBE fourier_eta_per_T=', eta_test
    write (*, '(a,es26.17)') 'RABE_PROBE fourier_local_radial_drift=', drift
    write (*, '(a,4(1x,es26.17))') 'RABE_PROBE fieldline_phi_theta0_phi0_iota=', &
        zeta_test, fieldline%theta_0, fieldline%phi_0, fieldline%iota
    write (*, '(a,es26.17)') 'RABE_PROBE fieldline_theta=', theta_along

end program rabe_cartesian_sign_probe
