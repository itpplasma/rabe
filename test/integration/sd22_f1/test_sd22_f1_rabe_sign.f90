! F1 (SD2.2 criterion v3, section 2): sign of rabe calc_lambda_LC on a synthetic
! tokamak-like field (M_pol=1, N_tor=0) against the sign of
! lambda_bB^tok = 1.46*sigma_sqrtg/(iota*sqrt(eps_t)) (SQuID_Goodman_2024 README),
! for sigma_sqrtg = +-1 and both signs of iota. s_R = sign(LC)/sign(lambda_tok).
program test_sd22_f1_rabe_sign
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use constants, only: pi
    use mock_field, only: mock_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use coefficients, only: calc_gradient_scaling_factor_r_eff
    use shaing_callen_mod, only: calc_lambda_LC
    implicit none
    integer, parameter :: n_fieldlines = 24, n_eta = 100
    real(dp), parameter :: eps = 0.01_dp, psi_edge = 1.0_dp
    type(mock_field_t) :: field
    type(flock_of_fieldlines_t) :: flock
    real(dp) :: xi_0(n_fieldlines), iota, dr, lc, tok, ratio, s_r
    integer :: isigma, iiota, sigma, k
    logical :: failed

    failed = .false.
    xi_0 = [(2.0_dp*pi*real(k - 1, dp)/real(n_fieldlines, dp), k=1, n_fieldlines)]
    ! mock B = B_0*(1 + eps_0*cos(chi)) with eps_t = |eps_0|; covariant B_phi = +1.
    call field%mock_field_init(1.0_dp, 0.0_dp, 1.0_dp, -eps)
    print '(a)', 'F1 rabe: sigma_sqrtg  iota_sign  LC  lambda_tok  LC/(dr*1.46*sqrt(eps)/iota)  s_R'
    do isigma = 1, 2
        sigma = merge(1, -1, isigma == 1)
        do iiota = 1, 2
            iota = merge(1.05_dp, -1.05_dp, iiota == 1)
           call make_flock_from_labels(flock, xi_0, iota, field, 1.0_dp, 0.0_dp, 1.0_dp)
            dr = calc_gradient_scaling_factor_r_eff(flock, psi_edge, sigma)
            lc = calc_lambda_LC(flock, field, n_eta, dr)
            tok = 1.46_dp*real(sigma, dp)/(iota*sqrt(eps))
            ! LC = dr_dAtheta*(B_phi/iota)*f_t with f_t = 1.46*sqrt(eps_t); the
            ! README normalises by eps_t (D31*), so only the sign is compared with
            ! tok, while the magnitude is compared with dr*1.46*sqrt(eps)/iota.
            ratio = lc/(dr*1.46_dp*sqrt(eps)/iota)
            s_r = sign(1.0_dp, lc)*sign(1.0_dp, tok)
            print '(2i4,4es14.5,f6.1)', sigma, iiota, lc, tok, ratio, s_r
            if (s_r /= 1.0_dp) then
                print *, 'F1 FAIL: s_R /= +1 for sigma, iota ', sigma, iota
                failed = .true.
            end if
            ! Magnitude check of the trapped fraction against 1.46*sqrt(eps).
            if (abs(ratio - 1.0_dp) > 0.05_dp) then
                print *, 'F1 FAIL: |LC| differs from the tokamak estimate by >5%'
                failed = .true.
            end if
        end do
    end do
    if (failed) error stop 'F1 rabe sign test failed'
    print *, 'F1 rabe: s_R = +1 for sigma_sqrtg = +-1 and both iota signs'
end program test_sd22_f1_rabe_sign
