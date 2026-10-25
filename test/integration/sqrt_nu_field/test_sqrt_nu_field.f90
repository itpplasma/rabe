program test_sqrt_nu_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines
    use coefficients, only: calc_offset_coefficients, calc_gradient_scaling_factor_r_eff
    use sqrt_nu_geometry_mod, only: sqrt_nu_geometry_t, compute_sqrt_nu_geometry
    use sqrt_nu_mod, only: sqrt_nu_result_t, evaluate_sqrt_nu_modes
    implicit none
    type(anti_sigma_field_t), target :: base
    type(mock_perturbed_field_t) :: field
    type(flock_of_fieldlines_t) :: flock
    type(sqrt_nu_geometry_t) :: coarse, fine, finest, more_lines, most_lines
    real(dp), parameter :: nu_star(2) = [0.001_dp, 0.002_dp]
    real(dp) :: lambda_A, lambda_B, dr, expected, anchor
    real(dp) :: value_33, value_65, value_129, anchor_16, anchor_32, anchor_64
    integer :: nmode, inu

    call check_reference_origin()
    call base%anti_sigma_field_init(0.0_dp, 1.0_dp, 1.0_dp, -0.2_dp, -0.002_dp)
    call field%mock_perturbed_field_init(base, 1.0_dp, 0.0_dp, 0.001_dp)
    call make_flock_of_fieldlines(flock, 16, 0.5_dp, field, 0.0_dp, 1.0_dp, 1.0_dp)
    nmode = (size(flock%fieldlines) - 1)/2
    call compute_sqrt_nu_geometry(field, flock, 3.0_dp, 1.0_dp, 1, 1.0_dp, &
                                  nu_star, 33, nmode, 0.001_dp, coarse)
    call compute_sqrt_nu_geometry(field, flock, 3.0_dp, 1.0_dp, 1, 1.0_dp, &
                                  nu_star, 65, nmode, 0.001_dp, fine)
    call compute_sqrt_nu_geometry(field, flock, 3.0_dp, 1.0_dp, 1, 1.0_dp, &
                                  nu_star, 129, nmode, 0.001_dp, finest)
    dr = calc_gradient_scaling_factor_r_eff(flock, 1.0_dp, 1)
    call calc_offset_coefficients(flock, 3.0_dp, dr, lambda_A, lambda_B)
    do inu = 1, size(nu_star)
        expected = lambda_A/sqrt(nu_star(inu)) + lambda_B/nu_star(inu)
        anchor = fine%normalization*sum(fine%g_modes(:, inu)*fine%iv_eta_b)
        if (abs(anchor - expected) > 1.e-11_dp*max(abs(expected), 1.e-12_dp)) &
            error stop 'field harmonic normalization does not reproduce native E77'
    end do
    if (.not. all(ieee_is_finite(fine%tau_b))) error stop 'nonfinite bounce table'
    if (any(fine%tau_b <= 0.0_dp)) error stop 'nonpositive bounce time'
    if (maxval(abs(fine%tau_d_per_nu - fine%tau_d_per_nu(1))) > 1.e-14_dp) &
        error stop 'FEdef boundary diffusion was not frozen'
    value_33 = evaluate_field(coarse)
    value_65 = evaluate_field(fine)
    value_129 = evaluate_field(finest)
    call check_refinement(abs(value_65 - value_33), abs(value_129 - value_65), &
                          abs(value_129), 'pitch points33/65/129')

    anchor_16 = fine%normalization*sum(fine%g_modes(:2, 1)*fine%iv_eta_b(:2))
    call make_flock_of_fieldlines(flock, 32, 0.5_dp, field, 0.0_dp, 1.0_dp, 1.0_dp)
    call compute_sqrt_nu_geometry(field, flock, 3.0_dp, 1.0_dp, 1, 1.0_dp, &
                                  nu_star, 33, 2, 0.001_dp, more_lines)
    anchor_32 = more_lines%normalization* &
                sum(more_lines%g_modes(:, 1)*more_lines%iv_eta_b)
    call make_flock_of_fieldlines(flock, 64, 0.5_dp, field, 0.0_dp, 1.0_dp, 1.0_dp)
    call compute_sqrt_nu_geometry(field, flock, 3.0_dp, 1.0_dp, 1, 1.0_dp, &
                                  nu_star, 33, 2, 0.001_dp, most_lines)
    anchor_64 = most_lines%normalization* &
                sum(most_lines%g_modes(:, 1)*most_lines%iv_eta_b)
    call check_refinement(abs(anchor_32 - anchor_16), abs(anchor_64 - anchor_32), &
                          abs(anchor_64), 'fieldlines16/32/64')
    call check_label_independence()
contains
    ! B2 limiting case: for an exactly label-independent QI field every label has
    ! the same tau_b and tau_d, so the per-label mean equals the B1 offset.
    subroutine check_label_independence()
        type(anti_sigma_field_t) :: qi
        type(flock_of_fieldlines_t) :: qi_flock
        type(sqrt_nu_geometry_t) :: geo
        type(sqrt_nu_result_t) :: line_result
        real(dp) :: b1, b2, omega(33), diffusion(33), dIv(33, 2), gsyn(2)
        integer :: k, n

        call qi%anti_sigma_field_init(0.0_dp, 1.0_dp, 1.0_dp, -0.2_dp, 0.0_dp)
        call make_flock_of_fieldlines(qi_flock, 16, 0.5_dp, qi, 0.0_dp, 1.0_dp, 1.0_dp)
        call compute_sqrt_nu_geometry(qi, qi_flock, 3.0_dp, 1.0_dp, 1, 1.0_dp, &
                                      nu_star, 33, 2, 0.001_dp, geo)
        n = size(geo%tau_d_per_nu_line)
        if (any(shape(geo%tau_b_line) /= [33, n])) error stop 'tau_b_line shape'
        if (maxval(abs(sum(geo%tau_b_line, dim=2)/real(n, dp) - geo%tau_b)) > &
            1.e-14_dp*maxval(geo%tau_b)) error stop 'tau_b_line mean is not tau_b'
        if (abs(sum(geo%tau_d_per_nu_line)/real(n, dp) - geo%tau_d_per_nu(1)) > &
            1.e-14_dp*geo%tau_d_per_nu(1)) error stop 'tau_d_line mean is not tau_d'
        omega = 0.05_dp
        ! An exact QI field has zero true drift, so supply a nonzero synthetic dI_V/deta
        ! to make the B1 and B2 offsets nonzero; only the tau tables come from the field.
        dIv(:, 1) = 1.0_dp + geo%eta
        dIv(:, 2) = -0.5_dp*geo%eta
        gsyn = [2.0_dp, 3.0_dp]
        call evaluate_sqrt_nu_modes(geo%eta, geo%tau_b, nu_star(1)*geo%tau_d_per_nu, &
                                    omega, dIv, gsyn, &
                                    0.5_dp, line_result)
        if (line_result%status /= 0) error stop 'label-independent B1 failed'
        b1 = line_result%lambda_off
        b2 = 0.0_dp
        do k = 1, n
            diffusion = nu_star(1)*geo%tau_d_per_nu_line(k)
            call evaluate_sqrt_nu_modes(geo%eta, geo%tau_b_line(:, k), &
                                        diffusion, omega, dIv, gsyn, &
                                        0.5_dp, line_result)
            if (line_result%status /= 0) error stop 'label-independent B2 failed'
            b2 = b2 + line_result%lambda_off/real(n, dp)
        end do
        if (abs(b1) < 1.e-8_dp) error stop 'label-independent offset is vacuously zero'
        print *, 'label-independent B1, B2, rel diff: ', b1, b2, &
            abs(b2 - b1)/max(abs(b1), 1.e-300_dp)
        if (abs(b2 - b1) > 1.e-12_dp*abs(b1)) error stop 'B2 /= B1 for QI field'
    end subroutine check_label_independence

    subroutine check_reference_origin()
        use constants, only: pi
        use field_checks, only: suspect_omnigenous_origin_not_minimum
        type(anti_sigma_field_t) :: minimum_field, maximum_field, qa_field
        type(flock_of_fieldlines_t) :: at_zero, at_reference
        type(sqrt_nu_geometry_t) :: zero_geometry, reference_geometry
        real(dp) :: theta, phi, first, second
        integer :: j

        ! These two analytic fields describe the same local well with opposite
        ! signs of cos(4*phi). The physical maximum_field is never modified.
        call minimum_field%anti_sigma_field_init(0.0_dp, 4.0_dp, &
                                                 1.0_dp, -0.2_dp, 0.0_dp)
        call maximum_field%anti_sigma_field_init(0.0_dp, 4.0_dp, &
                                                 1.0_dp, 0.2_dp, 0.0_dp)
        if (.not. suspect_omnigenous_origin_not_minimum(maximum_field, &
                                                        N_tor=4.0_dp, M_pol=0.0_dp)) &
            error stop 'the unchanged zero-origin check must reject a maximum'
        if (suspect_omnigenous_origin_not_minimum(maximum_field, &
                         N_tor=4.0_dp, M_pol=0.0_dp, reference_origin=[0.0_dp, pi/4])) &
            error stop 'analytic QI reference minimum rejected'
        call qa_field%anti_sigma_field_init(1.0_dp, 0.0_dp, &
                                            1.0_dp, 0.2_dp, 0.0_dp)
        if (suspect_omnigenous_origin_not_minimum(qa_field, &
                           N_tor=0.0_dp, M_pol=1.0_dp, reference_origin=[pi, 0.0_dp])) &
            error stop 'analytic QA reference minimum rejected'
        do j = 0, 12
            theta = real(j, dp)*0.17_dp
            phi = real(j, dp)*0.03_dp
            call minimum_field%compute_B_mod(theta, phi, first)
            call maximum_field%compute_B_mod(theta, phi + pi/4, second)
            if (abs(first - second) > 1.e-14_dp) &
                error stop 'independent cosine well identity failed'
        end do
        call make_flock_of_fieldlines(at_zero, 16, 0.5_dp, minimum_field, &
                                      0.0_dp, 4.0_dp, 4.0_dp)
        call make_flock_of_fieldlines(at_reference, 16, 0.5_dp, maximum_field, &
                                0.0_dp, 4.0_dp, 4.0_dp, reference_origin=[0.0_dp, pi/4])
        if (maxval(abs(at_reference%fieldlines%theta_0 - &
                       at_zero%fieldlines%theta_0)) > 1.e-14_dp) &
            error stop 'QI fieldline label changed under explicit phi reference'
        if (maxval(abs(at_reference%fieldlines%phi_0 - &
                       at_zero%fieldlines%phi_0 - pi/4)) > 1.e-14_dp) &
            error stop 'physical fieldline reference not retained'
        call compute_sqrt_nu_geometry(minimum_field, at_zero, 3.0_dp, &
                             1.0_dp, 1, 1.0_dp, nu_star, 33, 2, 0.001_dp, zero_geometry)
        call compute_sqrt_nu_geometry(maximum_field, at_reference, 3.0_dp, &
                        1.0_dp, 1, 1.0_dp, nu_star, 33, 2, 0.001_dp, reference_geometry)
        if (maxval(abs(zero_geometry%eta - reference_geometry%eta)) > 1.e-12_dp) &
            error stop 'pitch geometry changed under reference choice'
        if (maxval(abs(zero_geometry%tau_b - reference_geometry%tau_b)) > &
            1.e-9_dp*maxval(zero_geometry%tau_b)) &
            error stop 'bounce integrals changed under reference choice'
        if (abs(at_zero%I_ref - at_reference%I_ref) > 1.e-10_dp*at_zero%I_ref) &
            error stop 'boundary diffusion changed under reference choice'
    end subroutine check_reference_origin

    function evaluate_field(data) result(offset)
        type(sqrt_nu_geometry_t), intent(in) :: data
        real(dp) :: offset, omega(size(data%eta))
        type(sqrt_nu_result_t) :: result
        omega = 0.05_dp
        call evaluate_sqrt_nu_modes(data%eta, data%tau_b, &
                                   nu_star(1)*data%tau_d_per_nu, omega, data%dIv_deta, &
                                    data%g_modes(:, 1), data%normalization, result, &
                                    iv_deep=data%iv_deep)
        if (result%status /= 0) error stop 'synthetic geometry evaluation failed'
        offset = result%lambda_off
        if (.not. ieee_is_finite(offset)) error stop 'nonfinite synthetic field result'
    end function evaluate_field

    subroutine check_refinement(first_difference, second_difference, scale, label)
        real(dp), intent(in) :: first_difference, second_difference, scale
        character(len=*), intent(in) :: label
        print *, label, ': successive absolute differences ', &
            first_difference, second_difference
        ! Check numerical contraction on this synthetic case. This does not
        ! establish a production physics acceptance tolerance or WKB domain.
        if (second_difference > max(first_difference, 1.e-8_dp*scale)) &
            error stop 'synthetic refinement did not contract'
    end subroutine check_refinement
end program test_sqrt_nu_field
