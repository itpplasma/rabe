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
contains
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
