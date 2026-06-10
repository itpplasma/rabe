module coefficients
    use constants, only: dp, pi, machine_eps
    use fieldline_mod, only: flock_of_fieldlines_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

contains

    function calc_finite_boundary_layer_correction(flock, &
                                                   field, &
                                                   R, &
                                                   dr_dAtheta) &
        result(Lambda_S)
        use surface_average_mod, only: surface_average_t, calc_surface_averages
        use field_base, only: field_t
        type(flock_of_fieldlines_t), intent(in) :: flock
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: R, dr_dAtheta
        real(dp) :: Lambda_S

        real(dp) :: I_ref, M_pol, N_tor, iota, eta_b
        real(dp) :: B_theta_covariant, B_phi_covariant
        real(dp) :: I_ref_hat
        real(dp) :: helical_factor
        type(surface_average_t) :: average

        if (R < machine_eps) then
            print *, "error: major radius R must be positive."
            print *, "R =", R
            error stop
        end if
        if (ieee_is_nan(dr_dAtheta)) then
            print *, "error: dr_dAtheta is NaN."
            error stop
        end if

        call field%get_covariant_components(B_theta_covariant, B_phi_covariant)
        if (ieee_is_nan(B_theta_covariant)) then
            print *, "error: B_theta_covariant is NaN."
            error stop
        end if
        if (ieee_is_nan(B_phi_covariant)) then
            print *, "error: B_phi_covariant is NaN."
            error stop
        end if

        call calc_surface_averages(flock, average)

        M_pol = flock%M_pol
        N_tor = flock%N_tor
        iota = flock%iota
        I_ref = flock%I_ref
        eta_b = flock%eta_b

        helical_factor = (B_phi_covariant*M_pol + B_theta_covariant*N_tor) &
                         /(M_pol*iota - N_tor)

        I_ref_hat = I_ref/eta_b/(2.0_dp*pi*R)
        Lambda_S = sqrt(I_ref_hat)*average%B_squared*eta_b**2/average%lambda_b
        Lambda_S = 1.5_dp*0.855_dp/sqrt(pi)*Lambda_S
        Lambda_S = -helical_factor*dr_dAtheta*Lambda_S
        if (ieee_is_nan(Lambda_S)) then
            print *, "error: Lambda_S is NaN."
            error stop
        end if
    end function calc_finite_boundary_layer_correction

    function calc_gradient_scaling_factor_r_eff(flock, psi_tor_edge, sign_sqrtg) &
        result(dr_dAtheta)
        use surface_average_mod, only: surface_average_t, calc_surface_averages
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: psi_tor_edge
        integer, intent(in) :: sign_sqrtg
        real(dp) :: dr_dAtheta

        type(surface_average_t) :: average

        if (abs(sign_sqrtg) /= 1) then
            print *, "error: sign_sqrtg must be +1 or -1."
            print *, "sign_sqrtg =", sign_sqrtg
            error stop
        end if

        call calc_surface_averages(flock, average)

        if (abs(average%nabla_s*psi_tor_edge) < machine_eps) then
            print *, "error: nabla_s * psi_tor_edge is zero."
            print *, "nabla_s =", average%nabla_s
            print *, "psi_tor_edge =", psi_tor_edge
            error stop
        end if

        dr_dAtheta = real(sign_sqrtg, dp)/(average%nabla_s*psi_tor_edge)

        if (ieee_is_nan(dr_dAtheta)) then
            print *, "error: dr_dAtheta is NaN."
            error stop
        end if
    end function calc_gradient_scaling_factor_r_eff

    subroutine calc_offset_coefficients(flock, R, dr_dAtheta, &
                                        Lambda_A, Lambda_B)
        use deviation, only: calc_deviation
        use field_base, only: field_t
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: R, dr_dAtheta
        real(dp), intent(out) :: Lambda_A, Lambda_B

        real(dp) :: deviation_A, deviation_B

        call calc_deviation(flock, deviation_A, deviation_B)

        Lambda_A = deviation_A*dr_dAtheta*sqrt(0.5_dp*R*pi)
        Lambda_B = deviation_B*0.5_dp*R*pi*dr_dAtheta
    end subroutine calc_offset_coefficients

    function calc_nu_star_crit(flock, R) result(nu_star_crit)
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: R
        real(dp) :: nu_star_crit

        real(dp) :: eta_b, I_ref_hat, max_delta_eta

        if (R < machine_eps) then
            print *, "error: major radius R must be positive."
            print *, "R =", R
            error stop
        end if

        eta_b = flock%eta_b
        I_ref_hat = flock%I_ref/eta_b/(2.0_dp*pi*R)
        max_delta_eta = eta_b - 1.0_dp/minval(flock%fieldlines%B_max(1))

        nu_star_crit = 0.125_dp*(max_delta_eta/eta_b)**2.0_dp/I_ref_hat

        if (ieee_is_nan(nu_star_crit)) then
            print *, "error: Lambda_S is NaN."
            error stop
        end if
    end function calc_nu_star_crit

end module coefficients
