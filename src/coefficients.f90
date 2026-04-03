module coefficients
    use constants, only: dp, pi, machine_eps
    use fieldline_mod, only: fieldline_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

contains

    function calc_finite_boundary_layer_correction(fieldlines, &
                                                   R, &
                                                   dr_dAtheta, &
                                                   B_theta_covariant, &
                                                   B_phi_covariant) &
        result(Lambda_finite)
        use surface_average_mod, only: surface_average_t, calc_surface_averages
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: R, dr_dAtheta
        real(dp), intent(in) :: B_theta_covariant, B_phi_covariant
        real(dp) :: Lambda_finite

        real(dp) :: I_ref, M_pol, N_tor, iota, eta_b
        real(dp) :: I_ref_hat
        real(dp) :: helical_factor, covariant_factor
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
        if (ieee_is_nan(B_theta_covariant)) then
            print *, "error: B_theta_covariant is NaN."
            error stop
        end if
        if (ieee_is_nan(B_phi_covariant)) then
            print *, "error: B_phi_covariant is NaN."
            error stop
        end if

        call calc_surface_averages(fieldlines, average)

        M_pol = fieldlines(1)%M_pol
        N_tor = fieldlines(1)%N_tor
        iota = fieldlines(1)%iota
        I_ref = fieldlines(1)%I_ref
        eta_b = fieldlines(1)%eta_b

        covariant_factor = (B_phi_covariant + B_theta_covariant*iota)
        if (covariant_factor < machine_eps) then
            print *, "error: covariant factor must be positive."
            print *, "B_phi_covariant =", B_phi_covariant
            print *, "B_theta_covariant =", B_theta_covariant
            print *, "iota =", iota
            error stop
        end if

        helical_factor = (B_phi_covariant*M_pol + B_theta_covariant*N_tor) &
                         /(M_pol*iota - N_tor)

        I_ref_hat = I_ref/eta_b*covariant_factor/(2.0_dp*pi*R)
        Lambda_finite = sqrt(I_ref_hat)*average%B_squared*eta_b**2/average%lambda_b
        Lambda_finite = 1.5_dp*0.855_dp/sqrt(pi)*Lambda_finite
        Lambda_finite = -helical_factor*dr_dAtheta*Lambda_finite
        if (ieee_is_nan(Lambda_finite)) then
            print *, "error: Lambda_finite is NaN."
            error stop
        end if
    end function calc_finite_boundary_layer_correction

    function calc_nu_star_crit(fieldlines, &
                               R, &
                               B_theta_covariant, &
                               B_phi_covariant) result(nu_star_crit)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: R
        real(dp), intent(in) :: B_theta_covariant, B_phi_covariant
        real(dp) :: nu_star_crit

        real(dp) :: eta_b, iota, I_ref_hat, max_delta_eta
        real(dp) :: covariant_factor

        if (R < machine_eps) then
            print *, "error: major radius R must be positive."
            print *, "R =", R
            error stop
        end if
        if (ieee_is_nan(B_theta_covariant)) then
            print *, "error: B_theta_covariant is NaN."
            error stop
        end if
        if (ieee_is_nan(B_phi_covariant)) then
            print *, "error: B_phi_covariant is NaN."
            error stop
        end if

        eta_b = fieldlines(1)%eta_b
        iota = fieldlines(1)%iota
        covariant_factor = (B_phi_covariant + B_theta_covariant*iota)
        if (covariant_factor < machine_eps) then
            print *, "error: covariant factor must be positive."
            print *, "B_phi_covariant =", B_phi_covariant
            print *, "B_theta_covariant =", B_theta_covariant
            print *, "iota =", iota
            error stop
        end if
        I_ref_hat = 1.0_dp/(2.0_dp*pi*R*eta_b)*fieldlines(1)%I_ref*covariant_factor
        max_delta_eta = eta_b - 1.0_dp/minval(fieldlines%B_max(1))

        nu_star_crit = 0.125_dp*(max_delta_eta/eta_b)**2.0_dp/I_ref_hat

        if (ieee_is_nan(nu_star_crit)) then
            print *, "error: Lambda_finite is NaN."
            error stop
        end if
    end function calc_nu_star_crit

end module coefficients
