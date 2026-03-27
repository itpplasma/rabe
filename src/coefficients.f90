module coefficients
    use constants, only: dp, pi
    use fieldline_mod, only: fieldline_t
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

        call calc_surface_averages(fieldlines, average)

        M_pol = fieldlines(1)%M_pol
        N_tor = fieldlines(1)%N_tor
        iota = fieldlines(1)%iota
        I_ref = fieldlines(1)%I_ref

        helical_factor = (B_phi_covariant*M_pol + B_theta_covariant*N_tor) &
                         /(M_pol*iota - N_tor)
        covariant_factor = (B_phi_covariant*M_pol + B_theta_covariant*N_tor)

        I_ref_hat = I_ref/eta_b*covariant_factor/R
        Lambda_finite = sqrt(I_ref_hat)*average%B_squared*eta_b**2/average%lambda_b
        Lambda_finite = 0.75_dp*0.855_dp/pi*sqrt(2.0_dp)*Lambda_finite
        Lambda_finite = -helical_factor*dr_dAtheta*Lambda_finite
    end function calc_finite_boundary_layer_correction

end module coefficients
