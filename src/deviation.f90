module deviation
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

    type :: surface_average_t
        real(dp) :: normalization
        real(dp) :: B_squared
        real(dp) :: lambda_b
    end type surface_average_t

contains

    subroutine calc_deviation(fieldlines, deviation_A, deviation_B)
        use fieldline_integrals, only: fieldline_modes_t
        use fieldline_integrals, only: allocate_modes
        use fieldline_integrals, only: fourier_transform_over_label
        use misc, only: S_A, S_B

        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(out) :: deviation_A, deviation_B

        real(dp), parameter :: tol = 1e-12

        type(fieldline_modes_t) :: modes
        real(dp) :: iota_p, eta_b
        type(surface_average_t) :: average
        real(dp) :: symmetric_remainder
        real(dp) :: B_squared_sqrtg

        logical :: any_has_sin_part

        call fourier_transform_over_label(fieldlines, modes)

        any_has_sin_part = .false.
        if (has_sin_modes(modes%delta_aspect_ratio)) then
            print *, "error: non-vanishing sin part of delta aspect ratio: "
            print *, "sin part: ", sum(abs(modes%delta_aspect_ratio%sin_coeffs))
            print *, "cos part: ", sum(abs(modes%delta_aspect_ratio%cos_coeffs))
            any_has_sin_part = .true.
        end if
        if (has_sin_modes(modes%delta_eta)) then
            print *, "error: non-vanishing sin part of delta eta: "
            print *, "sin part: ", sum(abs(modes%delta_eta%sin_coeffs))
            print *, "cos part: ", sum(abs(modes%delta_eta%cos_coeffs))
            any_has_sin_part = .true.
        end if
        if (any_has_sin_part) error stop

        call calc_surface_averages(fieldlines, average)

        iota_p = fieldlines(1)%iota_p
        eta_b = fieldlines(1)%eta_b

        deviation_A = pi*sum(modes%radial_drift%sin_coeffs* &
                             modes%delta_aspect_ratio%cos_coeffs* &
                             S_A(iota_p*modes%delta_aspect_ratio%mode_numbers))

        symmetric_remainder = pi*sum(modes%radial_drift%cos_coeffs* &
                                     modes%delta_aspect_ratio%sin_coeffs* &
                                     S_A(iota_p*modes%delta_aspect_ratio%mode_numbers))

        if (abs(symmetric_remainder/deviation_A) > tol) then
            print *, "warning: non-vanishing symmetric part of deviation A: "
            print *, "symmetric: ", symmetric_remainder
            print *, "antisymmetric: ", deviation_A
            print *, "ratio: ", symmetric_remainder/deviation_A
        end if

        deviation_A = deviation_A*average%B_squared/average%lambda_b* &
                      sqrt(eta_b)*sqrt(fieldlines(1)%I_ref)/average%normalization

        deviation_B = pi*sum(modes%radial_drift%sin_coeffs* &
                             modes%delta_eta%cos_coeffs* &
                             S_B(iota_p*modes%delta_eta%mode_numbers))

        symmetric_remainder = pi*sum(modes%radial_drift%cos_coeffs* &
                                     modes%delta_eta%sin_coeffs* &
                                     S_B(iota_p*modes%delta_eta%mode_numbers))

        if (abs(symmetric_remainder/deviation_B) > tol) then
            print *, "warning: non-vanishing symmetric part of deviation B: "
            print *, "symmetric: ", symmetric_remainder
            print *, "antisymmetric: ", deviation_B
            print *, "ratio: ", symmetric_remainder/deviation_B
        end if

        deviation_B = deviation_B*average%B_squared/average%lambda_b*0.5_dp/ &
                      average%normalization

    end subroutine calc_deviation

    subroutine calc_surface_averages(fieldlines, surface_average)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(surface_average_t), intent(out) :: surface_average

        real(dp), dimension(size(fieldlines)) :: well_lengths

        integer :: n_fieldlines
        real(dp) :: dtheta_0

        n_fieldlines = size(fieldlines)
        dtheta_0 = (fieldlines(n_fieldlines)%theta_0 - fieldlines(1)%theta_0)/ &
                   real(n_fieldlines - 1, kind=dp)

        well_lengths = fieldlines%phi_max(2) - fieldlines%phi_max(1)
        surface_average%normalization = sum(fieldlines%integral_one_over_B_squared)* &
                                        dtheta_0
        surface_average%B_squared = sum(well_lengths)*dtheta_0/ &
                                    surface_average%normalization
        surface_average%lambda_b = sum(fieldlines%integral_lambda_b_over_B_squared)* &
                                   dtheta_0/surface_average%normalization
    end subroutine calc_surface_averages

    function has_sin_modes(modes)
        use fieldline_integrals, only: modes_t
        type(modes_t), intent(in) :: modes
        logical :: has_sin_modes

        real(dp), parameter :: tol = 1e-3, numerical_zero = 1e-8
        real(dp) :: sum_sin, sum_cos

        sum_sin = sum(abs(modes%sin_coeffs))
        sum_cos = sum(abs(modes%cos_coeffs))

        if (sum_cos > numerical_zero) then
            has_sin_modes = sum_sin/sum_cos > tol
        else
            has_sin_modes = sum_sin > numerical_zero
        end if
    end function has_sin_modes

    function get_B_squared_sqrtg_psi_edge(field) result(B_squared_sqrtg_psi_edge)
        use field_base, only: field_t

        class(field_t), intent(in) :: field

        real(dp) :: B_squared_sqrtg_psi_edge
        real(dp), parameter :: theta = 0.0_dp, phi = 0.0_dp
        real(dp) :: sqrtg, B_mod, dummy_dB_dx(3)

        call field%compute_B_sqrtg_dB_dx(theta, phi, B_mod, sqrtg, dummy_dB_dx)
        B_squared_sqrtg_psi_edge = abs(B_mod**2.0_dp*sqrtg)
    end function get_B_squared_sqrtg_psi_edge

end module deviation
