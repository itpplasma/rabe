module fieldline_integrals
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none

    type :: modes_t
        real(dp), dimension(:), allocatable :: cos_coeffs, sin_coeffs
        real(dp), dimension(:), allocatable :: mode_numbers
    end type modes_t

    type :: fieldline_modes_t
        type(modes_t) :: radial_drift
        type(modes_t) :: delta_eta
        type(modes_t) :: delta_aspect_ratio
        type(modes_t) :: g_off
    end type fieldline_modes_t

contains

    subroutine fourier_transform_over_label(field, fieldlines, fieldline_modes)
        use fourier, only: real_ft

        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        type(fieldline_modes_t), intent(out) :: fieldline_modes

        integer :: n_fieldlines, n_modes
        integer :: current
        real(dp) :: average_delta_aspect_ratio
        real(dp), dimension(size(fieldlines)) :: shifted_label

        n_fieldlines = size(fieldlines)
        n_modes = n_fieldlines/2 + 1

        do current = 1, n_fieldlines
            call calc_fieldline_integrals(field, fieldlines(current))
        end do

        ! I_ref can be chosen to be any e.g. I = I_1
        ! (I_ref/I_j)**0.5 - 1 = (max(I_j)/I_j)**0.5 -1 =
        ! ((I+delta)/(I+delta_j))**0.5 -1 ~ 0.5*(delta/I - delta_j/I)
        ! and the result in linear order only differs by a constant delta/I
        ! which does not enter the offset formula
        fieldlines(:)%delta_aspect_ratio = sqrt( &
                                fieldlines(1)%integral_lambda_b_over_B_squared/ &
                                fieldlines(:)%integral_lambda_b_over_B_squared &
                                           ) - 1
        ! average of delta_aspect ratio also does not enter offset formula
        ! can be set it to zero
        average_delta_aspect_ratio = sum(fieldlines(:)%delta_aspect_ratio)/n_fieldlines
        fieldlines(:)%delta_aspect_ratio = fieldlines(:)%delta_aspect_ratio - &
                                           average_delta_aspect_ratio

        call allocate_modes(fieldline_modes%radial_drift, n_modes)
        call allocate_modes(fieldline_modes%delta_aspect_ratio, n_modes)
        call allocate_modes(fieldline_modes%delta_eta, n_modes)

        call real_ft(fieldlines%theta_0, &
                     fieldlines%radial_drift, &
                     fieldline_modes%radial_drift%cos_coeffs, &
                     fieldline_modes%radial_drift%sin_coeffs)

        call real_ft(fieldlines%theta_0, &
                     fieldlines%delta_aspect_ratio, &
                     fieldline_modes%delta_aspect_ratio%cos_coeffs, &
                     fieldline_modes%delta_aspect_ratio%sin_coeffs)

        shifted_label = fieldlines%theta_0 - fieldlines%iota_p
        call real_ft(shifted_label, &
                     fieldlines%delta_eta, &
                     fieldline_modes%delta_eta%cos_coeffs, &
                     fieldline_modes%delta_eta%sin_coeffs)

    end subroutine fourier_transform_over_label

    subroutine calc_fieldline_integrals(field, fieldline)
        use integrate, only: integrate_1d_substituted
        type(fieldline_t), intent(inout) :: fieldline
        class(field_t), intent(in) :: field

        call integrate_1d_substituted(wrapper_local_radial_drift, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%radial_drift)

        call integrate_1d_substituted(wrapper_lambda_over_B_squared, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%integral_lambda_b_over_B_squared)

        call integrate_1d_substituted(wrapper_one_over_B_squared, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%integral_one_over_B_squared)

    contains

        function wrapper_lambda_over_B_squared(phi)
            use fieldline_integrands, only: lambda_over_B_squared

            real(dp), intent(in) :: phi
            real(dp) :: wrapper_lambda_over_B_squared

            real(dp) :: theta

            theta = fieldline%get_theta(phi)
            wrapper_lambda_over_B_squared = lambda_over_B_squared(field, &
                                                                  theta, &
                                                                  phi, &
                                                                  fieldline%eta_b)
        end function wrapper_lambda_over_B_squared

        function wrapper_local_radial_drift(phi)
            use fieldline_integrands, only: local_radial_drift

            real(dp), intent(in) :: phi
            real(dp) :: wrapper_local_radial_drift

            real(dp) :: theta

            theta = fieldline%get_theta(phi)
            wrapper_local_radial_drift = local_radial_drift(field, &
                                                            theta, &
                                                            phi, &
                                                            fieldline%eta_b)
        end function wrapper_local_radial_drift

        function wrapper_one_over_B_squared(phi)
            use fieldline_integrands, only: B_squared

            real(dp), intent(in) :: phi
            real(dp) :: wrapper_one_over_B_squared

            real(dp) :: theta

            theta = fieldline%get_theta(phi)
            wrapper_one_over_B_squared = 1.0_dp/B_squared(field, theta, phi)
        end function wrapper_one_over_B_squared

    end subroutine calc_fieldline_integrals

    subroutine allocate_modes(modes, n_modes)
        integer, intent(in) :: n_modes
        type(modes_t), intent(out) :: modes

        integer :: j

        allocate (modes%cos_coeffs(n_modes))
        allocate (modes%sin_coeffs(n_modes))
        allocate (modes%mode_numbers(n_modes))

        do j = 0, n_modes - 1
            modes%mode_numbers(j + 1) = j
        end do
    end subroutine allocate_modes

end module fieldline_integrals
