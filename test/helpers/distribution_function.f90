module distribution_function
    use constants, only: dp, pi
    use fieldline_mod, only: fieldline_t
    use fieldline_integrals, only: fourier_transform_over_label
    use fieldline_integrals, only: modes_t, allocate_modes
    use fieldline_integrals, only: fieldline_modes_t
    use deviation, only: surface_average_t, calc_surface_averages

    implicit none

contains

    function get_offset_from_distribution(g_modes, fieldlines) result(offset)
        type(modes_t), intent(in) :: g_modes
        type(fieldline_t), dimension(:), intent(in) :: fieldlines

        real(dp) :: offset

        type(fieldline_modes_t) :: modes
        type(surface_average_t) :: averages
        integer :: n_modes

        call fourier_transform_over_label(fieldlines, modes)
        call calc_surface_averages(fieldlines, averages)
        n_modes = size(modes%radial_drift%sin_coeffs)

        offset = pi*sum(modes%radial_drift%sin_coeffs* &
                        g_modes%sin_coeffs(1:n_modes))/ &
                 averages%normalization
    end function get_offset_from_distribution

    subroutine get_g_modes_from_fieldlines(fieldlines, field, nu_star, g_off_modes)
        use misc, only: S_A, S_B
        use neo_field, only: neo_field_t
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(neo_field_t), intent(in) :: field
        real(dp), intent(in) :: nu_star
        type(modes_t), intent(out) :: g_off_modes

        type(fieldline_modes_t) :: modes
        type(surface_average_t) :: averages

        real(dp) :: covariant_factor, prefactor_A
        real(dp) :: l_c

        integer :: max_mode
        real(dp), dimension(:), allocatable :: theta_mid, g_off

        character(len=1024) :: label

        call fourier_transform_over_label(fieldlines, modes)
        call calc_surface_averages(fieldlines, averages)

        max_mode = size(modes%delta_eta%cos_coeffs, dim=1)
        call allocate_modes(g_off_modes, max_mode)

        l_c = field%R*pi/(2.0_dp*nu_star)

        covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*field%iota)
        prefactor_A = 2.0_dp*sqrt(covariant_factor*fieldlines(1)%eta_b* &
                                  fieldlines(1)%I_ref/l_c)
        g_off_modes%sin_coeffs = prefactor_A*modes%delta_aspect_ratio%cos_coeffs* &
                         S_A(fieldlines(1)%iota_p*modes%delta_aspect_ratio%mode_numbers)
        g_off_modes%sin_coeffs = g_off_modes%sin_coeffs + &
                                 modes%delta_eta%cos_coeffs* &
                                 S_B(fieldlines(1)%iota_p*modes%delta_eta%mode_numbers)
        g_off_modes%sin_coeffs = g_off_modes%sin_coeffs* &
                                 averages%B_squared/averages%lambda_b* &
                                 l_c*0.5_dp
    end subroutine get_g_modes_from_fieldlines

    subroutine get_modes(x, y, modes)
        use fourier, only: real_ft
        real(dp), dimension(:), intent(in) :: x, y
        type(modes_t), intent(out) :: modes

        integer :: n_modes

        n_modes = size(x)/2 + 1

        call allocate_modes(modes, n_modes)

        call real_ft(x, y, modes%cos_coeffs, modes%sin_coeffs)
    end subroutine get_modes

end module distribution_function
