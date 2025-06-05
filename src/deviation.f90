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

    subroutine calc_surface_averages(fieldlines, surface_average)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(surface_average_t), intent(out) :: surface_average

        real(dp), dimension(size(fieldlines)) :: well_lengths

        well_lengths = fieldlines%phi_max(2) - fieldlines%phi_max(1)
        surface_average%normalization = sum(fieldlines%integral_one_over_B_squared)
        surface_average%B_squared = sum(well_lengths)/surface_average%normalization
        surface_average%lambda_b = sum(fieldlines%integral_lambda_b_over_B_squared)/ &
                                   surface_average%normalization
    end subroutine calc_surface_averages

end module deviation
