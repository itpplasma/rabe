module surface_average_mod
    use constants, only: dp, pi
    implicit none

    type :: surface_average_t
        real(dp) :: normalization
        real(dp) :: B_squared
        real(dp) :: lambda_b
        real(dp) :: sqrt_g11
    end type surface_average_t

contains

    subroutine calc_surface_averages(fieldlines, surface_average)
        use fieldline_mod, only: fieldline_t
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(surface_average_t), intent(out) :: surface_average

        real(dp), dimension(size(fieldlines)) :: well_lengths

        integer :: n_fieldlines
        real(dp) :: dxi_0

        n_fieldlines = size(fieldlines)
        dxi_0 = (fieldlines(n_fieldlines)%xi_0 - fieldlines(1)%xi_0)/ &
                real(n_fieldlines - 1, kind=dp)

        well_lengths = fieldlines%phi_max(2) - fieldlines%phi_max(1)
        surface_average%normalization = sum(fieldlines%integral_one_over_B_squared)* &
                                        dxi_0
        surface_average%B_squared = sum(well_lengths)*dxi_0/ &
                                    surface_average%normalization
        surface_average%lambda_b = sum(fieldlines%integral_lambda_b_over_B_squared)* &
                                   dxi_0/surface_average%normalization
        surface_average%sqrt_g11 = sum(fieldlines%integral_sqrt_g11_over_B_squared)* &
                                   dxi_0/surface_average%normalization
    end subroutine calc_surface_averages

end module surface_average_mod
