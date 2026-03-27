module surface_average_mod
    use constants, only: dp, pi, machine_eps
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
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(surface_average_t), intent(out) :: surface_average

        real(dp), dimension(size(fieldlines)) :: well_lengths

        integer :: n_fieldlines
        real(dp) :: dxi_0

        n_fieldlines = size(fieldlines)
        if (n_fieldlines < 2) then
   print *, "error: at least two fieldlines are required to calculate surface averages."
            error stop
        end if

        dxi_0 = (fieldlines(n_fieldlines)%xi_0 - fieldlines(1)%xi_0)/ &
                real(n_fieldlines - 1, kind=dp)

        well_lengths = fieldlines%phi_max(2) - fieldlines%phi_max(1)
        surface_average%normalization = sum(fieldlines%integral_one_over_B_squared)* &
                                        dxi_0
        if (abs(surface_average%normalization) < machine_eps) then
            print *, "error: surface average normalization is zero."
            error stop
        end if
        surface_average%B_squared = sum(well_lengths)*dxi_0/ &
                                    surface_average%normalization
        surface_average%lambda_b = sum(fieldlines%integral_lambda_b_over_B_squared)* &
                                   dxi_0/surface_average%normalization
        surface_average%sqrt_g11 = sum(fieldlines%integral_sqrt_g11_over_B_squared)* &
                                   dxi_0/surface_average%normalization

        if (surface_average%lambda_b < machine_eps) then
            print *, "error: average lambda_b <= 0."
            error stop
        end if
        if (surface_average%B_squared < machine_eps) then
            print *, "error: average B_squared <= 0."
            error stop
        end if
        if (ieee_is_nan(surface_average%sqrt_g11)) then
            print *, "error: average sqrt_g11 is NaN."
            error stop
        end if
    end subroutine calc_surface_averages

end module surface_average_mod
