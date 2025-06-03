module misc
    use constants, only: dp, pi
    implicit none

contains

    elemental function S_A(angle)
        real(dp), intent(in) :: angle
        real(dp) :: S_A

        real(dp) :: angle_in_period
        real(dp), parameter :: fit_fac = 0.26_dp

        angle_in_period = modulo(angle + pi, 2.0_dp*pi) - pi

        S_A = fit_fac*(angle_in_period - sign(0.5_dp*pi, angle_in_period))
    end function S_A

    elemental function S_B(angle)
        real(dp), intent(in) :: angle
        real(dp) :: S_B

        real(dp) :: angle_in_period
        real(dp), parameter :: fit_fac = 1.85_dp

        angle_in_period = modulo(angle + pi, 2.0_dp*pi) - pi

        S_B = sign(fit_fac, angle_in_period)
    end function S_B

end module misc
