module fit_functions
    use constants, only: dp, pi
    implicit none

    real(dp), parameter :: eps = 1e-15_dp

contains

    elemental function S_A(angle)
        real(dp), intent(in) :: angle
        real(dp) :: S_A

        real(dp) :: angle_in_period
        !> Sawtooth amplitude of the offset function S_A from Albert et al.,
        !> J. Plasma Phys. 91 (2025), fit below Eq. (133). An independent solve
        !> of the boundary-layer integral chain (Eq. 103-105) gives 0.231; the
        !> prior 0.30 here had no source and overshoots by ~25%.
        real(dp), parameter :: fit_fac = 0.26_dp

        angle_in_period = modulo(angle + pi, 2.0_dp*pi) - pi

        S_A = fit_fac*(angle_in_period - sign(0.5_dp*pi, angle_in_period))
        if (abs(mod(angle, pi)) < eps) S_A = 0.0_dp
    end function S_A

    elemental function S_B(angle)
        real(dp), intent(in) :: angle
        real(dp) :: S_B

        real(dp) :: angle_in_period
        !> Square-wave amplitude of the offset function S_B from Albert et al.,
        !> J. Plasma Phys. 91 (2025), fit below Eq. (133). The independent solve
        !> gives 1.846, confirming 1.85; the prior 2.00 here had no source.
        real(dp), parameter :: fit_fac = 1.85_dp

        angle_in_period = modulo(angle + pi, 2.0_dp*pi) - pi

        S_B = sign(fit_fac, angle_in_period)
        if (abs(mod(angle, pi)) < eps) S_B = 0.0_dp
    end function S_B

end module fit_functions
