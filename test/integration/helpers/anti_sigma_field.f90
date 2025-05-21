module anti_sigma_field
    use constants, only: dp, pi
    use field_base, only: field_t

    implicit none

    type, extends(field_t) :: anti_sigma_field_t
        real(dp) :: N_tor
        real(dp) :: B_0, eps_0, eps_1

        real(dp) :: av_B2_over_av_lambda_b
    contains
        procedure :: anti_sigma_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
    end type anti_sigma_field_t

contains

    subroutine anti_sigma_field_init(self, N_tor, B_0, eps_0, eps_1)
        class(anti_sigma_field_t), intent(out) :: self
        real(dp), intent(in) :: N_tor
        real(dp), intent(in) :: B_0, eps_0, eps_1

        self%N_tor = N_tor
        self%B_0 = B_0
        self%eps_0 = eps_0
        self%eps_1 = eps_1

        self%av_B2_over_av_lambda_b = B_0**2*pi/sqrt(8*eps_0)
    end subroutine anti_sigma_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(anti_sigma_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "The type anti_sigma_field does not provide sqrtg."// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(anti_sigma_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        call self%compute_B_mod(theta, phi, B_mod)
        dB_dx(1) = 0.0_dp
        dB_dx(2) = -self%B_0*self%eps_1*sin(theta)*(1 - cos(self%N_tor*phi))
        dB_dx(3) = -self%B_0*(self%N_tor*self%eps_0*sin(self%N_tor*phi) &
                              - self%N_tor*self%eps_1*cos(theta)*sin(self%N_tor*phi))
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(anti_sigma_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        B_mod = self%B_0*(1 + self%eps_0*cos(self%N_tor*phi) &
                          + self%eps_1*cos(theta)*(1 - cos(self%N_tor*phi)))
    end subroutine compute_B_mod

end module anti_sigma_field
