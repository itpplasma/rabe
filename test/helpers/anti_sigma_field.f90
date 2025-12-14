module anti_sigma_field
    use constants, only: dp, pi
    use field_base, only: field_t

    implicit none

    type, extends(field_t) :: anti_sigma_field_t
        real(dp) :: M_pol, N_tor
        real(dp) :: B_0, eps_0, eps_1
        real(dp) :: sign
    contains
        procedure :: anti_sigma_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_sqrt_g11
    end type anti_sigma_field_t

contains

    subroutine anti_sigma_field_init(self, M_pol, N_tor, B_0, eps_0, eps_1)
        class(anti_sigma_field_t), intent(out) :: self
        real(dp), intent(in) :: M_pol, N_tor
        real(dp), intent(in) :: B_0, eps_0, eps_1

        self%M_pol = M_pol
        self%N_tor = N_tor
        self%B_0 = B_0
        self%eps_0 = eps_0
        self%eps_1 = eps_1

        if (self%eps_0*self%eps_1 < 0.0_dp) then
            print *, "error anti_sigma_field_init: eps_0 and eps_1 must have same sign!"
            error stop
        else
            self%sign = sign(1.0_dp, eps_0)
        end if
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

        real(dp) :: chi

        call self%compute_B_mod(theta, phi, B_mod)

        chi = self%M_pol*theta - self%N_tor*phi
        dB_dx(1) = 0.0_dp
        dB_dx(2) = self%B_0*( &
                   -self%M_pol*self%eps_0*sin(chi) &
                   + self%sign*self%M_pol*self%eps_1*cos(theta)*sin(chi) &
                   - self%eps_1*sin(theta)*(1 - self%sign*cos(chi)) &
                   )
        dB_dx(3) = self%B_0*( &
                   +self%N_tor*self%eps_0*sin(chi) &
                   - self%N_tor*self%eps_1*cos(theta)*self%sign*sin(chi) &
                   )
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(anti_sigma_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: chi

        chi = self%M_pol*theta - self%N_tor*phi

        B_mod = self%B_0*(1 &
                          + self%eps_0*cos(chi) &
                          + self%eps_1*cos(theta)*(1 - self%sign*cos(chi)) &
                          )
    end subroutine compute_B_mod

    subroutine compute_sqrt_g11(self, theta, phi, sqrt_g11)
        class(anti_sigma_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: sqrt_g11

        sqrt_g11 = 1.0_dp

    end subroutine compute_sqrt_g11

end module anti_sigma_field
