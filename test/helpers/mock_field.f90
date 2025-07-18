module mock_field
    use constants, only: dp
    use field_base, only: field_t

    implicit none

    type, extends(field_t) :: mock_field_t
        real(dp) :: theta_mode, phi_mode
        real(dp) :: B_0, B_amplitude
    contains
        procedure :: mock_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
    end type mock_field_t

contains

    subroutine mock_field_init(self, theta_mode, phi_mode, B_0, B_amplitude)
        class(mock_field_t), intent(out) :: self
        real(dp), intent(in) :: theta_mode, phi_mode
        real(dp), intent(in) :: B_0, B_amplitude

        self%theta_mode = theta_mode
        self%phi_mode = phi_mode
        self%B_0 = B_0
        self%B_amplitude = B_amplitude
    end subroutine mock_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(mock_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "The type mock_field does not provide sqrtg."// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(mock_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        call self%compute_B_mod(theta, phi, B_mod)
        dB_dx(1) = 0.0_dp
        dB_dx(2) = -self%B_amplitude*self%theta_mode &
                   *sin(self%theta_mode*theta - self%phi_mode*phi)
        dB_dx(3) = +self%B_amplitude*self%phi_mode &
                   *sin(self%theta_mode*theta - self%phi_mode*phi)

    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(mock_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        B_mod = self%B_0 + &
                self%B_amplitude*cos(self%theta_mode*theta - self%phi_mode*phi)
    end subroutine compute_B_mod

end module mock_field
