module mock_pert_field
    use constants, only: dp
    use field_base, only: field_t
    use mock_field, only: mock_field_t

    implicit none

    type, extends(field_t) :: mock_pert_field_t
        real(dp) :: theta_mode, phi_mode
        real(dp) :: B_pert
        type(mock_field_t) :: mock_field
    contains
        procedure :: mock_pert_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_mod
    end type mock_pert_field_t

contains

    subroutine mock_pert_field_init(self, theta_mode, phi_mode, B_pert)
        class(mock_pert_field_t), intent(out) :: self
        real(dp), intent(in) :: theta_mode, phi_mode
        real(dp), intent(in) :: B_pert

        self%theta_mode = theta_mode
        self%phi_mode = phi_mode
        self%B_pert = B_pert
    end subroutine mock_pert_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(mock_pert_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "The type mock_pert_field does only provide B_mod."// &
            "Use compute_B_mod instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(mock_pert_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        call self%mock_field%compute_B_mod(theta, phi, B_mod)
        B_mod = B_mod + self%B_pert*cos(self%theta_mode*theta - self%phi_mode*phi)
    end subroutine compute_B_mod

end module mock_pert_field
