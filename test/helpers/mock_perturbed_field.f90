module mock_perturbed_field
    use constants, only: dp
    use field_base, only: field_t
    use mock_field, only: mock_field_t

    implicit none

    type, extends(field_t) :: mock_perturbed_field_t
        class(field_t), allocatable :: field
        type(mock_field_t) :: pert_field
    contains
        procedure :: mock_perturbed_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
    end type mock_perturbed_field_t

contains

    subroutine mock_perturbed_field_init(self, field, theta_mode, phi_mode, B_amplitude)
        class(mock_perturbed_field_t), intent(out) :: self
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta_mode, phi_mode
        real(dp), intent(in) :: B_amplitude

        allocate (self%field, source=field)
        call self%pert_field%mock_field_init(theta_mode, phi_mode, 0.0_dp, B_amplitude)
    end subroutine mock_perturbed_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(mock_perturbed_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "The type mock_perturbed_field does not provide sqrtg."// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(mock_perturbed_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: B_mod_pert, dB_dx_pert(3)

        call self%field%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)
        call self%pert_field%compute_B_and_dB_dx(theta, phi, B_mod_pert, dB_dx_pert)
        B_mod = B_mod + B_mod_pert
        dB_dx = dB_dx + dB_dx_pert

        call self%compute_B_mod(theta, phi, B_mod)

    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(mock_perturbed_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: B_pert_mod

        call self%field%compute_B_mod(theta, phi, B_mod)
        call self%pert_field%compute_B_mod(theta, phi, B_pert_mod)
        B_mod = B_mod + B_pert_mod
    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(mock_perturbed_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        nabla_s = 1.0_dp

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(mock_perturbed_field_t), intent(in) :: self

        rel_accuracy_B = self%field%rel_accuracy_B()
        rel_accuracy_B = rel_accuracy_B + 1e-14_dp
    end function rel_accuracy_B

end module mock_perturbed_field
