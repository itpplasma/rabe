module precession
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

contains

    subroutine compute_precession_correction(field, l_c, omega_hat, correction)
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(in) :: l_c
        real(dp), intent(in) :: omega_hat
        real(dp), intent(out) :: correction

        correction = 0.0_dp

    end subroutine compute_precession_correction

end module precession
