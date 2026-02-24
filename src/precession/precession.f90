module precession
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

contains

    subroutine compute_precession_correction(field, l_c, Omega_hat, correction)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: l_c
        real(dp), intent(in) :: Omega_hat
        real(dp), intent(out) :: correction

        correction = 0.0_dp

    end subroutine compute_precession_correction

end module precession
