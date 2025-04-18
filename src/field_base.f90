module field_base
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    type, abstract :: field_t
    contains
        procedure(compute_B_sqrtg_dB_dx), deferred :: compute_B_sqrtg_dB_dx
    end type field_t

    interface
        subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
            import :: field_t, dp
            class(field_t), intent(in) :: self
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)
        end subroutine
    end interface

end module field_base
