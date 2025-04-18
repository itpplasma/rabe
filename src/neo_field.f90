module neo_field
    use field_base, only: field_t
    use neo_magfie, only: neo_magfie_a

    implicit none
    integer, parameter :: dp = kind(1.0d0)

    type, extends(field_t) :: neo_field_t
    contains
        procedure :: neo_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: neo_change_stor
    end type neo_field_t

contains

    subroutine neo_field_init(self, bc_filename, stor)
        use neo_magfie, only: magfie_newspline

        class(neo_field_t), intent(out) :: self
        character(*), intent(in) :: bc_filename
        real(dp), intent(in) :: stor

        real(dp) :: x(3), dummy_1, dummy_2, dummy_3(3)

        if (magfie_newspline .ne. 1) then
            print *, "There can only be one neo_field "// &
                "and only for one fluxsurface at a time! "// &
                "Change the fluxsurface with neo_change_stor!"
            error stop
        end if

        x = (/stor, 0.0_dp, 0.0_dp/)
        call neo_magfie_a(x, dummy_1, dummy_2, dummy_3)

    end subroutine neo_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(neo_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        real(dp) :: x(3)

        x = (/0.0_dp, theta, phi/)
        call neo_magfie_a(x, B_mod, sqrtg, dB_dx)

    end subroutine compute_B_sqrtg_dB_dx

    subroutine neo_change_stor(self, stor)
        use neo_magfie, only: magfie_newspline
        use neo_input, only: es

        class(neo_field_t), intent(in) :: self
        real(dp), intent(in) :: stor

        real(dp) :: x(3), dummy_1, dummy_2, dummy_3(3)

        if (magfie_newspline .ne. 1) magfie_newspline = 1
        x = (/stor, 0.0_dp, 0.0_dp/)
        call neo_magfie_a(x, dummy_1, dummy_2, dummy_3)
    end subroutine

end module neo_field
