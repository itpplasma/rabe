module neo_field
    use constants, only: dp
    use field_base, only: field_t
    use neo_magfie, only: neo_magfie_a

    implicit none

    type, extends(field_t) :: neo_field_t
        real(dp) :: iota
    contains
        procedure :: neo_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: neo_change_stor
    end type neo_field_t

contains

    subroutine neo_field_init(self, bc_filename, stor)
        use neo_magfie, only: magfie_newspline
        use neo_control, only: in_file

        class(neo_field_t), intent(out) :: self
        character(*), intent(in) :: bc_filename
        real(dp), intent(in) :: stor

        real(dp) :: x(3), dummy_B_mod, dummy_sqrtg, dummy_dB_dx(3)

        if (magfie_newspline .ne. 1) then
            print *, "There can only be one neo_field "// &
                "and only for one fluxsurface at a time! "// &
                "Change the fluxsurface with neo_change_stor!"
            error stop
        end if

        in_file = bc_filename
        x = (/stor, 0.0_dp, 0.0_dp/)
        call neo_magfie_a(x, dummy_B_mod, dummy_sqrtg, dummy_dB_dx, self%iota)

    end subroutine neo_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(neo_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        real(dp) :: x(3), dB_dx_neo2(3), dummy_iota

        x = (/0.0_dp, phi, theta/) !neo convention of x=(r, phi, theta)
        call neo_magfie_a(x, B_mod, sqrtg, dB_dx_neo2, dummy_iota)
        dB_dx(1) = dB_dx_neo2(1)
        dB_dx(2) = dB_dx_neo2(3) !neo convention of x=(r, phi, theta)
        dB_dx(3) = dB_dx_neo2(2) !neo convention of x=(r, phi, theta)

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(neo_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: x(3), dB_dx_neo2(3), dummy_sqrtg, dummy_iota

        x = (/0.0_dp, phi, theta/) !neo convention of x=(r, phi, theta)
        call neo_magfie_a(x, B_mod, dummy_sqrtg, dB_dx_neo2, dummy_iota)
        dB_dx(1) = dB_dx_neo2(1)
        dB_dx(2) = dB_dx_neo2(3) !neo convention of x=(r, phi, theta)
        dB_dx(3) = dB_dx_neo2(2) !neo convention of x=(r, phi, theta)

    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(neo_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: x(3), dummy_iota, dummy_sqrtg, dummy_dB_dx(3)

        x = (/0.0_dp, phi, theta/) !neo convention of x=(r, phi, theta)
        call neo_magfie_a(x, B_mod, dummy_sqrtg, dummy_dB_dx, dummy_iota)

    end subroutine compute_B_mod

    subroutine neo_change_stor(self, stor)
        use neo_magfie, only: magfie_newspline
        use neo_input, only: es

        class(neo_field_t), intent(inout) :: self
        real(dp), intent(in) :: stor

        real(dp) :: x(3), dummy_B_mod, dummy_sqrtg, dummy_dB_dx(3)

        if (magfie_newspline .ne. 1) magfie_newspline = 1
        x = (/stor, 0.0_dp, 0.0_dp/)
        call neo_magfie_a(x, dummy_B_mod, dummy_sqrtg, dummy_dB_dx, self%iota)
    end subroutine

end module neo_field
