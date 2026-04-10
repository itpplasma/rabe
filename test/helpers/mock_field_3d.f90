module mock_field_3d
    use constants, only: dp, pi
    use field_base, only: field_t
    use field_base, only: field_3D_t

    implicit none

    type, extends(field_3D_t) :: mock_field_3d_t
        class(field_t), allocatable :: field_2D
    contains
        procedure :: mock_field_3d_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: evaluate
        procedure :: get_iota
    end type mock_field_3d_t

contains

    subroutine mock_field_3d_init(self, field_2D_in)
        class(mock_field_3d_t), intent(out) :: self
        class(field_t), intent(in) :: field_2D_in

        allocate (self%field_2D, source=field_2D_in)

    end subroutine mock_field_3d_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        call self%field_2D%compute_B_sqrtg_dB_dx(theta, phi, B_mod, sqrtg, dB_dx)

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        call self%field_2D%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)

    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        call self%field_2D%compute_B_mod(theta, phi, B_mod)

    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        call self%field_2D%compute_nabla_s(theta, phi, nabla_s)

    end subroutine compute_nabla_s

    subroutine evaluate(self, x, bmod, sqrtg, bder, &
                        hcovar, hctrvr, hcurl)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        use constants, only: machine_eps
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp) :: B_mod, dB_dx(3)

        if (any(ieee_is_nan(x))) then
            print *, "Input x is NaN!"
            print *, "x = ", x
            error stop
        end if

        call self%field_2D%compute_B_and_dB_dx(x(2), x(3), B_mod, dB_dx)

        sqrtg = 1.0_dp
        bmod = B_mod
        bder = dB_dx/B_mod
        hcovar = 0.0_dp
        hcovar(3) = 1.0_dp
        hctrvr = 0.0_dp
        hctrvr(2) = 1.0_dp
        hctrvr(3) = 1.0_dp
        hcurl = 0.0_dp

        if (ieee_is_nan(bmod)) then
            print *, "bmod is NaN!"
            error stop
        end if
        if (ieee_is_nan(sqrtg)) then
            print *, "sqrtg is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(bder))) then
            print *, "bder is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(hcovar))) then
            print *, "hcovar is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(hctrvr))) then
            print *, "hctrvr is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(hcurl))) then
            print *, "hcurl is NaN!"
            error stop
        end if
        if (abs(bmod) < machine_eps) then
            print *, "bmod is zero!"
            error stop
        end if
        if (abs(sqrtg) < machine_eps) then
            print *, "sqrtg is zero!"
            error stop
        end if
    end subroutine evaluate

    subroutine get_iota(self, s_tor, iota)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: s_tor
        real(dp), intent(out) :: iota

        real(dp) :: dummy(11), hctrvr(3)
        real(dp) :: x(3)

        x(1) = s_tor
        x(2) = 0.0_dp
        x(3) = 0.0_dp

        call self%evaluate(x, dummy(1), dummy(2), dummy(3:5), &
                           dummy(6:8), hctrvr, dummy(9:11))

        iota = hctrvr(2)/hctrvr(3)

    end subroutine get_iota

end module mock_field_3d
