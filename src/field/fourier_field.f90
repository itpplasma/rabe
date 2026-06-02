module fourier_field

    use constants, only: dp, pi
    use field_base, only: field_t
    use interpolate, only: SplineData2D, construct_splines_2d, evaluate_splines_2d_der

    implicit none

    integer, parameter :: n_grid_default = 200

    type, extends(field_t) :: fourier_field_t
        type(SplineData2D) :: spl
        real(dp) :: nfp
        integer :: n_grid, mn_max
    contains
        procedure :: fourier_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
    end type fourier_field_t

contains

    ! Build 2D spline from flat Fourier mode lists.
    ! B = sum_k B_mn(k) * cos(m(k)*theta - n(k)*phi)
    subroutine fourier_field_init(self, m, n, B_mn, nfp, n_grid)
        class(fourier_field_t), intent(out) :: self
        integer, intent(in) :: m(:), n(:)
        real(dp), intent(in) :: B_mn(:)
        real(dp), intent(in), optional :: nfp
        integer, intent(in), optional :: n_grid

        integer :: n_theta, n_phi, i_theta, i_phi
        real(dp) :: theta, phi, dtheta, dphi
        real(dp), allocatable :: grid_B(:, :)

        if (present(nfp)) then
            self%nfp = nfp
        else
            self%nfp = 1.0_dp
        end if

        n_theta = n_grid_default
        n_phi = n_grid_default
        if (present(n_grid)) then
            n_theta = n_grid
            n_phi = n_grid
        end if
        self%n_grid = n_theta
        self%mn_max = max(maxval(abs(m)), maxval(abs(n)))

        dtheta = 2.0_dp*pi/real(n_theta - 1, dp)
        dphi = 2.0_dp*pi/(self%nfp*real(n_phi - 1, dp))

        allocate (grid_B(n_theta, n_phi))

        do i_phi = 1, n_phi
            phi = real(i_phi - 1, dp)*dphi
            do i_theta = 1, n_theta
                theta = real(i_theta - 1, dp)*dtheta
                grid_B(i_theta, i_phi) = sum(B_mn* &
                                      cos(real(m, dp)*theta - real(n, dp)*self%nfp*phi))
            end do
        end do

        call construct_splines_2d( &
            x_min=[0.0_dp, 0.0_dp], &
            x_max=[2.0_dp*pi, 2.0_dp*pi/self%nfp], &
            y=grid_B, &
            order=[5, 5], &
            periodic=[.true., .true.], &
            spl=self%spl)

    end subroutine fourier_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "fourier_field_t does not provide sqrtg. "// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: x(2), dy(2)

        x = [theta, phi]
        call evaluate_splines_2d_der(self%spl, x, B_mod, dy)
        dB_dx(1) = 0.0_dp
        dB_dx(2) = dy(1)
        dB_dx(3) = dy(2)

    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: dummy(3)

        call self%compute_B_and_dB_dx(theta, phi, B_mod, dummy)

    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        nabla_s = 1.0_dp

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(fourier_field_t), intent(in) :: self

        rel_accuracy_B = (real(self%mn_max, dp)/real(self%n_grid, dp))**6
    end function rel_accuracy_B

end module fourier_field
