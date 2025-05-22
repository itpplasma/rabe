module fourier
    use constants, only: dp, pi
    implicit none

contains

    subroutine real_dft(f, f_cos, f_sin)
        real(dp), dimension(:), intent(in) :: f
        real(dp), dimension(size(f)/2 + 1), intent(out) :: f_cos, f_sin

        integer :: N, k, j
        complex(dp) :: phase_step, phase, sum
        complex(dp), dimension(size(f)/2 + 1) :: f_exp
        complex(dp), parameter :: Icmplx = (0.0_dp, 1.0_dp)

        N = size(f)

        do k = 0, N/2
            phase_step = exp(-2.0d0*pi*Icmplx*real(k, kind=dp)/real(N, kind=dp))
            phase = (1.0_dp, 0.0_dp)
            sum = (0.0_dp, 0.0_dp)
            do j = 0, N - 1
                sum = sum + f(j + 1)*phase
                phase = phase*phase_step
            end do
            f_exp(k + 1) = sum
        end do

        f_cos(1) = real(f_exp(1))
        f_sin(1) = 0.0_dp
        do k = 1, N/2
            f_cos(k + 1) = 2.0_dp*real(f_exp(k + 1))
            f_sin(k + 1) = -2.0_dp*aimag(f_exp(k + 1))
        end do

        if (mod(N, 2) == 0) then
            f_cos(N/2 + 1) = f_cos(N/2 + 1)/2.0_dp
            f_sin(N/2 + 1) = 0.0_dp
        end if

        f_cos = f_cos/real(N, kind=dp)
        f_sin = f_sin/real(N, kind=dp)

    end subroutine real_dft

end module fourier
