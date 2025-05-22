program test_real_dft
    use constants, only: dp, pi
    use utils, only: is_same, linspace
    use fourier, only: real_ft

    implicit none

    integer, parameter :: N = 8
    real(dp), parameter :: reltol = 1e-15

    real(dp), dimension(N) :: f, k, x
    real(dp), dimension(N/2 + 1) :: found_f_cos, found_f_sin, f_cos, f_sin
    complex(dp), dimension(N) :: f_exp

    real(dp) :: const
    integer :: k0, k1, k2
    integer :: n0
    integer :: j

    do j = 0, N - 1
        x(j + 1) = 2.0_dp*pi*real(j, kind=dp)/real(N, kind=dp)
    end do

    do k0 = 0, N/2
        f = cos(x*k0)
        call real_ft(x, f, found_f_cos, found_f_sin)
        f_cos(:) = 0.0_dp
        f_cos(k0 + 1) = 1.0_dp
        f_sin(:) = 0.0_dp

        if (is_same(f_cos, found_f_cos, reltol) .or. &
            is_same(f_sin, found_f_sin, reltol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_real_dft failed for case cos(2pi*k0*n/N), k0 = ", k0
            print *, "found f_cos: ", found_f_cos
            print *, "f_cos: ", f_cos
            print *, "found f_sin: ", found_f_sin
            print *, "f_sin: ", f_sin
            error stop
        end if
    end do

    call linspace(0.0_dp, real(N - 1, kind=dp), N, k)
    do n0 = 1, N - 1
        f(:) = 0.0_dp
        f(n0 + 1) = 1.0_dp
        call real_ft(x, f, found_f_cos, found_f_sin)
        f_exp = exp(-2.0_dp*pi*(0.0_dp, 1.0_dp)/real(N, kind=dp) &
                    *real(n0, kind=dp)*k)
        f_cos(:) = real(f_exp(1:N/2 + 1))/real(N, kind=dp)*2.0_dp
        f_sin(:) = -aimag(f_exp(1:N/2 + 1))/real(N, kind=dp)*2.0_dp
        f_cos(1) = f_cos(1)/2.0_dp
        if (mod(N, 2) == 0) f_cos(N/2 + 1) = f_cos(N/2 + 1)/2.0_dp

        if (is_same(f_cos, found_f_cos, reltol) .or. &
            is_same(f_sin, found_f_sin, reltol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_real_dft failed for case kronecker_(n,n0), n0 = ", n0
            print *, "found f_cos: ", found_f_cos
            print *, "f_cos: ", f_cos
            print *, "found f_sin: ", found_f_sin
            print *, "f_sin: ", f_sin
            error stop
        end if
    end do

    const = 2.0_dp
    f(:) = const
    call real_ft(x, f, found_f_cos, found_f_sin)
    f_cos(:) = 0.0_dp
    f_sin(:) = 0.0_dp
    f_cos(1) = const

    if (is_same(f_cos, found_f_cos, reltol) .or. &
        is_same(f_sin, found_f_sin, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_real_dft failed for case const = ", const
        print *, "found f_cos: ", found_f_cos
        print *, "f_cos: ", f_cos
        print *, "found f_sin: ", found_f_sin
        print *, "f_sin: ", f_sin
        error stop
    end if

    k1 = 2
    k2 = 4
    f = 2.0_dp*sin(x*k1) - 3.0_dp*cos(x*k2)
    call real_ft(x, f, found_f_cos, found_f_sin)
    f_cos(:) = 0.0_dp
    f_cos(k2 + 1) = -3.0_dp
    f_sin(:) = 0.0_dp
    f_sin(k1 + 1) = 2.0_dp

    if (is_same(f_cos, found_f_cos, reltol) .or. &
        is_same(f_sin, found_f_sin, reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_real_dft failed for case sum of cos and sin"
        print *, "found f_cos: ", found_f_cos
        print *, "f_cos: ", f_cos
        print *, "found f_sin: ", found_f_sin
        print *, "f_sin: ", f_sin
        error stop
    end if

end program test_real_dft
