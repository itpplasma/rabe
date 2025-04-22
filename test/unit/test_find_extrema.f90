program test_find_extrema
    use constants, only: dp, pi

    implicit none

    call test_find_local_minima()
    call test_find_local_maxima()
    call test_find_global_extrema()

contains

    subroutine test_find_local_minima()
        use find_extrema, only: find_local_minima

        real(dp), parameter :: tol = 1e-2
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 2.0_dp*pi/)
        real(dp), dimension(2), parameter :: expect_locs = (/2.1386_dp, 5.64832_dp/)
        real(dp) :: found_loc(1), found_locs(2)

        call find_local_minima(negative_sincos_func, interval, found_loc)
        call find_local_minima(negative_sincos_func, interval, found_locs)

        if (any(abs(found_loc/expect_locs(1) - 1.0_dp) > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_minima"
            print *, "found 1st location: ", found_loc
            print *, "expected 1st location: ", expect_locs(1)
            error stop
        end if

        if (any(abs(found_locs/expect_locs - 1.0_dp) > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_minima"
            print *, "found locations: ", found_locs
            print *, "expected locations: ", expect_locs
            error stop
        end if
    end subroutine test_find_local_minima

    subroutine negative_sincos_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        call sincos_func(x, value)
        value = -value
    end subroutine negative_sincos_func

    subroutine test_find_local_maxima()
        use find_extrema, only: find_local_maxima

        real(dp), parameter :: tol = 1e-2
        real(dp), dimension(2), parameter :: interval = (/0.0_dp, 2.0_dp*pi/)
        real(dp), dimension(2), parameter :: expect_locs = (/2.1386_dp, 5.64832_dp/)
        real(dp) :: found_loc(1), found_locs(2)

        call find_local_maxima(sincos_func, interval, found_loc)
        call find_local_maxima(sincos_func, interval, found_locs)

        if (any(abs(found_loc/expect_locs(1) - 1.0_dp) > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_maxima"
            print *, "found 1st location: ", found_loc
            print *, "expected 1st location: ", expect_locs(1)
            error stop
        end if

        if (any(abs(found_locs/expect_locs - 1.0_dp) > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_maxima"
            print *, "found locations: ", found_locs
            print *, "expected locations: ", expect_locs
            error stop
        end if

    end subroutine test_find_local_maxima

    subroutine sincos_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        value = cos(x) - sin(2*x)
    end subroutine sincos_func

    subroutine test_find_global_extrema
        use find_extrema, only: find_global_extrema

        real(dp) :: extrema(2), expected_extrema(2)
        real(dp), parameter :: interval(2) = [0.0_dp, 2.0_dp*pi]
        real(dp), parameter :: tol = 1e-3

        expected_extrema = (/-1.0_dp, 1.0_dp/)
        extrema = find_global_extrema(sin_func, interval)

        if (any(abs(extrema/expected_extrema - 1) > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_global_extrema"
            error stop
        end if
    end subroutine test_find_global_extrema

    subroutine sin_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        value = sin(x)
    end subroutine sin_func

end program test_find_extrema
