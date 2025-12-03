module test_cumint_mod
    use constants, only: dp
    implicit none

    abstract interface
        function func(x)
            import :: dp
            real(dp), intent(in) :: x
            real(dp) :: func
        end function func
    end interface

    procedure(func), pointer, private :: trial_func_to_test => null()

contains

    function cumint_not_equal_antideriv(trial_func, interval, antiderivative_func)
        use utils, only: not_same, linspace
        use shaing_callen_integration, only: cumint
        procedure(func) :: trial_func, antiderivative_func
        real(dp), intent(in) :: interval(2)
        logical :: cumint_not_equal_antideriv

        integer, parameter :: n_x = 10
        real(dp), dimension(:), allocatable :: x
        real(dp), dimension(:), allocatable :: antiderivative, found_antiderivative
        real(dp), parameter :: reltol = 1e-8, abstol = 1e-15
        integer :: current

        cumint_not_equal_antideriv = .false.

        allocate (x(n_x), found_antiderivative(n_x), antiderivative(n_x))
        call linspace(interval(1), interval(2), n_x, x)
        trial_func_to_test => trial_func
        found_antiderivative = cumint(x, defined_integral)
        do current = 1, size(x)
            antiderivative(current) = antiderivative_func(x(current))
        end do
        antiderivative = antiderivative - antiderivative(1)

        if (not_same(antiderivative, found_antiderivative, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_cumint failed:"
            print *, "found: ", found_antiderivative
            print *, "analytic: ", antiderivative
            print *, "maximum error: ", maxval(abs(found_antiderivative - &
                                                   antiderivative))
            cumint_not_equal_antideriv = .true.
        end if

        deallocate (x, found_antiderivative, antiderivative)

    end function cumint_not_equal_antideriv

    function defined_integral(x_start, x_end)
        use integrate, only: integrate_1d_substituted
        real(dp), intent(in) :: x_start, x_end
        real(dp) :: defined_integral

        call integrate_1d_substituted(trial_func_to_test, &
                                      x_start, &
                                      x_end, &
                                      defined_integral)
    end function defined_integral

    function trial_func_0(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_0

        trial_func_0 = 3.0_dp*x**2
    end function trial_func_0

    function antiderivative_0(x)
        real(dp), intent(in) :: x
        real(dp) :: antiderivative_0

        antiderivative_0 = x**3
    end function antiderivative_0

    function trial_func_1(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_1

        trial_func_1 = sin(x)
    end function trial_func_1

    function antiderivative_1(x)
        real(dp), intent(in) :: x
        real(dp) :: antiderivative_1

        antiderivative_1 = -cos(x)
    end function antiderivative_1

    function trial_func_2(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_2

        trial_func_2 = 0.5_dp/sqrt(1.0_dp + x)
    end function trial_func_2

    function antiderivative_2(x)
        real(dp), intent(in) :: x
        real(dp) :: antiderivative_2

        antiderivative_2 = sqrt(1.0_dp + x)
    end function antiderivative_2

    function trial_func_3(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_3

        trial_func_3 = exp(x)
    end function trial_func_3

    function antiderivative_3(x)
        real(dp), intent(in) :: x
        real(dp) :: antiderivative_3

        antiderivative_3 = exp(x) - 1.0_dp
    end function antiderivative_3

end module test_cumint_mod
