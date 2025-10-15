program test_cumint
    use constants, only: dp, pi
    use integrate, only: sum_trapez_1d
    use test_cumint_mod, only: trial_func_0, antiderivative_0
    use test_cumint_mod, only: trial_func_1, antiderivative_1
    use test_cumint_mod, only: trial_func_2, antiderivative_2
    use test_cumint_mod, only: trial_func_3, antiderivative_3
    use test_cumint_mod, only: cumint_not_equal_antideriv

    implicit none

    real(dp), dimension(2), parameter :: interval_0 = [-1.0_dp, 1.0_dp]
    real(dp), dimension(2), parameter :: interval_1 = [-pi, pi]
    real(dp), dimension(2), parameter :: interval_2 = [-1.0_dp, 1.0_dp]
    real(dp), dimension(2), parameter :: interval_3 = [-1.0_dp, 1.0_dp]

    logical :: failed_test

    failed_test = .false.

    print *, "trial function 0:"
    if (cumint_not_equal_antideriv(trial_func_0, interval_0, antiderivative_0)) then
        failed_test = .true.
    end if
    print *, "trial function 1:"
    if (cumint_not_equal_antideriv(trial_func_1, interval_1, antiderivative_1)) then
        failed_test = .true.
    end if
    print *, "trial function 2:"
    if (cumint_not_equal_antideriv(trial_func_2, interval_2, antiderivative_2)) then
        failed_test = .true.
    end if
    print *, "trial function 3:"
    if (cumint_not_equal_antideriv(trial_func_3, interval_3, antiderivative_3)) then
        failed_test = .true.
    end if

    if (failed_test) error stop

contains

end program test_cumint
