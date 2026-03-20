program test_linspace
    implicit none

    logical :: test_failed

    test_failed = .false.

    call test_linspace_with_endpoint(test_failed)
    call test_linspace_without_endpoint(test_failed)
    call test_linspace_default_endpoint(test_failed)

    if (test_failed) error stop

contains

    subroutine test_linspace_with_endpoint(test_failed)
        use constants, only: dp
        use utils, only: linspace, not_same

        logical, intent(inout) :: test_failed
        integer, parameter :: n = 5
        real(dp) :: x(n)
        real(dp), parameter :: expected(5) = &
                               [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]

        call linspace(0.0_dp, 1.0_dp, n, x, include_endpoint=.true.)

        if (not_same(x, expected, abstol_in=1.0e-15_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_linspace failed: with endpoint"
            print *, "found: ", x
            print *, "expected: ", expected
            test_failed = .true.
        end if
    end subroutine test_linspace_with_endpoint

    subroutine test_linspace_without_endpoint(test_failed)
        use constants, only: dp
        use utils, only: linspace, not_same

        logical, intent(inout) :: test_failed
        integer, parameter :: n = 4
        real(dp) :: x(n)
        real(dp), parameter :: expected(4) = &
                               [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp]

        call linspace(0.0_dp, 1.0_dp, n, x, include_endpoint=.false.)

        if (not_same(x, expected, abstol_in=1.0e-15_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_linspace failed: without endpoint"
            print *, "found: ", x
            print *, "expected: ", expected
            test_failed = .true.
        end if
    end subroutine test_linspace_without_endpoint

    subroutine test_linspace_default_endpoint(test_failed)
        use constants, only: dp
        use utils, only: linspace, not_same

        logical, intent(inout) :: test_failed
        integer, parameter :: n = 3
        real(dp) :: x(n)
        real(dp), parameter :: expected(3) = [1.0_dp, 2.0_dp, 3.0_dp]

        call linspace(1.0_dp, 3.0_dp, n, x)

        if (not_same(x, expected, abstol_in=1.0e-15_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_linspace failed: default (should include endpoint)"
            print *, "found: ", x
            print *, "expected: ", expected
            test_failed = .true.
        end if
    end subroutine test_linspace_default_endpoint

end program test_linspace
