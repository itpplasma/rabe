program test_error_handling
    use error_handling, only: set_unsafe_mode, failed_sanity_check, &
                              reset_failed_check_counter, did_fail_any_sanity_check
    implicit none

    logical :: test_failed

    test_failed = .false.

    call test_no_failure_after_reset(test_failed)
    call test_failure_detected(test_failed)
    call test_multiple_failures_detected(test_failed)
    call test_reset_clears_counter(test_failed)

    if (test_failed) error stop

contains

    subroutine test_no_failure_after_reset(test_failed)
        logical, intent(inout) :: test_failed

        call reset_failed_check_counter()

        if (did_fail_any_sanity_check()) then
            print *, "-------------------------------------------------------------"
            print *, "test_error_handling failed: did_fail_any_sanity_check()" &
                //" should be .false. after reset"
            test_failed = .true.
        end if
    end subroutine test_no_failure_after_reset

    subroutine test_failure_detected(test_failed)
        logical, intent(inout) :: test_failed

        call reset_failed_check_counter()
        call set_unsafe_mode(.true.)
        call failed_sanity_check()
        call set_unsafe_mode(.false.)

        if (.not. did_fail_any_sanity_check()) then
            print *, "-------------------------------------------------------------"
            print *, "test_error_handling failed: did_fail_any_sanity_check()" &
                //" should be .true. after failed_sanity_check() in unsafe mode"
            test_failed = .true.
        end if
    end subroutine test_failure_detected

    subroutine test_multiple_failures_detected(test_failed)
        logical, intent(inout) :: test_failed

        call reset_failed_check_counter()
        call set_unsafe_mode(.true.)
        call failed_sanity_check()
        call failed_sanity_check()
        call failed_sanity_check()
        call set_unsafe_mode(.false.)

        if (.not. did_fail_any_sanity_check()) then
            print *, "-------------------------------------------------------------"
            print *, "test_error_handling failed: did_fail_any_sanity_check()" &
                //" should be .true. after multiple failed_sanity_check() calls"
            test_failed = .true.
        end if
    end subroutine test_multiple_failures_detected

    subroutine test_reset_clears_counter(test_failed)
        logical, intent(inout) :: test_failed

        call set_unsafe_mode(.true.)
        call failed_sanity_check()
        call set_unsafe_mode(.false.)
        call reset_failed_check_counter()

        if (did_fail_any_sanity_check()) then
            print *, "-------------------------------------------------------------"
            print *, "test_error_handling failed: did_fail_any_sanity_check()" &
                //" should be .false. after reset_failed_check_counter()"
            test_failed = .true.
        end if
    end subroutine test_reset_clears_counter

end program test_error_handling
