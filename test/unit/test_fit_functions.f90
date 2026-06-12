program test_fit_functions
    use constants, only: dp, pi
    use fit_functions, only: S_A, S_B
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-15_dp

    real(dp), parameter :: positive_angle = 0.25_dp*pi
    real(dp), parameter :: negative_angle = -0.25_dp*pi
    real(dp), parameter :: positive_out_of_period_angle = 1.3_dp*pi
    real(dp), parameter :: negative_out_of_period_angle = -1.25_dp*pi

    real(dp), parameter :: S_A_positive_angle = 0.26_dp*(-0.25_dp*pi)
    real(dp), parameter :: S_A_negative_angle = 0.26_dp*(+0.25_dp*pi)
    real(dp), parameter :: S_A_positive_out_of_period_angle = 0.26_dp*(-0.20_dp*pi)
    real(dp), parameter :: S_A_negative_out_of_period_angle = 0.26_dp*(+0.25_dp*pi)

    real(dp), parameter :: S_B_positive_angle = 1.85_dp
    real(dp), parameter :: S_B_negative_angle = -1.85_dp
    real(dp), parameter :: S_B_positive_out_of_period_angle = -1.85_dp
    real(dp), parameter :: S_B_negative_out_of_period_angle = 1.85_dp

    call check_S_A(positive_angle, S_A(positive_angle), S_A_positive_angle)
    call check_S_A(negative_angle, S_A(negative_angle), S_A_negative_angle)
    call check_S_A(negative_out_of_period_angle, &
                   S_A(negative_out_of_period_angle), &
                   S_A_negative_out_of_period_angle)
    call check_S_A(positive_out_of_period_angle, &
                   S_A(positive_out_of_period_angle), &
                   S_A_positive_out_of_period_angle)

    call check_S_B(positive_angle, S_B(positive_angle), S_B_positive_angle)
    call check_S_B(negative_angle, S_B(negative_angle), S_B_negative_angle)
    call check_S_B(negative_out_of_period_angle, &
                   S_B(negative_out_of_period_angle), &
                   S_B_negative_out_of_period_angle)
    call check_S_B(positive_out_of_period_angle, &
                   S_B(positive_out_of_period_angle), &
                   S_B_positive_out_of_period_angle)

contains

    subroutine check_S_A(x, value, expected_value)
        real(dp) :: x, value, expected_value

        if (not_same(expected_value, value, reltol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_fit_functions failed: for S_A at angle ", x
            print *, "found: ", value
            print *, "expected: ", expected_value
            error stop
        end if
    end subroutine check_S_A

    subroutine check_S_B(x, value, expected_value)
        real(dp) :: x, value, expected_value

        if (not_same(expected_value, value, reltol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_fit_functions failed: for S_B at angle ", x
            print *, "found: ", value
            print *, "expected: ", expected_value
            error stop
        end if
    end subroutine check_S_B

end program test_fit_functions
