program test_misc
    use constants, only: dp, pi
    use misc, only: S_A
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-15_dp
    integer, parameter :: n_points = 10
    real(dp), dimension(n_points) :: angles
    real(dp), parameter :: positive_angle = 0.25_dp*pi
    real(dp), parameter :: negative_angle = -0.25_dp*pi
    real(dp), parameter :: negative_out_of_period_angle = -1.25_dp*pi
    real(dp), parameter :: positive_out_of_period_angle = 1.3_dp*pi

    real(dp) :: S_A_positive_angle = 0.26_dp*(-0.25_dp*pi)
    real(dp) :: S_A_negative_angle = 0.26_dp*(+0.25_dp*pi)
    real(dp) :: S_A_positive_out_of_period_angle = 0.26_dp*(-0.20_dp*pi)

    if (not_same(S_A_positive_angle, &
                 S_A(positive_angle), &
                 reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_misc failed: for angle ", positive_angle
        print *, "found: ", S_A(positive_angle)
        print *, "expected: ", S_A_positive_angle
        error stop
    end if

    if (not_same(S_A_negative_angle, &
                 S_A(negative_angle), &
                 reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_misc failed: for angle ", negative_angle
        print *, "found: ", S_A(negative_angle)
        print *, "expected: ", S_A_negative_angle
        error stop
    end if

    if (not_same(S_A_negative_angle, &
                 S_A(negative_out_of_period_angle), &
                 reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_misc failed: for angle ", negative_out_of_period_angle
        print *, "found: ", S_A(negative_out_of_period_angle)
        print *, "expected: ", S_A_negative_angle
        error stop
    end if

    if (not_same(S_A_positive_out_of_period_angle, &
                 S_A(positive_out_of_period_angle), &
                 reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_misc failed: for angle ", positive_out_of_period_angle
        print *, "found: ", S_A(positive_out_of_period_angle)
        print *, "expected: ", S_A_positive_out_of_period_angle
        error stop
    end if

end program test_misc
