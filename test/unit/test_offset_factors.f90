program test_offset_factors
    use constants, only: dp, pi
    use offset_factors, only: S_A, S_B, set_offset_factor_mode, &
                              base_wiener_hopf_residual
    use utils, only: not_same

    implicit none

    !> Discretization tolerance: the exact factors are read off a finite grid,
    !! so match the Python oracle table (points=121, extent=7.0) to ~1e-6.
    real(dp), parameter :: factor_abstol = 1.0e-6_dp
    real(dp), parameter :: zero_abstol = 1.0e-9_dp

    real(dp), parameter :: angles(5) = [0.0_dp, pi/6.0_dp, pi/4.0_dp, &
                                        pi/3.0_dp, pi/2.0_dp]
    real(dp), parameter :: expected_S_A(5) = [ &
                           0.0_dp, &
                           -1.5827902093804535e-02_dp, &
                           -1.1862821398844337e-02_dp, &
                           -7.8745528453801287e-03_dp, &
                           0.0_dp]
    real(dp), parameter :: expected_S_B(5) = [ &
                           0.0_dp, &
                           1.2766385382097867e-01_dp, &
                           1.3329636074423146e-01_dp, &
                           1.3683900710093186e-01_dp, &
                           1.3948146565032457e-01_dp]

    real(dp) :: got_S_A(5), got_S_B(5)
    real(dp) :: residual
    integer :: i

    call set_offset_factor_mode(.true.)

    residual = base_wiener_hopf_residual()
    if (residual > 1.0e-10_dp) then
        print *, "test_offset_factors failed: base Wiener-Hopf residual ", residual
        error stop
    end if

    got_S_A = S_A(angles)
    got_S_B = S_B(angles)

    do i = 1, size(angles)
        call check("S_A", angles(i), got_S_A(i), expected_S_A(i))
        call check("S_B", angles(i), got_S_B(i), expected_S_B(i))
    end do

    ! Invariants documented by the oracle.
    if (abs(got_S_A(1)) > zero_abstol .or. abs(got_S_B(1)) > zero_abstol) then
        print *, "test_offset_factors failed: zero-angle symmetry"
        print *, "S_A(0) = ", got_S_A(1), " S_B(0) = ", got_S_B(1)
        error stop
    end if
    if (abs(got_S_A(5)) > zero_abstol) then
        print *, "test_offset_factors failed: S_A(pi/2) not zero: ", got_S_A(5)
        error stop
    end if

    ! mode=fit must restore the analytic fit factors unchanged.
    call set_offset_factor_mode(.false.)
    call check_fit_default()

    print *, "test_offset_factors passed: residual ", residual

contains

    subroutine check(name, angle, value, expected_value)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: angle, value, expected_value

        if (not_same(expected_value, value, reltol_in=0.0_dp, &
                     abstol_in=factor_abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_offset_factors failed: for ", name, " at angle ", angle
            print *, "found: ", value
            print *, "expected: ", expected_value
            error stop
        end if
    end subroutine check

    subroutine check_fit_default()
        use fit_functions, only: S_A_ref => S_A, S_B_ref => S_B

        real(dp) :: probe(3)
        real(dp), parameter :: tight = 1.0e-15_dp

        probe = [0.25_dp*pi, -0.25_dp*pi, 1.3_dp*pi]
        if (not_same(S_A(probe), S_A_ref(probe), tight)) then
            print *, "test_offset_factors failed: fit S_A changed under mode=fit"
            error stop
        end if
        if (not_same(S_B(probe), S_B_ref(probe), tight)) then
            print *, "test_offset_factors failed: fit S_B changed under mode=fit"
            error stop
        end if
    end subroutine check_fit_default

end program test_offset_factors
