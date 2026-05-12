program test_handle_multiple_maxima
    use constants, only: dp
    use make_fieldline, only: maxima_t, pick_maximum_on_each_side
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none

    integer, parameter :: n_maxima = 4

    call test_clear_winner()
    call test_tied_by_error_pick_closest()
    call test_tied_by_sym_tol_pick_closest()
    call test_no_false_tie()

contains

    subroutine assert_phi_max(label, phi_max, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: phi_max(2), expected(2)

        if (any(phi_max /= expected)) then
            print *, "-------------------------------------------------------------"
            print *, "FAILED: ", label
            print *, "phi_max(1): got ", phi_max(1), " expected ", expected(1)
            print *, "phi_max(2): got ", phi_max(2), " expected ", expected(2)
            error stop
        end if
    end subroutine assert_phi_max

    !> Initialize a maxima_t: unused phi entries set to NaN (as in production code),
    !> so they are excluded by both strict masks regardless of their value
    subroutine init_maxima(maxima)
        type(maxima_t), intent(out) :: maxima

        maxima%n = n_maxima
        maxima%phi = ieee_value(0.0_dp, ieee_quiet_nan)
        maxima%B = 0.0_dp
        maxima%error = 0.0_dp
    end subroutine init_maxima

    !> Two maxima on each side with clearly different B: higher B wins
    !> even when it is farther from phi_0
    subroutine test_clear_winner()
        type(maxima_t) :: maxima
        real(dp) :: phi_max(2)
        real(dp), parameter :: phi_0 = 0.0_dp, sym_tol = 0.0_dp

        call init_maxima(maxima)
        maxima%phi(1:n_maxima) = [-3.0_dp, -1.0_dp, 1.0_dp, 3.0_dp]
        maxima%B(1:n_maxima) = [2.0_dp, 1.0_dp, 1.0_dp, 2.5_dp]

        call pick_maximum_on_each_side(maxima, phi_0, phi_max, sym_tol)

        call assert_phi_max("clear_winner", phi_max, [-3.0_dp, 3.0_dp])
    end subroutine test_clear_winner

    !> Two maxima on each side whose B values differ by less than error_i + error_j:
    !> intervals overlap so the one closest to phi_0 wins
    subroutine test_tied_by_error_pick_closest()
        type(maxima_t) :: maxima
        real(dp) :: phi_max(2)
        real(dp), parameter :: phi_0 = 0.0_dp, sym_tol = 0.0_dp
        ! B difference = 0.01, combined error = 0.0 + 0.02 = 0.02 -> tied

        call init_maxima(maxima)
        maxima%phi(1:n_maxima) = [-3.0_dp, -1.0_dp, 1.0_dp, 3.0_dp]
        maxima%B(1:n_maxima) = [2.0_dp, 1.99_dp, 1.99_dp, 2.0_dp]
        maxima%error(2:3) = 0.02_dp

        call pick_maximum_on_each_side(maxima, phi_0, phi_max, sym_tol)

        call assert_phi_max("tied_by_error_pick_closest", phi_max, [-1.0_dp, 1.0_dp])
    end subroutine test_tied_by_error_pick_closest

    !> With zero errors, B difference within 2*symmetry_violation*B: tied,
    !> closest to phi_0 wins
    subroutine test_tied_by_sym_tol_pick_closest()
        type(maxima_t) :: maxima
        real(dp) :: phi_max(2)
        real(dp), parameter :: phi_0 = 0.0_dp
        ! B_max = 2.0, B difference = 0.01
        ! need 2*symmetry_violation*B_max > 0.01 -> symmetry_violation > 0.0025
        real(dp), parameter :: symmetry_violation = 0.003_dp

        call init_maxima(maxima)
        maxima%phi(1:n_maxima) = [-3.0_dp, -1.0_dp, 1.0_dp, 3.0_dp]
        maxima%B(1:n_maxima) = [2.0_dp, 1.99_dp, 1.99_dp, 2.0_dp]

        call pick_maximum_on_each_side(maxima, phi_0, phi_max, &
                                       symmetry_violation)

        call assert_phi_max("tied_by_sym_tol_pick_closest", phi_max, [-1.0_dp, 1.0_dp])
    end subroutine test_tied_by_sym_tol_pick_closest

    !> B difference just exceeds combined tolerance: no tie, higher B wins
    !> even when it is farther from phi_0
    subroutine test_no_false_tie()
        type(maxima_t) :: maxima
        real(dp) :: phi_max(2)
        real(dp), parameter :: phi_0 = 0.0_dp, sym_tol = 0.0_dp
        ! B difference = 0.02, combined error = 0.009 + 0.0 = 0.009 -> not tied

        call init_maxima(maxima)
        maxima%phi(1:n_maxima) = [-3.0_dp, -1.0_dp, 1.0_dp, 3.0_dp]
        maxima%B(1:n_maxima) = [2.0_dp, 1.98_dp, 1.98_dp, 2.0_dp]
        maxima%error(2:3) = 0.009_dp

        call pick_maximum_on_each_side(maxima, phi_0, phi_max, sym_tol)

        call assert_phi_max("no_false_tie", phi_max, [-3.0_dp, 3.0_dp])
    end subroutine test_no_false_tie

end program test_handle_multiple_maxima
