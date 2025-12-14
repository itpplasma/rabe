program test_neo_sqrt_g11
    use neo_field, only: neo_field_t
    use constants, only: dp, pi
    use utils, only: not_same

    implicit none

    real(dp) :: reltol = 1e-14
    character(len=*), parameter :: bc_filename = "input/quasi_helical.bc"

    type(neo_field_t) :: field
    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: sqrt_g11(n_cases)
    real(dp) :: found

    integer :: case
    logical :: test_failed

    test_failed = .false.

    stor = [0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp]
    theta = [-1.00_dp, 4.1_dp*pi, 0.25_dp*pi, -1.0_dp*pi]
    phi = [1.00_dp, -10.1_dp*pi, 0.65_dp*pi, 2.0_dp*pi]

    sqrt_g11 = [2.0053230989059028_dp, &
                7.4603992604910152_dp, &
                11.978734591413621_dp, &
                19.901348996309611_dp]

    call field%neo_field_init(bc_filename, stor(1))
    do case = 1, n_cases
        call field%compute_sqrt_g11(theta(case), phi(case), found)
        if (not_same(sqrt_g11(case), &
                     found, &
                     reltol_in=reltol, &
                     abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_neo_sqrt_g11 failed: sqrt(g11) for case ", case
            print *, "expected: ", sqrt_g11(case)
            print *, "found: ", found
            print *, "relative difference: ", &
                abs(found/sqrt_g11(case) - 1.0_dp)
            test_failed = .true.
        end if
        if (case + 1 .le. n_cases) call field%neo_change_stor(stor(case + 1))
    end do

    if (test_failed) error stop

end program test_neo_sqrt_g11
