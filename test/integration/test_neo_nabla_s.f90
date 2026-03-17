program test_neo_nabla_s
    use neo_field, only: neo_field_t
    use constants, only: dp, pi
    use utils, only: not_same

    implicit none

    real(dp) :: reltol = 1e-14
    character(len=*), parameter :: bc_filename = "input/quasi_helical.bc"

    type(neo_field_t) :: field
    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: nabla_psi(n_cases)
    real(dp) :: found

    integer :: case
    logical :: test_failed

    test_failed = .false.

    stor = [0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp]
    theta = [-1.00_dp, 4.1_dp*pi, 0.25_dp*pi, -1.0_dp*pi]
    phi = [1.00_dp, -10.1_dp*pi, 0.65_dp*pi, 2.0_dp*pi]

    nabla_psi = [2.0053230989059028_dp, &
                 7.4603992604910152_dp, &
                 11.978734591413621_dp, &
                 19.901348996309611_dp]

    call field%neo_field_init(bc_filename, stor(1))
    do case = 1, n_cases
        call field%compute_nabla_s(theta(case), phi(case), found)
        found = found*field%psi_tor_edge
        if (not_same(nabla_psi(case), &
                     found, &
                     reltol_in=reltol, &
                     abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_neo_nabla_s failed: nabla psi for case ", case
            print *, "expected: ", nabla_psi(case)
            print *, "found: ", found
            print *, "relative difference: ", &
                abs(found/nabla_psi(case) - 1.0_dp)
            test_failed = .true.
        end if
        if (case + 1 .le. n_cases) call field%neo_change_stor(stor(case + 1))
    end do

    if (test_failed) error stop

end program test_neo_nabla_s
