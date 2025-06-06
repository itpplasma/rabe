program test_against_neo2_anti_sigma
    use neo_field, only: neo_field_t
    use constants, only: dp, pi

    implicit none

    real(dp) :: retol = 1e-11
    character(len=*), parameter :: bc_filename = "input/helical_anti.bc"

    type(neo_field_t) :: field
    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: bmod
    real(dp) :: bmod_analytic(n_cases), bmod_neo2(n_cases)
    real(dp), parameter :: B_min = 0.875_dp, B_max = 1.05_dp
    integer :: case

    logical :: test_failed

    test_failed = .false.

    stor = (/0.50_dp, 0.70_dp, 0.80_dp, 0.9999_dp/)
    theta = (/0.0_dp, 0.25_dp*pi, 2.0_dp*pi, 0.0_dp/)
    phi = (/0.0_dp, 0.15_dp*pi, 0.0_dp, -0.1_dp*pi/)

    bmod_analytic = (/B_min, &
                      B_max, &
                      B_min, &
                      B_max/)

    bmod_neo2 = (/0.86660052969026313_dp, &
                  1.0488269770770471_dp, &
                  0.87920244657041735_dp, &
                  1.0499999999555609_dp/)

    call field%neo_field_init(bc_filename, stor(1))
    do case = 1, n_cases
        call field%neo_change_stor(stor(case))
        call field%compute_B_mod(theta(case), phi(case), bmod)
        if (abs(bmod/bmod_analytic(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_neo2_field failed: B analytic"
            print *, "B: ", bmod
            print *, "analytic: ", bmod_analytic(case)
            print *, "relative error: ", 1.0_dp - bmod/bmod_analytic(case)
            test_failed = .true.
        end if
        if (abs(bmod/bmod_neo2(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_neo2_field failed: B NEO-2"
            print *, "B: ", bmod
            print *, "NEO-2: ", bmod_neo2(case)
            print *, "relative error: ", 1.0_dp - bmod/bmod_neo2(case)
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_neo2_anti_sigma
