program test_against_neo2_field
    use neo_field, only: neo_field_t
    use constants, only: dp

    implicit none

    real(dp) :: retol = 1e-11
    character(len=*), parameter :: bc_filename = "input/quasi_helical.bc"

    type(neo_field_t) :: field
    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: bmod, sqrtg, dB_dx(3)
    real(dp) :: bmod_neo2(n_cases), sqrtg_neo2(n_cases), dB_dx_neo2(n_cases, 3)
    integer :: case

    logical :: test_failed

    test_failed = .false.

    stor = [0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp]
    theta = [-1.00_dp, 0.00_dp, 4.3498975203550980_dp, 1.00_dp]
    phi = [1.00_dp, -1.00_dp, 0.56947933720590893_dp, 0.00_dp]

    bmod_neo2 = [5.8461732541782538_dp, &
                 6.3747018649171556_dp, &
                 6.3402245795215624_dp, &
                 5.4298976321806043_dp]
    sqrtg_neo2 = [-19930779.196453653_dp, &
                  -16763094.050708901_dp, &
                  -16945898.197301224_dp, &
                  -23104179.862998683_dp]
    dB_dx_neo2(1,:) = [-0.92227365739728728_dp,       0.13109042596457662_dp,      -0.52436170385830649_dp]
    dB_dx_neo2(2,:) = [0.53462460260148859_dp,      -0.55835449392655978_dp,        2.2334179757062391_dp]
  dB_dx_neo2(3, :) = [0.340429548051560_dp, 0.784044027926972_dp, -3.136176111707888_dp]
    dB_dx_neo2(4,:) = [-0.20339965461773155_dp,       0.77633188889709603_dp,       -3.1053275555883841_dp]

    call field%neo_field_init(bc_filename, stor(1))
    do case = 1, n_cases
        call field%neo_change_stor(stor(case))
        call field%compute_B_sqrtg_dB_dx(theta(case), phi(case), bmod, sqrtg, dB_dx)
        if (abs(bmod/bmod_neo2(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_neo2_field failed: B"
            print *, "B: ", bmod
            print *, "NEO-2: ", bmod_neo2(case)
            test_failed = .true.
        end if
        if (abs(sqrtg/sqrtg_neo2(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_neo2_field failed: sqrtg"
            print *, "sqrtg: ", sqrtg
            print *, "NEO-2: ", sqrtg_neo2(case)
            test_failed = .true.
        end if
        if (any(abs(dB_dx/dB_dx_neo2(case, :) - 1) > retol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_neo2_field failed: dB_dx"
            print *, "dB_dx: ", dB_dx
            print *, "NEO-2: ", dB_dx_neo2(case, :)
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_neo2_field
