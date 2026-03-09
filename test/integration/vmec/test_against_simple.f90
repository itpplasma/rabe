!> Test that the field implementation matches the SIMPLE code for a few points in the plasma.
program test_against_simple
    use constants, only: dp
    use new_vmec_stuff_mod, only: netcdffile, ns_s, ns_tp, multharm
    use spline_vmec_sub, only: spline_vmec_data
    use boozer_sub, only: get_boozer_coordinates
    use magfie_sub, only: init_magfie, BOOZER

    implicit none

    real(dp) :: retol = 1e-11
    character(len=*), parameter :: nc_filename = "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: bmod, sqrtg, dB_dx(3)
    real(dp) :: bmod_simple(n_cases), sqrtg_simple(n_cases), dB_dx_simple(n_cases, 3)
    integer :: case

    logical :: test_failed

    ! 1. Init VMEC splines
    netcdffile = nc_filename
    ns_s = 4
    ns_tp = 4
    multharm = 6
    call spline_vmec_data()
    call get_boozer_coordinates()
    call init_magfie(BOOZER) test_failed = .false.

    stor = [0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp]
    theta = [-1.00_dp, 0.00_dp, 4.3498975203550980_dp, 1.00_dp]
    phi = [1.00_dp, -1.00_dp, 0.56947933720590893_dp, 0.00_dp]

    bmod_simple = [5.8461732541782538_dp, &
                   6.3747018649171556_dp, &
                   6.3402245795215624_dp, &
                   5.4298976321806043_dp]
    sqrtg_simple = [-19930779.196453653_dp, &
                    -16763094.050708901_dp, &
                    -16945898.197301224_dp, &
                    -23104179.862998683_dp]
    dB_dx_simple(1,:) = [-0.92227365739728728_dp,       0.13109042596457662_dp,      -0.52436170385830649_dp]
    dB_dx_simple(2,:) = [0.53462460260148859_dp,      -0.55835449392655978_dp,        2.2334179757062391_dp]
dB_dx_simple(3, :) = [0.340429548051560_dp, 0.784044027926972_dp, -3.136176111707888_dp]
    dB_dx_simple(4,:) = [-0.20339965461773155_dp,       0.77633188889709603_dp,       -3.1053275555883841_dp]

    call field%neo_field_init(nc_filename)
    do case = 1, n_cases
        if (abs(bmod/bmod_simple(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_simple failed: B"
            print *, "B: ", bmod
            print *, "SIMPLE: ", bmod_simple(case)
            test_failed = .true.
        end if
        if (abs(sqrtg/sqrtg_simple(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_simple failed: sqrtg"
            print *, "sqrtg: ", sqrtg
            print *, "SIMPLE: ", sqrtg_simple(case)
            test_failed = .true.
        end if
        if (any(abs(dB_dx/dB_dx_simple(case, :) - 1) > retol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_simple failed: dB_dx"
            print *, "dB_dx: ", dB_dx
            print *, "SIMPLE: ", dB_dx_simple(case, :)
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_simple
