program test_neo_covariant_components
    use neo_field, only: neo_field_t
    use constants, only: dp, pi
    use utils, only: not_same

    implicit none

    real(dp) :: reltol = 1e-3, abstol = 1e-16
    real(dp) :: reltol_psi = 1e-6
    character(len=*), parameter :: bc_filename = "input/quasi_helical.bc"
    real(dp), parameter :: N = 4.0_dp

    type(neo_field_t) :: field
    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases)
    real(dp) :: B_theta_covariant(n_cases), B_phi_covariant(n_cases), iota(n_cases)
    real(dp) :: psi_tor_edge

    integer :: case
    logical :: test_failed

    test_failed = .false.

    stor = (/0.02_dp, 0.50_dp, 0.75_dp, 0.98_dp/)
    B_theta_covariant = (/7.45e-10, 3.1e-09, -2.84e-09, 1.17e-08/)*2.0_dp*1e-7
    B_phi_covariant = (/-1.277e8, -1.277e8, -1.277e8, -1.277e8/)*N*2.0_dp*1e-7
    iota = (/-1.238_dp, -1.243_dp, -1.244_dp, -1.245_dp/)
    psi_tor_edge = 41.86388_dp/(2.0_dp*pi)

    call field%neo_field_init(bc_filename, stor(1))
    do case = 1, n_cases - 1
        if (not_same(B_theta_covariant(case), &
                     field%B_theta_covariant, &
                     reltol_in=0.0_dp, &
                     abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_neo_covariant_components failed: B_theta_covariant"
            print *, "expected: ", B_theta_covariant(case)
            print *, "got: ", field%B_theta_covariant
            print *, "relative difference: ", &
                field%B_theta_covariant/B_theta_covariant(case) - 1.0_dp
            test_failed = .true.
        end if
        if (not_same(B_phi_covariant(case), &
                     field%B_phi_covariant, &
                     reltol_in=reltol, &
                     abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_neo_covariant_components failed: B_phi_covariant"
            print *, "expected: ", B_phi_covariant(case)
            print *, "got: ", field%B_phi_covariant
            print *, "relative difference: ", &
                field%B_phi_covariant/B_phi_covariant(case) - 1.0_dp
            test_failed = .true.
        end if
        if (not_same(iota(case), &
                     field%iota, &
                     reltol_in=reltol, &
                     abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_neo_covariant_components failed: iota"
            print *, "expected: ", iota(case)
            print *, "got: ", field%iota
            print *, "relative difference: ", &
                field%iota/iota(case) - 1.0_dp
            test_failed = .true.
        end if
        if (not_same(psi_tor_edge, &
                     field%psi_tor_edge, &
                     reltol_in=reltol_psi, &
                     abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_neo_covariant_components failed: psi_tor_edge"
            print *, "expected: ", psi_tor_edge
            print *, "got: ", field%psi_tor_edge
            print *, "relative difference: ", &
                field%psi_tor_edge/psi_tor_edge - 1.0_dp
            test_failed = .true.
        end if
        call field%neo_change_stor(stor(case + 1))
    end do

    if (test_failed) error stop

end program test_neo_covariant_components
