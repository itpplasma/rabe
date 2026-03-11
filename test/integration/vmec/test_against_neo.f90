!> Test that boozer_field_t produces the same results as neo_field_t
!> for the same equilibrium on multiple flux surfaces.
!> Also writes all neo reference values to neo_reference_values.dat.
program test_against_neo
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use neo_field, only: neo_field_t
    use utils, only: not_same

    implicit none

    real(dp) :: reltol = 1e-4, abstol = 1e-8
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
    character(len=*), parameter :: bc_filename = &
                                   "input/landreman_paul_qa.bc"

    type(boozer_field_t) :: bfield
    type(neo_field_t) :: nfield

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)

    real(dp) :: bmod_b, sqrtg_b, dB_dx_b(3)
    real(dp) :: bmod_n, sqrtg_n, dB_dx_n(3)
    real(dp) :: iota_b, B_theta_b, B_phi_b

    integer :: case, u
    logical :: test_failed

    call bfield%boozer_field_init(nc_filename, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=3)
    call nfield%neo_field_init(bc_filename, stor=0.5_dp)

    test_failed = .false.
    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta = [0.0_dp, 1.0_dp, 3.14_dp, 0.5_dp, 2.0_dp]
    phi = [0.0_dp, 0.5_dp, 1.57_dp, 2.0_dp, 3.0_dp]

    open (newunit=u, file='neo_reference_values.dat', status='replace')
    write (u, '(A)') '# Neo reference values for test_against_neo'
    write (u, '(A)') '# psi_tor_edge, nfp (global)'
    write (u, '(A, ES25.17)') 'psi_tor_edge = ', nfield%psi_tor_edge
    write (u, '(A, ES25.17)') 'nfp          = ', nfield%nfp
    write (u, '(A)') ''

    do case = 1, n_cases
        call nfield%neo_change_stor(stor(case))
        call bfield%fix_to_surface(stor(case))

        call bfield%get_iota_and_covariant_components( &
            stor(case), iota_b, B_theta_b, B_phi_b)

        call bfield%compute_B_sqrtg_dB_dx(theta(case), phi(case), &
                                          bmod_b, sqrtg_b, dB_dx_b)
        call nfield%compute_B_sqrtg_dB_dx(theta(case), phi(case), &
                                          bmod_n, sqrtg_n, dB_dx_n)

        write (u, '(A, I1)') '# case ', case
        write (u, '(A, F4.1, A, F5.2, A, F5.2)') &
            '# stor=', stor(case), &
            ' theta=', theta(case), ' phi=', phi(case)
        write (u, '(A, ES25.17)') 'iota               = ', nfield%iota
        write (u, '(A, ES25.17)') 'B_theta_covariant  = ', &
            nfield%B_theta_covariant
        write (u, '(A, ES25.17)') 'B_phi_covariant    = ', &
            nfield%B_phi_covariant
        write (u, '(A, ES25.17)') 'bmod               = ', bmod_n
        write (u, '(A, ES25.17)') 'sqrtg              = ', sqrtg_n
        write (u, '(A, 3ES25.17)') 'dB_dx              = ', dB_dx_n
        write (u, '(A)') ''

        if (not_same(iota_b, nfield%iota, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "iota mismatch at case ", case
            print *, "  boozer: ", iota_b
            print *, "  neo:    ", nfield%iota
            test_failed = .true.
        end if

        if (not_same(B_theta_b, nfield%B_theta_covariant, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_theta_covariant mismatch at case ", case
            print *, "  boozer: ", B_theta_b
            print *, "  neo:    ", nfield%B_theta_covariant
            test_failed = .true.
        end if

        if (not_same(B_phi_b, nfield%B_phi_covariant, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_phi_covariant mismatch at case ", case
            print *, "  boozer: ", B_phi_b
            print *, "  neo:    ", nfield%B_phi_covariant
            test_failed = .true.
        end if

        if (not_same(bfield%psi_tor_edge, nfield%psi_tor_edge, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "psi_tor_edge mismatch at case ", case
            print *, "  boozer: ", bfield%psi_tor_edge
            print *, "  neo:    ", nfield%psi_tor_edge
            test_failed = .true.
        end if

        if (not_same(bfield%nfp, nfield%nfp, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "nfp mismatch at case ", case
            print *, "  boozer: ", bfield%nfp
            print *, "  neo:    ", nfield%nfp
            test_failed = .true.
        end if

        if (not_same(bmod_b, bmod_n, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_mod mismatch at case ", case
            print *, "  boozer: ", bmod_b
            print *, "  neo:    ", bmod_n
            test_failed = .true.
        end if

        if (not_same(sqrtg_b, sqrtg_n, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "sqrtg mismatch at case ", case
            print *, "  boozer: ", sqrtg_b
            print *, "  neo:    ", sqrtg_n
            test_failed = .true.
        end if

        if (not_same(dB_dx_b, dB_dx_n, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "dB_dx mismatch at case ", case
            print *, "  boozer: ", dB_dx_b
            print *, "  neo:    ", dB_dx_n
            test_failed = .true.
        end if
    end do

    close (u)
    print *, "Neo reference values written to neo_reference_values.dat"

    if (test_failed) error stop

end program test_against_neo
