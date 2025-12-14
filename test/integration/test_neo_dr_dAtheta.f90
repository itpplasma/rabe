program test_neo_dr_dAtheta
    use neo_field, only: neo_field_t
    use constants, only: dp, pi
    use utils, only: not_same, linspace
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use surface_average_mod, only: surface_average_t, calc_surface_averages

    implicit none

    real(dp), parameter :: reltol = 2e-4
    real(dp), parameter :: phi_tol = 6e-7
    character(len=*), parameter :: bc_filename = "input/quasi_helical.bc"

    type(neo_field_t) :: field
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 4.0_dp
    real(dp), parameter :: stor = 0.25_dp
    integer, parameter :: n_fieldlines = 100

    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    real(dp) :: iota, nfp
    real(dp) :: xi_0(n_fieldlines), temp(n_fieldlines + 1)
    type(surface_average_t) :: average
    real(dp) :: dr_dAtheta, neo_dr_dAtheta

    real(dp), parameter :: neo_ds_dr = 0.00807615_dp ![1/cm]
    real(dp), parameter :: angle_screw_sign = -1.0_dp

    real(dp) :: sqrt_g11
    logical :: test_failed

    test_failed = .false.

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    nfp = field%nfp
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  xi_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  nfp, &
                                  phi_tol)

    call calc_surface_averages(fieldlines, average)

    neo_dr_dAtheta = angle_screw_sign/(field%psi_tor_edge*neo_ds_dr*100.0_dp)

    dr_dAtheta = angle_screw_sign*sign(1.0_dp, field%psi_tor_edge)/average%sqrt_g11

    if (not_same(neo_dr_dAtheta, &
                 dr_dAtheta, &
                 reltol_in=reltol, &
                 abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_neo_dr_dAtheta failed: dr_dAtheta"
        print *, "expected: ", neo_dr_dAtheta
        print *, "found: ", dr_dAtheta
        print *, "relative difference: ", &
            abs(dr_dAtheta/neo_dr_dAtheta - 1.0_dp)
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_neo_dr_dAtheta
