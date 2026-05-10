!> Test that boozer_field_t produces the same results as the NEO code
!> for a .bc file that is generated directly from boozer_field_t
program test_against_neo_on_boozer_field
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use neo_field, only: neo_field_t
    use bc_file, only: write_field_B_mod_to_bc, delete_bc_file
    use utils, only: not_same
    use fieldline_mod, only: fieldline_t
    use fieldline_labels, only: get_labels
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    implicit none

    !> reltol for the final coefficients is way bigger than the actually difference
    !> in the field i.e. B (difference ~1e-10), because the coefficients themselves
    !> scale by the small violations of omnigeneity -> amplification of differences
    real(dp), parameter :: reltol = 1.e-4, abstol = 0.0_dp
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 0.0_dp
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield
    type(neo_field_t) :: nfield

    integer, parameter :: n_stor = 8
    integer :: idx
    real(dp), dimension(n_stor) :: stors
    character(len=*), parameter :: bc_filename = "test.bc"

    real(dp) :: stor

    integer, parameter :: max_n_fieldlines = 20
    type(fieldline_t), dimension(:), allocatable :: fieldlines
    real(dp), dimension(:), allocatable :: xi_0
    real(dp) :: iota, approx_iota, B_theta_covariant, B_phi_covariant
    real(dp) :: covariant_factor, R
    real(dp), parameter :: dr_dAtheta = 1.0_dp ! dummy value, same for both field case
    real(dp), dimension(n_stor) :: neo_Lambda_bl, neo_Lambda_lm
    real(dp), dimension(n_stor) :: bfield_Lambda_bl, bfield_Lambda_lm

    real(dp) :: deviation_A, deviation_B

    logical :: test_failed

    test_failed = .false.

    do idx = 1, n_stor
        stors(idx) = (real(idx, dp) - 0.5_dp)/real(n_stor, dp)
    end do

    call bfield%boozer_field_init(nc_filename, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=12)
    call write_field_B_mod_to_bc(bfield, stors, m_max=31, n_max=32, &
                                 filename=bc_filename)

    call nfield%neo_field_init(bc_filename, stors(1))
    ! We only compare the inner surfaces (away from axis and edge)
    ! this restricts evaluation to regions of not too strong omniegenity violation
    ! so that the calculation can execute without crashing
    do idx = 2, 6
        stor = stors(idx)
        print *, "Testing surface with s_tor =", stor
        print *, "Boozer field"
        call bfield%fix_to_surface(stor)

        call bfield%get_iota_and_covariant_components(stor, &
                                                      iota, &
                                                      B_theta_covariant, &
                                                      B_phi_covariant)
        call get_labels(max_n_fieldlines, iota, M_pol, N_tor, bfield%nfp, &
                        xi_0, approx_iota)
        allocate (fieldlines(size(xi_0)))
        call make_flock_of_fieldlines(fieldlines, &
                                      xi_0, &
                                      approx_iota, &
                                      bfield, &
                                      M_pol, &
                                      N_tor, &
                                      bfield%nfp)
        call calc_deviation(fieldlines, deviation_A, deviation_B)
        R = bfield%R
        covariant_factor = B_phi_covariant + approx_iota*B_theta_covariant
        bfield_Lambda_bl(idx) = deviation_A*dr_dAtheta* &
                                sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
        bfield_Lambda_lm(idx) = deviation_B*0.5*R*pi*dr_dAtheta
        deallocate (fieldlines)

        print *, "NEO field"
        call nfield%neo_change_stor(stor)
        iota = nfield%iota
        B_theta_covariant = nfield%B_theta_covariant
        B_phi_covariant = nfield%B_phi_covariant
        call get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfield%nfp, &
                        xi_0, approx_iota)
        allocate (fieldlines(size(xi_0)))
        call make_flock_of_fieldlines(fieldlines, &
                                      xi_0, &
                                      approx_iota, &
                                      nfield, &
                                      M_pol, &
                                      N_tor, &
                                      nfield%nfp)
        call calc_deviation(fieldlines, deviation_A, deviation_B)
        R = nfield%R
        covariant_factor = B_phi_covariant + approx_iota*B_theta_covariant
        neo_Lambda_bl(idx) = deviation_A*dr_dAtheta* &
                             sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
        neo_Lambda_lm(idx) = deviation_B*0.5*R*pi*dr_dAtheta
        deallocate (fieldlines)
    end do

    do idx = 2, 6
        if (not_same(bfield_Lambda_bl(idx), neo_Lambda_bl(idx), reltol, abstol)) then
            print *, "Lambda_bl mismatch at s_tor =", stors(idx)
            print *, "bfield Lambda_bl =", bfield_Lambda_bl(idx)
            print *, "neo Lambda_bl =", neo_Lambda_bl(idx)
            print *, "Relative difference =", abs(bfield_Lambda_bl(idx) &
                                                  - neo_Lambda_bl(idx)) &
                /abs(neo_Lambda_bl(idx))
            test_failed = .true.
        end if
        if (not_same(bfield_Lambda_lm(idx), neo_Lambda_lm(idx), reltol, abstol)) then
            print *, "Lambda_lm mismatch at s_tor =", stors(idx)
            print *, "bfield Lambda_lm =", bfield_Lambda_lm(idx)
            print *, "neo Lambda_lm =", neo_Lambda_lm(idx)
            print *, "Relative difference =", abs(bfield_Lambda_lm(idx) - &
                                                  neo_Lambda_lm(idx)) &
                /abs(neo_Lambda_lm(idx))
            test_failed = .true.
        end if
    end do

    call delete_bc_file(bc_filename)

    if (test_failed) error stop

end program test_against_neo_on_boozer_field
