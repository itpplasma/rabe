!> Test that boozer_field_t produces the same results as the NEO code
!> The reference values were obtained from a field evalution of a .bc file
!> converted from the original .nc file
program test_against_neo_on_boozer_field
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use neo_field, only: neo_field_t
    use bc_file, only: write_field_B_mod_to_bc, delete_bc_file
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.e-9, abstol = 1.e-16
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield
    type(neo_field_t) :: nfield

    integer, parameter :: n_stor = 8
    integer :: idx
    real(dp), dimension(n_stor) :: stors
    character(len=*), parameter :: bc_filename = "test.bc"

    real(dp) :: stor, theta, phi
    real(dp) :: B_mod_bfield, B_mod_nfield
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

    theta = 0.3_dp
    phi = 0.7_dp
    call nfield%neo_field_init(bc_filename, stors(1))
    do idx = 2, n_stor
        stor = stors(idx)
        call bfield%fix_to_surface(stor)
        call bfield%compute_B_mod(theta, phi, B_mod_bfield)

        call nfield%neo_change_stor(stor)
        call nfield%compute_B_mod(theta, phi, B_mod_nfield)
        if (not_same(B_mod_bfield, B_mod_nfield, reltol, abstol)) then
            print *, "B_mod mismatch at s_tor =", stor, "theta =", theta, "phi =", phi
            print *, "B_mod_bfield =", B_mod_bfield
            print *, "B_mod_nfield =", B_mod_nfield
            print *, "Relative difference =", abs(B_mod_bfield - B_mod_nfield) &
                /abs(B_mod_nfield)
            test_failed = .true.
        end if
    end do

    call delete_bc_file(bc_filename)

    if (test_failed) error stop

end program test_against_neo_on_boozer_field
