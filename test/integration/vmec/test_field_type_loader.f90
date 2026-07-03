!> Re-initialization consistency check for boozer_field_t.
!>
!> The Boozer splines live in module-global state (singleton). This test loads
!> the same VMEC wout file into two boozer_field_t instances and asserts that
!> the second init reproduces the metadata (nfp, psi_tor_edge, R) of the first.
program test_field_type_loader
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    type(boozer_field_t) :: field_vmec
    type(boozer_field_t) :: field_vmec_ref

    character(len=*), parameter :: vmec_nc_file = &
        'input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc'

    real(dp), parameter :: abstol = 1e-9_dp, reltol = 1e-9_dp
    logical :: test_failed

    test_failed = .false.

    call field_vmec%init_from_vmec(vmec_nc_file, grid_refinement=6)
    if (.not. field_vmec%initialized) then
        print *, "-------------------------------------------------------------"
        print *, "test_field_type_loader failed: first field not initialized"
        test_failed = .true.
    end if

    call field_vmec_ref%init_from_vmec(vmec_nc_file, grid_refinement=6)
    if (.not. field_vmec_ref%initialized) then
        print *, "-------------------------------------------------------------"
        print *, "test_field_type_loader failed: second field not initialized"
        test_failed = .true.
    end if

    if (not_same(field_vmec%nfp, field_vmec_ref%nfp, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_field_type_loader failed: nfp mismatch"
        print *, "first init: ", field_vmec%nfp
        print *, "second init: ", field_vmec_ref%nfp
        test_failed = .true.
    end if

    if (not_same(field_vmec%psi_tor_edge, field_vmec_ref%psi_tor_edge, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_field_type_loader failed: psi_tor_edge mismatch"
        print *, "first init: ", field_vmec%psi_tor_edge
        print *, "second init: ", field_vmec_ref%psi_tor_edge
        test_failed = .true.
    end if

    if (not_same(field_vmec%R, field_vmec_ref%R, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_field_type_loader failed: R mismatch"
        print *, "first init: ", field_vmec%R
        print *, "second init: ", field_vmec_ref%R
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_field_type_loader
