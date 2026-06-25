!> Cross-check test for booz_xform and chartmap field loaders.
!>
!> Both paths read a field generated from circ.bc via libneo Python converters.
!> The test asserts that:
!>   1. boozer_field_t initializes without error via field_type='chartmap'
!>   2. Bmod at three reference points agrees with the Python-computed grid values
!>      to within 1e-4 (fixture is FP-accurate after the off-by-one fix).
!>   3. boozer_field_t initializes without error via field_type='booz_xform'
!>   4. Bmod from the booz_xform path agrees with the chartmap path at the same
!>      points to within 1e-4 (both paths use the fixed converter).
!>
!> All reference values were computed by bc_to_booz_xform + booz_xform_to_boozer_chartmap
!> from circ.bc (a circular cross-section tokamak example).  The chartmap carries
!> nrho=20, ntheta=24, nzeta=24.
program test_booz_xform_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    type(boozer_field_t) :: field_chart, field_booz

    character(len=*), parameter :: chartmap_file = &
        'input/circ_chartmap.nc'
    character(len=*), parameter :: boozmn_file = &
        'input/circ_boozmn.nc'

    ! Reference Bmod values in Tesla at (rho, theta, zeta) grid points ir=10,it=13,iz=13
    ! and ir=5,it=7,iz=7 and ir=15,it=19,iz=7 of the 20x24x24 grid.
    ! Source: circ.bc -> circ_boozmn.nc -> circ_chartmap.nc via fixed converter.
    real(dp), parameter :: Bmod_ref_1 = 2.2254159391e+00_dp  ! T (ir=10,it=13,iz=13)
    real(dp), parameter :: Bmod_ref_2 = 1.9672547079e+00_dp  ! T (ir=5,it=7,iz=7)
    real(dp), parameter :: Bmod_ref_3 = 2.0602177739e+00_dp  ! T (ir=15,it=19,iz=7)

    ! Coordinates: (stor, theta, zeta) where stor = rho^2
    real(dp), parameter :: stor_1 = 0.224876_dp
    real(dp), parameter :: theta_1 = 3.1416_dp
    real(dp), parameter :: zeta_1 = 3.1416_dp

    real(dp), parameter :: stor_2 = 0.044654_dp
    real(dp), parameter :: theta_2 = 1.5708_dp
    real(dp), parameter :: zeta_2 = 1.5708_dp

    real(dp), parameter :: stor_3 = 0.543324_dp
    real(dp), parameter :: theta_3 = 4.7124_dp
    real(dp), parameter :: zeta_3 = 1.5708_dp

    ! 1e-4 relative tolerance: on-grid chartmap evaluation from 20-point rho grid
    real(dp), parameter :: reltol_chart = 1.0e-4_dp
    ! 1e-4 relative tolerance: boozmn vs chartmap cross-check (both paths use
    ! the same fixed converter; off-by-one fix brings agreement to ~1e-8)
    real(dp), parameter :: reltol_booz = 1.0e-4_dp
    real(dp), parameter :: abstol = 1.0e-10_dp

    real(dp) :: Bmod_chart_1, Bmod_chart_2, Bmod_chart_3
    real(dp) :: Bmod_booz_1, Bmod_booz_2, Bmod_booz_3
    logical :: test_failed

    test_failed = .false.

    ! ----------------------------------------------------------------
    ! Load chartmap field
    ! ----------------------------------------------------------------
    call field_chart%boozer_field_init(chartmap_file, field_type='chartmap')

    if (.not. field_chart%initialized) then
        print *, "FAIL: chartmap field not initialized"
        test_failed = .true.
    end if

    call field_chart%fix_to_surface(stor_1)
    call field_chart%compute_B_mod(theta_1, zeta_1, Bmod_chart_1)

    call field_chart%fix_to_surface(stor_2)
    call field_chart%compute_B_mod(theta_2, zeta_2, Bmod_chart_2)

    call field_chart%fix_to_surface(stor_3)
    call field_chart%compute_B_mod(theta_3, zeta_3, Bmod_chart_3)

    if (not_same(Bmod_chart_1, Bmod_ref_1, reltol_in=reltol_chart, abstol_in=abstol)) then
        print *, "------------------------------------------------------------"
        print *, "FAIL: chartmap Bmod_1 mismatch"
        print *, "  computed: ", Bmod_chart_1
        print *, "  reference:", Bmod_ref_1
        print *, "  reldiff:  ", abs(Bmod_chart_1 - Bmod_ref_1)/Bmod_ref_1
        test_failed = .true.
    end if

    if (not_same(Bmod_chart_2, Bmod_ref_2, reltol_in=reltol_chart, abstol_in=abstol)) then
        print *, "------------------------------------------------------------"
        print *, "FAIL: chartmap Bmod_2 mismatch"
        print *, "  computed: ", Bmod_chart_2
        print *, "  reference:", Bmod_ref_2
        print *, "  reldiff:  ", abs(Bmod_chart_2 - Bmod_ref_2)/Bmod_ref_2
        test_failed = .true.
    end if

    if (not_same(Bmod_chart_3, Bmod_ref_3, reltol_in=reltol_chart, abstol_in=abstol)) then
        print *, "------------------------------------------------------------"
        print *, "FAIL: chartmap Bmod_3 mismatch"
        print *, "  computed: ", Bmod_chart_3
        print *, "  reference:", Bmod_ref_3
        print *, "  reldiff:  ", abs(Bmod_chart_3 - Bmod_ref_3)/Bmod_ref_3
        test_failed = .true.
    end if

    print *, "chartmap Bmod_1:", Bmod_chart_1, " ref:", Bmod_ref_1
    print *, "chartmap Bmod_2:", Bmod_chart_2, " ref:", Bmod_ref_2
    print *, "chartmap Bmod_3:", Bmod_chart_3, " ref:", Bmod_ref_3

    ! ----------------------------------------------------------------
    ! Load booz_xform field
    ! ----------------------------------------------------------------
    call field_booz%boozer_field_init(boozmn_file, field_type='booz_xform')

    if (.not. field_booz%initialized) then
        print *, "FAIL: booz_xform field not initialized"
        test_failed = .true.
    end if

    call field_booz%fix_to_surface(stor_1)
    call field_booz%compute_B_mod(theta_1, zeta_1, Bmod_booz_1)

    call field_booz%fix_to_surface(stor_2)
    call field_booz%compute_B_mod(theta_2, zeta_2, Bmod_booz_2)

    call field_booz%fix_to_surface(stor_3)
    call field_booz%compute_B_mod(theta_3, zeta_3, Bmod_booz_3)

    if (not_same(Bmod_booz_1, Bmod_ref_1, reltol_in=reltol_booz, abstol_in=abstol)) then
        print *, "------------------------------------------------------------"
        print *, "FAIL: booz_xform Bmod_1 mismatch vs chartmap reference"
        print *, "  computed: ", Bmod_booz_1
        print *, "  reference:", Bmod_ref_1
        print *, "  reldiff:  ", abs(Bmod_booz_1 - Bmod_ref_1)/Bmod_ref_1
        test_failed = .true.
    end if

    if (not_same(Bmod_booz_2, Bmod_ref_2, reltol_in=reltol_booz, abstol_in=abstol)) then
        print *, "------------------------------------------------------------"
        print *, "FAIL: booz_xform Bmod_2 mismatch vs chartmap reference"
        print *, "  computed: ", Bmod_booz_2
        print *, "  reference:", Bmod_ref_2
        print *, "  reldiff:  ", abs(Bmod_booz_2 - Bmod_ref_2)/Bmod_ref_2
        test_failed = .true.
    end if

    if (not_same(Bmod_booz_3, Bmod_ref_3, reltol_in=reltol_booz, abstol_in=abstol)) then
        print *, "------------------------------------------------------------"
        print *, "FAIL: booz_xform Bmod_3 mismatch vs chartmap reference"
        print *, "  computed: ", Bmod_booz_3
        print *, "  reference:", Bmod_ref_3
        print *, "  reldiff:  ", abs(Bmod_booz_3 - Bmod_ref_3)/Bmod_ref_3
        test_failed = .true.
    end if

    print *, "booz_xform Bmod_1:", Bmod_booz_1, " ref:", Bmod_ref_1
    print *, "booz_xform Bmod_2:", Bmod_booz_2, " ref:", Bmod_ref_2
    print *, "booz_xform Bmod_3:", Bmod_booz_3, " ref:", Bmod_ref_3

    if (test_failed) error stop

end program test_booz_xform_chartmap
