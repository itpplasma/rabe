!> Round trip: VMEC wout -> boozer_field_t -> write .bc -> bc_field_t.
!> Asserts that B, its derivatives, sqrt(g), iota, and the covariant
!> components from the .bc file agree with the originating Boozer field.
!> The .bc writer emits dummy R/Z/v harmonics, so nabla_s is not compared.
program test_bc_roundtrip
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use boozer_field, only: boozer_field_t
    use bc_field, only: bc_field_t
    use bc_file, only: write_field_B_mod_to_bc, delete_bc_file

    implicit none

    character(len=*), parameter :: nc_filename = &
        "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    character(len=*), parameter :: bc_filename = "roundtrip.bc"

    integer, parameter :: n_surf = 9
    integer, parameter :: m_max = 16, n_max = 16
    real(dp), parameter :: s_test = 0.5_dp
    integer, parameter :: n_theta = 17, n_phi = 13

    type(boozer_field_t) :: field_nc
    type(bc_field_t) :: field_bc

    real(dp) :: s_tors(n_surf)
    real(dp) :: theta, phi
    real(dp) :: B_nc, B_bc, sqrtg_nc, sqrtg_bc
    real(dp) :: dB_dx_nc(3), dB_dx_bc(3)
    real(dp) :: iota_nc, iota_bc
    real(dp) :: B_theta_cov_nc, B_phi_cov_nc, B_theta_cov_bc, B_phi_cov_bc
    real(dp) :: max_rel_diff_B, max_rel_diff_sqrtg
    real(dp) :: max_diff_dB_ang, max_abs_dB_ang
    real(dp) :: max_diff_dB_ds, max_abs_dB_ds
    integer :: i_theta, i_phi
    logical :: test_failed

    real(dp), parameter :: tol_B = 1e-5_dp
    real(dp), parameter :: tol_sqrtg = 1e-5_dp
    real(dp), parameter :: tol_dB_ang = 1e-3_dp !< relative to max |dB_ang|
    real(dp), parameter :: tol_dB_ds = 1e-2_dp !< relative to max |dB/ds|
    real(dp), parameter :: tol_flux_funcs = 1e-7_dp

    test_failed = .false.

    call field_nc%boozer_field_init(nc_filename, &
                                    radial_spline_order=5, &
                                    angular_spline_order=5, &
                                    grid_refinement=3)
    call linspace(0.1_dp, 0.9_dp, n_surf, s_tors)
    call delete_bc_file(bc_filename)
    call write_field_B_mod_to_bc(field_nc, s_tors, m_max, n_max, bc_filename)

    call field_bc%bc_field_init(bc_filename)

    if (not_same(field_bc%psi_tor_edge, field_nc%psi_tor_edge, &
                 reltol_in=tol_flux_funcs)) then
        print *, "psi_tor_edge mismatch: ", field_bc%psi_tor_edge, &
            field_nc%psi_tor_edge
        test_failed = .true.
    end if
    if (not_same(field_bc%nfp, field_nc%nfp, reltol_in=1e-14_dp)) then
        print *, "nfp mismatch: ", field_bc%nfp, field_nc%nfp
        test_failed = .true.
    end if
    if (not_same(field_bc%R, field_nc%R, reltol_in=tol_flux_funcs)) then
        print *, "R mismatch: ", field_bc%R, field_nc%R
        test_failed = .true.
    end if

    call field_nc%fix_to_surface(s_test)
    call field_bc%fix_to_surface(s_test)

    call field_nc%get_iota(s_test, iota_nc)
    call field_bc%get_iota(s_test, iota_bc)
    if (not_same(iota_bc, iota_nc, reltol_in=tol_flux_funcs)) then
        print *, "iota mismatch: ", iota_bc, iota_nc
        test_failed = .true.
    end if

    call field_nc%get_covariant_components(B_theta_cov_nc, B_phi_cov_nc)
    call field_bc%get_covariant_components(B_theta_cov_bc, B_phi_cov_bc)
    if (not_same(B_theta_cov_bc, B_theta_cov_nc, reltol_in=tol_flux_funcs) &
        .or. not_same(B_phi_cov_bc, B_phi_cov_nc, &
                      reltol_in=tol_flux_funcs)) then
        print *, "covariant components mismatch: ", &
            B_theta_cov_bc, B_theta_cov_nc, B_phi_cov_bc, B_phi_cov_nc
        test_failed = .true.
    end if

    max_rel_diff_B = 0.0_dp
    max_rel_diff_sqrtg = 0.0_dp
    max_diff_dB_ang = 0.0_dp
    max_abs_dB_ang = 0.0_dp
    max_diff_dB_ds = 0.0_dp
    max_abs_dB_ds = 0.0_dp

    do i_phi = 0, n_phi - 1
        do i_theta = 0, n_theta - 1
            theta = real(i_theta, dp)*2.0_dp*pi/real(n_theta, dp)
            phi = real(i_phi, dp)*2.0_dp*pi/(field_nc%nfp*real(n_phi, dp))

            call field_nc%compute_B_sqrtg_dB_dx(theta, phi, B_nc, &
                                                sqrtg_nc, dB_dx_nc)
            call field_bc%compute_B_sqrtg_dB_dx(theta, phi, B_bc, &
                                                sqrtg_bc, dB_dx_bc)

            max_rel_diff_B = max(max_rel_diff_B, abs(B_bc - B_nc)/abs(B_nc))
            max_rel_diff_sqrtg = max(max_rel_diff_sqrtg, &
                                     abs(sqrtg_bc - sqrtg_nc)/abs(sqrtg_nc))
            max_diff_dB_ang = max(max_diff_dB_ang, &
                                  maxval(abs(dB_dx_bc(2:3) - dB_dx_nc(2:3))))
            max_abs_dB_ang = max(max_abs_dB_ang, maxval(abs(dB_dx_nc(2:3))))
            max_diff_dB_ds = max(max_diff_dB_ds, &
                                 abs(dB_dx_bc(1) - dB_dx_nc(1)))
            max_abs_dB_ds = max(max_abs_dB_ds, abs(dB_dx_nc(1)))
        end do
    end do

    print *, "max rel diff B_mod:           ", max_rel_diff_B
    print *, "max rel diff sqrtg:           ", max_rel_diff_sqrtg
    print *, "max diff dB/d(theta,phi):     ", max_diff_dB_ang, &
        " of max |dB_ang| ", max_abs_dB_ang
    print *, "max diff dB/ds:               ", max_diff_dB_ds, &
        " of max |dB/ds| ", max_abs_dB_ds

    if (max_rel_diff_B > tol_B) then
        print *, "B_mod round trip exceeds tolerance ", tol_B
        test_failed = .true.
    end if
    if (max_rel_diff_sqrtg > tol_sqrtg) then
        print *, "sqrtg round trip exceeds tolerance ", tol_sqrtg
        test_failed = .true.
    end if
    if (max_diff_dB_ang > tol_dB_ang*max_abs_dB_ang) then
        print *, "dB/d(theta,phi) round trip exceeds tolerance"
        test_failed = .true.
    end if
    if (max_diff_dB_ds > tol_dB_ds*max_abs_dB_ds) then
        print *, "dB/ds round trip exceeds tolerance"
        test_failed = .true.
    end if

    if (test_failed) error stop "test_bc_roundtrip failed"
    print *, "test_bc_roundtrip passed"

end program test_bc_roundtrip
