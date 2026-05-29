program test_find_global_maximum
    use constants, only: dp, pi
    use utils, only: not_same, linspace
    use mock_field, only: mock_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines

    implicit none

    real(dp), parameter :: theta_mode = -1.0_dp, phi_mode = 4.0_dp
    real(dp), parameter :: nfp = max(1.0_dp, abs(phi_mode))
    real(dp), parameter :: B_0 = 1.0_dp, B_amplitude = -0.5_dp, B_pert = -0.15_dp
    real(dp), parameter :: global_B_max = B_0 + abs(B_amplitude) + abs(B_pert)
    type(mock_field_t) :: field
    type(mock_perturbed_field_t) :: perturbed_field

    real(dp) :: tol_phi_max
    real(dp) :: B_reltol

    integer, parameter :: n_fieldlines = 100
    real(dp), parameter :: iota = -1.0_dp

    real(dp), dimension(n_fieldlines) :: theta_0
    type(flock_of_fieldlines_t) :: flock

    real(dp) :: found_global_B_max, fieldline_B_max

    call field%mock_field_init(theta_mode, phi_mode, B_0, B_amplitude)
    call perturbed_field%mock_perturbed_field_init(field, theta_mode, 0.0_dp, B_pert)

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)
    call make_flock_of_fieldlines(flock, theta_0, iota, perturbed_field, &
                                  theta_mode, phi_mode, nfp)

    tol_phi_max = max(maxval(flock%fieldlines%phi_max_error(1)), &
                      maxval(flock%fieldlines%phi_max_error(2)))
    B_reltol = max((tol_phi_max*abs(phi_mode))**2.0_dp, &
                   (2.0_dp*pi/n_fieldlines)**2.0_dp*1e-3)

    found_global_B_max = 1.0_dp/flock%eta_b
    if (not_same(global_B_max, found_global_B_max, B_reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_find_global_maximum failed: global_B_max"
        print *, "found: ", found_global_B_max
        print *, "expected: ", global_B_max
        print *, "relative error: ", 1.0_dp - found_global_B_max/global_B_max
        error stop
    end if

end program test_find_global_maximum
