program test_anti_sigma_analytic
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: surface_average_t, calc_surface_averages
    use fieldline_integrals, only: fourier_transform_over_label
    use fieldline_integrals, only: fieldline_modes_t
    use utils, only: linspace
    use utils, only: not_same

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = 0.0125_dp, eps_1 = 0.0005_dp
    real(dp), parameter :: I_v_1 = -8.0_dp*eps_1*sqrt(2*eps_0)/(B_0**2*N_tor)
    real(dp), parameter :: delta_A_o_1 = 0.25_dp*eps_1/eps_0
    real(dp), parameter :: average_lambda_b = sqrt(8.0_dp*eps_0)/pi
    real(dp), parameter :: phi_0 = pi
    type(anti_sigma_field_t) :: field

    real(dp), parameter :: phi_tol = 4e-6
    real(dp), parameter :: reltol_radial_drift = 2.0_dp*eps_0
    real(dp), parameter :: reltol_aspect_ratio = 2.0_dp*eps_1/eps_0
    real(dp), parameter :: reltol_average_lambda_b = 2.0_dp*eps_0
    real(dp), parameter :: reltol_average_B_squared = 2.0_dp*eps_0**2
    real(dp), parameter :: abstol = 1e-13
    integer, parameter :: n_fieldlines = 60

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.0_dp ! analytic formula for small iota
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    type(fieldline_modes_t) :: fieldline_modes

    type(surface_average_t) :: surface_average

    integer :: current
    real(dp) :: B_mod
    real(dp), dimension(:), allocatable :: zeros

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, theta_0, iota, field, M_pol, N_tor, &
                                  phi_tol)

    do current = 1, n_fieldlines
        if (not_same(phi_0, modulo(fieldlines(current)%phi_0, 2.0_dp*pi), &
                     abstol_in=2.0_dp*phi_tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_anti_sigma_analytic failed: phi_0"
            print *, "found: ", fieldlines(current)%phi_0
            print *, "expected: ", phi_0
            error stop
        end if
    end do

    call fourier_transform_over_label(fieldlines, fieldline_modes)

    if (not_same(I_v_1, fieldline_modes%radial_drift%sin_coeffs(2), &
                 reltol_in=reltol_radial_drift, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: 1st radial drift sin mode"
        print *, "found: ", fieldline_modes%radial_drift%sin_coeffs(2)
        print *, "expected: ", I_v_1
        print *, "ratio: ", fieldline_modes%radial_drift%sin_coeffs(2)/I_v_1
        error stop
    end if

    if (not_same(delta_A_o_1, fieldline_modes%delta_aspect_ratio%cos_coeffs(2), &
                 reltol_in=reltol_aspect_ratio, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: 1st delta_aspect_ratio cos mode"
        print *, "found: ", fieldline_modes%delta_aspect_ratio%cos_coeffs(2)
        print *, "expected: ", delta_A_o_1
        print *, "ratio: ", fieldline_modes%delta_aspect_ratio%cos_coeffs(2)/delta_A_o_1
        error stop
    end if

    allocate (zeros(size(fieldline_modes%radial_drift%cos_coeffs)))
    zeros = 0.0_dp

    if (not_same(zeros, fieldline_modes%radial_drift%cos_coeffs, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: radial drift cos modes"
        print *, "found: ", fieldline_modes%radial_drift%cos_coeffs
        print *, "expected: all ", zeros(1)
        error stop
    end if

    if (not_same(zeros, fieldline_modes%delta_aspect_ratio%sin_coeffs, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: delta_aspect_ratio sin modes"
        print *, "found: ", fieldline_modes%delta_aspect_ratio%sin_coeffs
        print *, "expected: all ", zeros(1)
        error stop
    end if

    if (not_same(zeros, fieldline_modes%delta_eta%cos_coeffs, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: delta_eta cos modes"
        print *, "found: ", fieldline_modes%delta_eta%cos_coeffs
        print *, "expected: all ", zeros(1)
        error stop
    end if

    if (not_same(zeros, fieldline_modes%delta_eta%sin_coeffs, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: delta_eta sin modes"
        print *, "found: ", fieldline_modes%delta_eta%sin_coeffs
        print *, "expected: all ", zeros(1)
        error stop
    end if

    deallocate (zeros)

    call calc_surface_averages(fieldlines, surface_average)

    if (not_same(average_lambda_b, surface_average%lambda_b, &
                 reltol_in=reltol_average_lambda_b, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: surface average lambda_b"
        print *, "found: ", surface_average%lambda_b
        print *, "expected: ", average_lambda_b
        print *, "ratio: ", surface_average%lambda_b/average_lambda_b
        error stop
    end if

    if (not_same(B_0**2, surface_average%B_squared, &
                 reltol_in=reltol_average_B_squared, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_anti_sigma_analytic failed: surface average lambda_b"
        print *, "found: ", surface_average%B_squared
        print *, "expected: ", B_0**2
        print *, "ratio: ", surface_average%B_squared/B_0**2
        error stop
    end if

end program test_anti_sigma_analytic
