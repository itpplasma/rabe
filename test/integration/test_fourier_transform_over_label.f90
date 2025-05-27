program test_fourier_transform_over_label
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_mod, only: make_flock_of_fieldlines
    use fieldline_integrals, only: fourier_transform_over_label
    use fieldline_integrals, only: fieldline_modes_t
    use utils, only: linspace
    use utils, only: not_same

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = 0.00125_dp, eps_1 = 0.0005_dp
    real(dp), parameter :: I_v_1 = -8.0_dp*eps_1*sqrt(2*eps_0)/(B_0**2*N_tor)
    real(dp), parameter :: phi_0 = pi
    type(anti_sigma_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-4
    real(dp), parameter :: reltol = 1e-2, abstol = 1e-15
    real(dp), parameter :: stor = 0.5_dp
    integer, parameter :: n_fieldlines = 50

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.0_dp ! formula I_v_1 for small iota
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    type(fieldline_modes_t) :: fieldline_modes

    integer :: current
    real(dp) :: B_mod
    real(dp), dimension(:), allocatable :: zeros

    call field%anti_sigma_field_init(N_tor, B_0, eps_0, eps_1)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, theta_0, iota, field, M_pol, N_tor, &
                                  phi_tol)

    do current = 1, n_fieldlines
        if (not_same(phi_0, modulo(fieldlines(current)%phi_0, 2.0_dp*pi), &
                     abstol_in=2.0_dp*phi_tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_fourier_transform_over_label failed: phi_0"
            print *, "found: ", fieldlines(current)%phi_0
            print *, "expected: ", phi_0
            error stop
        end if
    end do

    call fourier_transform_over_label(field, fieldlines, fieldline_modes)

    if (not_same(I_v_1, fieldline_modes%radial_drift%sin_coeffs(2), &
                 reltol_in=reltol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_fourier_transform_over_label failed: 1st radial drift sin mode"
        print *, "found: ", fieldline_modes%radial_drift%sin_coeffs(2)
        print *, "expected: ", I_v_1
        print *, "ratio: ", fieldline_modes%radial_drift%sin_coeffs(2)/I_v_1
        stop
    end if

    allocate (zeros(size(fieldline_modes%radial_drift%cos_coeffs)))
    zeros = 0.0_dp
    if (not_same(zeros, fieldline_modes%radial_drift%cos_coeffs, &
                 abstol_in=abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_fourier_transform_over_label failed: radial drift cos modes"
        print *, "found: ", fieldline_modes%radial_drift%cos_coeffs
        print *, "expected: all ", zeros(1)
        stop
    end if
    deallocate (zeros)

end program test_fourier_transform_over_label
