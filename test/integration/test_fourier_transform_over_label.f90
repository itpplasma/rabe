program test_fourier_transform_over_label
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_mod, only: make_flock_of_fieldlines
    use fieldline_integrals, only: fourier_transform_over_label
    use utils, only: linspace
    use utils, only: is_same

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = 0.125_dp, eps_1 = 0.0
    real(dp), parameter :: phi_0 = pi
    type(anti_sigma_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-4
    real(dp), parameter :: abstol = 2*phi_tol
    real(dp), parameter :: stor = 0.5_dp
    integer, parameter :: n_fieldlines = 10

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), parameter :: iota = 0.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    integer :: current
    real(dp) :: B_mod

    call field%anti_sigma_field_init(N_tor, B_0, eps_0, eps_1)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)

    call make_flock_of_fieldlines(fieldlines, theta_0, iota, field, M_pol, N_tor, &
                                  phi_tol)

    do current = 1, n_fieldlines
        if (is_same(phi_0, modulo(fieldlines(current)%phi_0, 2.0_dp*pi), &
                    abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_fourier_transform_over_label failed: phi_0"
            print *, "found: ", fieldlines(current)%phi_0
            print *, "expected: ", phi_0
            error stop
        end if
    end do

    call fourier_transform_over_label(field, fieldlines)

end program test_fourier_transform_over_label
