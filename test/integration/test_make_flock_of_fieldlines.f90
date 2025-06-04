program test_make_flock_of_fieldlines
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use utils, only: linspace
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: bc_filename = "input/single_mode_m_2_n_minus4.bc"
    real(dp), parameter :: M_pol = 2.0_dp, N_tor = -4.0_dp
    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 1e-4
    real(dp), parameter :: abstol = 2*phi_tol
    real(dp), parameter :: stor = 0.5_dp
    integer, parameter :: n_fieldlines = 10

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), parameter :: iota = -1.0_dp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    integer :: current
    real(dp) :: phi_max(2)

    call field%neo_field_init(bc_filename, stor)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)

    call make_flock_of_fieldlines(fieldlines, theta_0, iota, field, M_pol, N_tor, &
                                  phi_tol)

    do current = 1, n_fieldlines
        phi_max = (/-0.5*pi, 0.5*pi/) + fieldlines(current)%phi_0
        if (not_same(phi_max, fieldlines(current)%phi_max, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_make_flock_of_fieldlines failed: phi_max"
            print *, "found: ", fieldlines(current)%phi_max
            print *, "expected: ", phi_max
            error stop
        end if
    end do

end program test_make_flock_of_fieldlines
