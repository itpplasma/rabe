program rabe
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use fieldline, only: guess_alpha_at_minimum, find_maxima_along_fieldline

    implicit none

    character(len=*), parameter :: bc_filename = "test/integration/input/"// &
                                   "quasi_helical.bc"
    real(dp), parameter :: M = 1.0_dp, theta_0 = 0.0_dp
    real(dp), dimension(2), parameter :: interval = (/0.0_dp, 2.0_dp*pi/)

    type(neo_field_t) :: field
    real(dp) :: bmod, sqrtg, dB_dx(3)
    real(dp) :: found_phi_max(2)
    real(dp) :: iota, alpha_at_min

    call field%neo_field_init(bc_filename, stor=0.5_dp)
    iota = field%iota
    print *, iota
    call find_maxima_along_fieldline(field, iota, theta_0, interval, found_phi_max)
    print *, found_phi_max
    call guess_alpha_at_minimum(field, alpha_at_min, M)
    print *, alpha_at_min

end program rabe
