program rabe
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use find_extrema, only: find_local_maxima
    use fieldline_labels, only: guess_alpha_at_minimum

    implicit none

    character(len=*), parameter :: bc_filename = "test/integration/input/"// &
                                   "quasi_helical.bc"
    real(dp), parameter :: M = 1.0_dp
    real(dp), parameter :: nfp = 4.0_dp, scan_n_periods = 2.0_dp
    real(dp), dimension(2), parameter :: interval = (/0.0_dp, &
                                                      2.0_dp*pi/nfp*scan_n_periods/)

    type(neo_field_t) :: field
    real(dp) :: bmod, sqrtg, dB_dx(3)
    real(dp) :: found_phi_max(2)
    real(dp) :: iota, alpha_at_min

    call field%neo_field_init(bc_filename, stor=0.5_dp)
    iota = field%iota
    print *, iota
    call find_local_maxima(B_mod_along_fieldline, interval, found_phi_max)
    print *, found_phi_max
    call guess_alpha_at_minimum(field, alpha_at_min, M)
    print *, alpha_at_min

contains

    subroutine B_mod_along_fieldline(phi, B_mod)
        real(dp), dimension(:), intent(in) :: phi
        real(dp), dimension(:), intent(out) :: B_mod

        real(dp), dimension(size(phi, 1)) :: theta
        integer :: idx

        theta = phi*iota

        do idx = 1, size(phi, 1)
            call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
        end do
    end subroutine B_mod_along_fieldline

end program rabe
