program rabe
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: set_fieldline_phi_0_to_mode_minimum
    use make_fieldline, only: find_maxima_along_fieldline

    implicit none

    real(dp), parameter :: stor = 0.5_dp
    integer, parameter :: n_fieldlines = 10
    real(dp), parameter :: theta_mode = 1.0_dp, phi_mode = -4.0_dp
    character(len=*), parameter :: bc_filename = "test/integration/input/"// &
                                   "quasi_helical.bc"

    type(neo_field_t) :: field
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines
    real(dp), dimension(n_fieldlines) :: theta_0

    integer :: current
    real(dp) :: interval(2)

    call field%neo_field_init(bc_filename, stor)

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, theta_0)
    fieldlines(:)%theta_0 = theta_0(:)
    fieldlines(:)%iota = field%iota

    call set_fieldline_phi_0_to_mode_minimum(field, theta_mode, phi_mode, &
                                             fieldlines)
    do current = 1, n_fieldlines
        interval = (/0.0_dp, 2*pi/) + fieldlines(current)%phi_0
        call find_maxima_along_fieldline(field, fieldlines(current), &
                                         interval)
    end do

end program rabe
