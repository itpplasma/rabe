program rabe
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use fieldline_mod, only: guess_alpha_over_M_at_minimum, find_maxima_along_fieldline

    implicit none

    character(len=*), parameter :: bc_filename = "test/integration/input/"// &
                                   "quasi_helical.bc"
    type(neo_field_t) :: field

    call field%neo_field_init(bc_filename, stor=0.5_dp)

end program rabe
