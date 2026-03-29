module make_fieldline
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

    integer, parameter :: maximum_possible_n = 10
    type :: maxima_t
        integer :: n
        real(dp), dimension(maximum_possible_n) :: phi
        real(dp), dimension(maximum_possible_n) :: B
    end type maxima_t

contains

    subroutine make_flock_of_fieldlines(fieldlines, xi_0, iota, &
                                        field, M_pol, N_tor, nfp, &
                                        phi_tol, &
                                        err_flag)
        use fieldline_integrals, only: calc_fieldline_integrals
        use fieldline_labels, only: calc_iota_p
        use fieldline_labels, only: suspect_omnigenous_origin_not_minimum
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        real(dp), dimension(:), intent(in) :: xi_0
        real(dp), intent(in) :: iota
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), intent(in) :: phi_tol
        integer, intent(out), optional :: err_flag

        real(dp) :: interval(2)
        type(maxima_t) :: maxima
        logical :: more_than_two_maxima
        real(dp) :: I_ref
        integer :: n_fieldlines
        integer :: current

        if (present(err_flag)) then
            err_flag = 0
        end if

        call check_if_valid_input(M_pol, N_tor, nfp, iota)

        n_fieldlines = size(fieldlines)

        fieldlines%xi_0 = xi_0
        fieldlines%iota = iota
        fieldlines%M_pol = M_pol
        fieldlines%N_tor = N_tor
        fieldlines%nfp = nfp

        if (suspect_omnigenous_origin_not_minimum(field, M_pol, N_tor, phi_tol)) then
            print *, "error: The origin of the IDEAL omnigenous configuration"
            print *, "(theta=phi=0) must be a global and local minimum!"
            print *, "Origin of provided field suggests that this is not the case!"
            error stop
        end if
        !> if the origin of the ideal omnigenous field is a minimum (above)
        !> we can put the labels along the chi = 0 line
        fieldlines%theta_0 = N_tor*fieldlines%xi_0/nfp
        fieldlines%phi_0 = M_pol*fieldlines%xi_0/nfp

        fieldlines%iota_p = calc_iota_p(iota, M_pol, N_tor, nfp)

        more_than_two_maxima = .false.
        do current = 1, n_fieldlines
            interval = [-1.5_dp*pi, 1.5_dp*pi]/abs(N_tor - iota*M_pol) + &
                       fieldlines(current)%phi_0
            call find_maxima_along_fieldline(field, fieldlines(current), &
                                             interval, maxima, phi_tol)
            if (maxima%n < 2) then
                print *, "---------------------------------------------------------"
                print *, "---------------------------------------------------------"
                print *, "---------------------------------------------------------"
                print *, "Found less than two maxima in provided interval!"
                print *, "theta_0: ", fieldlines(current)%theta_0
                print *, "phi_0: ", fieldlines(current)%phi_0
                print *, "interval: ", interval
                print *, "phi_max: ", maxima%phi
                print *, "B_max: ", maxima%B
                print *, "---------------------------------------------------------"
                print *, "---------------------------------------------------------"
                print *, "---------------------------------------------------------"
                error stop
            elseif (maxima%n > 2) then
                fieldlines(current)%phi_max = get_biggest_maxima_on_each_side(maxima, &
                                                              fieldlines(current)%phi_0)
                more_than_two_maxima = .true.
            else
                fieldlines(current)%phi_max = maxima%phi(1:2)
            end if

            ! To ensure that the there are no maxima in between found phi_max
            ! we move phi_max inside the well by the maximal potential error
            call nudge_maxima_inward(field, fieldlines(current), phi_tol)
        end do

        if (more_than_two_maxima) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "warning in make_flock_of_fieldlines: "
            print *, "The provided field violates omnigeneity too strongly!"
            print *, "-> Found more than two local maxima per period", &
                " for at least one fieldline!"
            print *, "Calculation done with biggest maximum in each half period!"
            print *, "Final result for bootstrap deviation can not be trusted!"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            if (present(err_flag)) then
                err_flag = 1
            end if
        end if

        fieldlines%eta_b = (1.0_dp)/get_global_B_max(fieldlines)

        do current = 1, n_fieldlines
            call calc_fieldline_integrals(field, fieldlines(current))
        end do

        fieldlines%delta_eta = 1.0_dp/fieldlines(:)%B_max(1) - fieldlines(:)%eta_b
        ! I_ref can be chosen to be any I
        ! (I_ref/I_j)**0.5 - 1 = (max(I_j)/I_j)**0.5 -1 =
        ! ((I+delta)/(I+delta_j))**0.5 -1 ~ 0.5*(delta/I - delta_j/I)
        ! and the result in linear order only differs by a constant delta/I
        ! which does not enter the offset formula
        ! We choose I_ref so that the average of delta_aspect is zero
        fieldlines%I_ref = ( &
                           n_fieldlines/ &
                       (sum(1.0_dp/sqrt(fieldlines%integral_lambda_b_over_B_squared))) &
                           )**2.0_dp
        fieldlines%delta_aspect_ratio = sqrt( &
                                        fieldlines%I_ref/ &
                                        fieldlines%integral_lambda_b_over_B_squared &
                                        ) - 1.0_dp

    end subroutine make_flock_of_fieldlines

    subroutine check_if_valid_input(M_pol, N_tor, nfp, iota)
        use utils, only: not_same
        real(dp), intent(in) :: M_pol, N_tor, nfp, iota

        real(dp), parameter :: tol = 1e-15
        logical :: is_valid

        is_valid = .true.

        if (is_not_integer(M_pol, tol)) then
            print *, "M_pol must be integer"
            is_valid = .false.
        end if
        if (is_not_integer(N_tor, tol)) then
            print *, "N_tor must be integer"
            is_valid = .false.
        end if
        if (is_not_integer(nfp, tol)) then
            print *, "nfp must be integer"
            is_valid = .false.
        end if
        if (nint(nfp) <= 0) then
            print *, "nfp must be positiv"
            is_valid = .false.
        end if
        if (nint(N_tor) /= 0) then
            if (not_same(N_tor, nfp, reltol_in=tol, abstol_in=0.0_dp)) then
                print *, "nonzero N_tor must be equal nfp"
                is_valid = .false.
            end if
        else
            if (nint(M_pol) /= 1) then
                is_valid = .false.
                print *, "M_pol must be 1 if N_tor=0"
            end if
        end if

        if (abs(M_pol*iota - N_tor) < tol) then
            print *, "Error: (M_pol*iota - N_tor) must not be (close) zero."
            print *, "abs(M_pol*iota - N_tor) = ", abs(M_pol*iota - N_tor)
            is_valid = .false.
        end if

        if (.not. is_valid) then
            print *, "Error: not valid input:"
            print *, "M_pol: ", M_pol
            print *, "N_tor: ", N_tor
            print *, "nfp: ", nfp
            print *, "iota: ", iota
            error stop
        end if

    end subroutine check_if_valid_input

    logical function is_not_integer(x, tol)
        real(dp), intent(in) :: x, tol

        is_not_integer = abs(x - nint(x)) > tol
    end function is_not_integer

    subroutine find_maxima_along_fieldline(field, &
                                           fieldline, &
                                           interval, &
                                           maxima, &
                                           phi_tol)
        use find_extrema, only: find_local_maxima
        use, intrinsic :: ieee_arithmetic

        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline
        real(dp), intent(in) :: interval(2)
        real(dp), intent(in), optional :: phi_tol
        type(maxima_t) :: maxima

        call find_local_maxima(B_mod_along_fieldline, interval, &
                               maxima%phi, phi_tol)

        maxima%n = count(.not. ieee_is_nan(maxima%phi))
        call B_mod_along_fieldline(maxima%phi(1:maxima%n), maxima%B(1:maxima%n))

    contains
        subroutine B_mod_along_fieldline(phi, B_mod)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(phi)) :: theta
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_fieldline
    end subroutine find_maxima_along_fieldline

    function get_biggest_maxima_on_each_side(maxima, phi_0) result(phi_max)
        type(maxima_t), intent(in) :: maxima
        real(dp), intent(in) :: phi_0
        real(dp), dimension(2) :: phi_max

        integer :: of_biggest_B

        of_biggest_B = maxloc(maxima%B, mask=(maxima%phi < phi_0), dim=1)
        phi_max(1) = maxima%phi(of_biggest_B)
        of_biggest_B = maxloc(maxima%B, mask=(maxima%phi > phi_0), dim=1)
        phi_max(2) = maxima%phi(of_biggest_B)
    end function get_biggest_maxima_on_each_side

    subroutine nudge_maxima_inward(field, fieldline, phi_tol)
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline
        real(dp), intent(in) :: phi_tol

        fieldline%phi_max(1) = fieldline%phi_max(1) + phi_tol
        fieldline%phi_max(2) = fieldline%phi_max(2) - phi_tol

        call B_mod_along_fieldline(fieldline%phi_max, fieldline%B_max)

    contains
        subroutine B_mod_along_fieldline(phi, B_mod)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(phi)) :: theta
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_fieldline
    end subroutine nudge_maxima_inward

    function get_global_B_max(fieldlines) result(global_B_max)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: global_B_max

        integer :: current
        real(dp) :: B_max_1, B_max_2

        B_max_1 = maxval(fieldlines(:)%B_max(1))
        B_max_2 = maxval(fieldlines(:)%B_max(2))
        global_B_max = max(B_max_1, B_max_2)
    end function get_global_B_max

end module make_fieldline
