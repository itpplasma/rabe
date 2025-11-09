module make_fieldline
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

contains

    subroutine make_flock_of_fieldlines(fieldlines, xi_0, iota, &
                                        field, M_pol, N_tor, nfp, phi_tol)
        use fieldline_integrals, only: calc_fieldline_integrals
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        real(dp), dimension(:), intent(in) :: xi_0
        real(dp), intent(in) :: iota
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), intent(in), optional :: phi_tol

        real(dp) :: interval(2)
        real(dp) :: normalization
        real(dp) :: I_ref
        integer :: n_fieldlines
        integer :: current
        logical :: more_than_2_maxima, too_strong_violation

        call check_if_valid_input(M_pol, N_tor, nfp)

        n_fieldlines = size(fieldlines)

        fieldlines%xi_0 = xi_0
        fieldlines%iota = iota
        fieldlines%M_pol = M_pol
        fieldlines%N_tor = N_tor
        fieldlines%nfp = nfp

        call set_fieldline_labels_along_chi_min(field, M_pol, N_tor, nfp, fieldlines, &
                                                phi_tol)

        normalization = N_tor**2.0_dp + M_pol**2.0_dp
        fieldlines%iota_p = sign(pi, iota*M_pol - N_tor)/normalization* &
                            (M_pol + &
                             nfp*(N_tor*iota + M_pol)/(iota*M_pol - N_tor))

        too_strong_violation = .false.
        do current = 1, n_fieldlines
            interval = [-1.5_dp*pi, 1.5_dp*pi]/abs(N_tor - iota*M_pol) + &
                       fieldlines(current)%phi_0
            call find_maxima_along_fieldline(field, fieldlines(current), &
                                             interval, phi_tol, more_than_2_maxima)
            if (more_than_2_maxima) too_strong_violation = .true.
        end do

        if (too_strong_violation) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "warning in make_flock_of_fieldlines: "
            print *, "The provided field violates omnigeneity too strongly!"
            print *, "-> Found more than one local maximum per half period", &
                " for at least one fieldline!"
            print *, "Calculation done with biggest maximum in each half period!"
            print *, "Final result for bootstrap deviation can not be trusted!"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
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
                                        ) - 1

    end subroutine make_flock_of_fieldlines

    subroutine check_if_valid_input(M_pol, N_tor, nfp)
        use utils, only: not_same
        real(dp), intent(in) :: M_pol, N_tor, nfp

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

        if (.not. is_valid) then
            print *, "Error: not valid input:"
            print *, "M_pol: ", M_pol
            print *, "N_tor: ", N_tor
            print *, "nfp: ", nfp
            error stop
        end if

    end subroutine check_if_valid_input

    logical function is_not_integer(x, tol)
        real(dp), intent(in) :: x, tol

        is_not_integer = abs(x - nint(x)) > tol
    end function is_not_integer

    subroutine set_fieldline_labels_along_chi_min(field, M_pol, N_tor, nfp, &
                                                  fieldlines, phi_tol)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), intent(in), optional :: phi_tol
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines

        real(dp) :: chi_min, tol

        call guess_chi_min(field, chi_min, N_tor, M_pol, phi_tol)

        if (present(phi_tol)) then
            tol = phi_tol*3.0_dp
        else
            tol = 3.0_dp*1e-2
        end if
        if (not_multiple_of_2pi(chi_min, tol)) then
            print *, "error: found chi_min is not multiple of 2pi"
            print *, "chi_min: ", chi_min/pi, "[pi]"
            print *, "The minima contour of the ideal omnigenous configuration"
            print *, "must pass through (theta=0,phi=0)!"
            error stop
        end if

        fieldlines%theta_0 = N_tor*fieldlines%xi_0/nfp
        fieldlines%phi_0 = M_pol*fieldlines%xi_0/nfp
    end subroutine set_fieldline_labels_along_chi_min

    subroutine guess_chi_min(field, chi_min, N_tor, M_pol, tol)
        use find_extrema, only: find_local_minima

        class(field_t), intent(in) :: field
        real(dp), intent(out) :: chi_min
        real(dp), intent(in) :: N_tor, M_pol
        real(dp), intent(in), optional :: tol

        ! chi = M*theta - N*phi
        ! as f~f(chi) = sum c_j*cos(j*chi) with 1<j<j_max
        ! periodic at least after 2pi -> minimum must be in e.g [-pi, 2pi]
        real(dp), dimension(2), parameter :: interval = [0.0_dp, 3.0_dp*pi]
        real(dp) :: location(1)

        ! If f(theta,phi) approx f(chi = M*theta - N*phi) one can estiamte
        ! f(chi/N) by choosing 1 specific theta-phi combination for that
        ! chi value e.g.
        ! - phi=-chi/N and theta=0 or
        ! - phi=0 and theta=chi/M

        if (nint(N_tor) /= 0) then
            call find_local_minima(B_mod_along_phi_axis, interval, location, tol)
        elseif (nint(M_pol) /= 0) then
            call find_local_minima(B_mod_along_theta_axis, interval, location, tol)
        else
            print *, "error in guess_chi_min: M_pol=N_tor=0"
            print *, "M_pol and N_tor must not be both zero!"
            error stop
        end if

        chi_min = location(1)

    contains

        subroutine B_mod_along_phi_axis(chi, B_mod)
            real(dp), dimension(:), intent(in) :: chi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(chi, 1)) :: phi, theta
            integer :: idx

            phi = -chi/N_tor
            theta = 0.0_dp

            do idx = 1, size(chi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_phi_axis

        subroutine B_mod_along_theta_axis(chi, B_mod)
            real(dp), dimension(:), intent(in) :: chi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(chi, 1)) :: phi, theta
            integer :: idx

            theta = chi/M_pol
            phi = 0.0_dp

            do idx = 1, size(chi, 1)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_theta_axis

    end subroutine guess_chi_min

    subroutine find_maxima_along_fieldline(field, &
                                           fieldline, &
                                           interval, &
                                           phi_tol, &
                                           more_than_2_maxima)
        use find_extrema, only: find_local_maxima
        use, intrinsic :: ieee_arithmetic

        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline
        real(dp), intent(in) :: interval(2)
        real(dp), intent(in), optional :: phi_tol
        logical, intent(out), optional :: more_than_2_maxima

        integer, parameter :: n_max = 10
        real(dp), dimension(n_max) :: phi_max, B_max
        integer :: found_maxima
        integer :: of_biggest_B

        if (present(more_than_2_maxima)) more_than_2_maxima = .false.

        call find_local_maxima(B_mod_along_fieldline, interval, &
                               phi_max, phi_tol)

        found_maxima = count(.not. ieee_is_nan(phi_max))

        if (found_maxima < 2) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "error in find_maxima_along_fieldline: "
            print *, "Found less than two maxima in provided interval!"
            print *, "theta_0: ", fieldline%theta_0
            print *, "phi_0: ", fieldline%phi_0
            print *, "interval: ", interval
            print *, "found maxima: ", phi_max(1:found_maxima)
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            error stop
        elseif (found_maxima > 2) then
            if (present(more_than_2_maxima)) more_than_2_maxima = .true.
            call B_mod_along_fieldline(phi_max(1:found_maxima), B_max(1:found_maxima))
            of_biggest_B = maxloc(B_max, mask=(phi_max < fieldline%phi_0), dim=1)
            fieldline%phi_max(1) = phi_max(of_biggest_B)
            of_biggest_B = maxloc(B_max, mask=(phi_max > fieldline%phi_0), dim=1)
            fieldline%phi_max(2) = phi_max(of_biggest_B)
        else
            fieldline%phi_max = phi_max(1:found_maxima)
        end if

        ! To ensure that the there are no maxima in between found phi_max
        ! we move phi_max inside the well by the maximal potential error
        if (present(phi_tol)) then
            fieldline%phi_max(1) = fieldline%phi_max(1) + phi_tol
            fieldline%phi_max(2) = fieldline%phi_max(2) - phi_tol
        end if

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
    end subroutine find_maxima_along_fieldline

    function get_global_B_max(fieldlines) result(global_B_max)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: global_B_max

        integer :: current
        real(dp) :: B_max_1, B_max_2

        B_max_1 = maxval(fieldlines(:)%B_max(1))
        B_max_2 = maxval(fieldlines(:)%B_max(2))
        global_B_max = max(B_max_1, B_max_2)
    end function get_global_B_max

    function not_multiple_of_2pi(angle, tol)
        real(dp), intent(in) :: angle
        real(dp), intent(in) :: tol
        logical :: not_multiple_of_2pi

        real(dp) :: remainder

        remainder = abs(mod(angle, 2.0_dp*pi))
        not_multiple_of_2pi = remainder > tol .and. abs(remainder - 2.0*pi) > tol
    end function not_multiple_of_2pi

end module make_fieldline
