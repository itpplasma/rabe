module make_fieldline
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none

    type :: maxima_t
        integer :: n
        real(dp), dimension(:), allocatable :: phi
        real(dp), dimension(:), allocatable :: B
        real(dp), dimension(:), allocatable :: B_error
        real(dp), dimension(:), allocatable :: phi_error
    end type maxima_t

contains

    subroutine make_flock_of_fieldlines(fieldlines, xi_0, iota, &
                                        field, M_pol, N_tor, nfp, &
                                        split_maxima)
        use fieldline_integrals, only: calc_fieldline_integrals
        use fieldline_labels, only: calc_iota_p
        use fieldline_labels, only: suspect_omnigenous_origin_not_minimum
        use constants, only: machine_eps
        type(fieldline_t), dimension(:), intent(inout) :: fieldlines
        real(dp), dimension(:), intent(in) :: xi_0
        real(dp), intent(in) :: iota
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp
        integer, intent(out), optional :: split_maxima

        real(dp) :: interval(2)
        type(maxima_t) :: maxima
        logical :: more_than_two_maxima
        real(dp) :: I_ref
        integer :: n_fieldlines
        integer :: current
        real(dp) :: symmetry_violation

        if (present(split_maxima)) then
            split_maxima = 0
        end if

        call check_if_valid_input(M_pol, N_tor, nfp, iota)

        n_fieldlines = size(fieldlines)

        fieldlines%xi_0 = xi_0
        fieldlines%iota = iota
        fieldlines%M_pol = M_pol
        fieldlines%N_tor = N_tor
        fieldlines%nfp = nfp

        symmetry_violation = estimate_symmetry_violation(field, iota, nfp)
        if (symmetry_violation > field%rel_accuracy_B()) then
            print *, "error: provided field violates stellarator symmetry too strongly!"
            print *, "symmetry violation (max|B(theta, phi) - B(-theta, -phi)|/B): ", &
                symmetry_violation
            error stop
        end if

        if (suspect_omnigenous_origin_not_minimum(field, M_pol, N_tor, &
                                                  symmetry_violation)) then
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
                                             interval, maxima)
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
                call pick_maximum_on_each_side(maxima, &
                                               fieldlines(current)%phi_0, &
                                               symmetry_violation, &
                                               fieldlines(current)%phi_max, &
                                               fieldlines(current)%phi_max_error)
                more_than_two_maxima = .true.
            else
                fieldlines(current)%phi_max = maxima%phi(1:2)
                fieldlines(current)%phi_max_error = maxima%phi_error(1:2)
            end if

            ! To ensure that there are no maxima in between found phi_max
            ! we move phi_max inside the well by the maximal potential error
            call nudge_maxima_inward(field, fieldlines(current))
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
            if (present(split_maxima)) then
                split_maxima = 1
            end if
        end if
        !> If eta_b is chosen exactly to be 1/B_max, one relies that
        !> B_max/B_max = 1 so that the pitchparameter does not become negative on
        !> the field line hitting the global maximum. However,
        !> B_max/Bmax = 1 +- ULP (~2.2e-16)
        !> due to floating point arithmetic. To not be dependend on machine/compiler
        !> specifics, we add a small buffer that does not change the physics.
        fieldlines%eta_b = (1.0_dp - 2.0_dp*machine_eps) &
                           /get_global_B_max(fieldlines)

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

    function estimate_symmetry_violation(field, iota, nfp) result(symmetry_violation)
        use utils, only: linspace
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: nfp
        real(dp) :: symmetry_violation

        integer, parameter :: n_points = 1000
        integer :: idx
        real(dp), dimension(n_points) :: theta, phi, B_ref, B_sym

        call linspace(0.0_dp, 2.0_dp*pi/nfp, n_points, phi)
        theta = iota*phi
        do idx = 1, n_points
            call field%compute_B_mod(theta(idx), phi(idx), B_ref(idx))
        end do
        phi = -phi
        theta = -theta
        do idx = 1, n_points
            call field%compute_B_mod(theta(idx), phi(idx), B_sym(idx))
        end do
        symmetry_violation = maxval(abs(B_ref - B_sym)/B_ref)
    end function estimate_symmetry_violation

    subroutine find_maxima_along_fieldline(field, &
                                           fieldline, &
                                           interval, &
                                           maxima)
        use find_extrema, only: find_local_maxima

        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline
        real(dp), intent(in) :: interval(2)
        type(maxima_t), intent(out) :: maxima

        call find_local_maxima(B_mod_along_fieldline, interval, &
                               maxima%phi, maxima%phi_error)

        maxima%n = size(maxima%phi)
        allocate (maxima%B(maxima%n), maxima%B_error(maxima%n))
        call B_mod_along_fieldline(maxima%phi, maxima%B)
        maxima%B_error = 0.0_dp
        call dB_dphi_along_fieldline(maxima%phi, maxima%B_error)
        maxima%B_error(1:maxima%n) = abs(maxima%B_error)*maxima%phi_error

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

        subroutine dB_dphi_along_fieldline(phi, dB_dphi)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: dB_dphi

            real(dp), dimension(size(phi)) :: theta
            real(dp) :: B_mod, dB_dx(3)
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_and_dB_dx(theta(idx), phi(idx), &
                                               B_mod, dB_dx)
                dB_dphi(idx) = dB_dx(3) + fieldline%iota*dB_dx(2)
            end do
        end subroutine dB_dphi_along_fieldline
    end subroutine find_maxima_along_fieldline

    !> Pick the biggest maximum on each side of phi_0, with ties broken by
    !> proximity to phi_0 to respect stellarator symmetry.
    subroutine pick_maximum_on_each_side(maxima, phi_0, symmetry_violation, &
                                         phi_max, phi_max_error)
        type(maxima_t), intent(in) :: maxima
        real(dp), intent(in) :: phi_0
        real(dp), intent(in) :: symmetry_violation
        real(dp), dimension(2), intent(out) :: phi_max
        real(dp), dimension(2), intent(out) :: phi_max_error

        integer :: idx

        idx = pick_maximum(maxima%phi, maxima%B, phi_0, maxima%B_error, &
                           symmetry_violation, mask=maxima%phi < phi_0)
        phi_max(1) = maxima%phi(idx)
        phi_max_error(1) = maxima%phi_error(idx)

        idx = pick_maximum(maxima%phi, maxima%B, phi_0, maxima%B_error, &
                           symmetry_violation, mask=maxima%phi > phi_0)
        phi_max(2) = maxima%phi(idx)
        phi_max_error(2) = maxima%phi_error(idx)
    end subroutine pick_maximum_on_each_side

    !> Pick the biggest maximum from the masked set. If multiple maxima are
    !> equivalent to the biggest, i.e. due to
    !> (a) their difference being within their error bounds or
    !> (b) them being symmetric mirrors within the symmetry violation -
    !> the one closest to phi_0 is chosen.
    function pick_maximum(phi, B, phi_0, error, symmetry_violation, mask) result(idx)
        real(dp), dimension(:), intent(in) :: phi, B
        real(dp), intent(in) :: phi_0
        real(dp), dimension(:), intent(in) :: error
        real(dp), intent(in) :: symmetry_violation
        logical, dimension(:), intent(in) :: mask
        integer :: idx

        real(dp), dimension(size(phi)) :: tol
        real(dp) :: biggest_B, error_of_biggest
        logical, dimension(size(phi)) :: equal_to_biggest

        idx = maxloc(B, mask=mask, dim=1)
        biggest_B = B(idx)
        error_of_biggest = error(idx)
        tol = error_of_biggest + error + 2.0_dp*symmetry_violation*B
        equal_to_biggest = mask .and. (abs(B - biggest_B) <= tol)
        idx = minloc(abs(phi - phi_0), mask=equal_to_biggest, dim=1)
    end function pick_maximum

    subroutine nudge_maxima_inward(field, fieldline)
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(inout) :: fieldline

        fieldline%phi_max(1) = fieldline%phi_max(1) + fieldline%phi_max_error(1)
        fieldline%phi_max(2) = fieldline%phi_max(2) - fieldline%phi_max_error(2)

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
