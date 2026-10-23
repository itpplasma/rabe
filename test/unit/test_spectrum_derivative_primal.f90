program test_spectrum_derivative_primal
    !! The smooth primal of the derivative mode (exact maxima, fixed
    !! Gauss-Legendre) must reproduce the standard rabe pipeline (bracketed and
    !! nudged maxima, adaptive QAG at rel tol 1e-6) on the same exact spectral
    !! field, and it must be converged in the number of Gauss points.
    use constants, only: dp
    use spectral_surface, only: spectral_surface_t, spectral_surface_init
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, n_outputs, &
                                   pack_parameters, find_structure, evaluate, &
                                   off_iota
    use synthetic_spectra, only: synthetic_qh, synthetic_qa, synthetic_qi

    implicit none

    type(surface_spectrum_t) :: spec
    type(derivative_settings_t) :: settings
    integer :: n_failed

    n_failed = 0
    call synthetic_qh(spec, settings)
    call compare("QH", spec, settings, n_failed)
    call synthetic_qa(spec, settings)
    call compare("QA", spec, settings, n_failed)
    call synthetic_qi(spec, settings)
    call compare("QI", spec, settings, n_failed)
    if (n_failed > 0) then
        print *, "FAILED checks: ", n_failed
        error stop 1
    end if
    print *, "test_spectrum_derivative_primal: all checks passed"

contains

    subroutine compare(name, spec, settings, n_failed)
        use fieldline_mod, only: flock_of_fieldlines_t
        use make_fieldline, only: make_flock_of_fieldlines
        use coefficients, only: calc_nu_star_crit, &
                                calc_finite_boundary_layer_correction, &
                                calc_gradient_scaling_factor_r_eff, &
                                calc_offset_coefficients
        character(len=*), intent(in) :: name
        type(surface_spectrum_t), intent(in) :: spec
        type(derivative_settings_t), intent(in) :: settings
        integer, intent(inout) :: n_failed

        ! The standard pipeline nudges maxima by the bracket error and
        ! integrates to rel tol 1e-6, so agreement is limited to about that
        ! level on B-dominated quantities and somewhat worse on the small
        ! symmetry-breaking deviations.
        real(dp), parameter :: tol(n_outputs) = [1e-4_dp, 1e-4_dp, 1e-6_dp, &
                                                 1e-6_dp, 1e-8_dp]
        real(dp), parameter :: tol_gauss = 1e-10_dp

        type(spectral_surface_t) :: field
        type(flock_of_fieldlines_t) :: flock
        type(derivative_structure_t) :: structure
        type(derivative_settings_t) :: fine, lattice
        real(dp), allocatable :: params(:)
        real(dp) :: standard(n_outputs), smooth(n_outputs), smooth_fine(n_outputs)
        real(dp) :: dr, rel(n_outputs), rel_gauss(n_outputs)
        integer :: split, K

        ! The standard pipeline needs a minimum at the origin, so it runs on
        ! the field with the origin chosen by find_structure (QI: half a
        ! period away from the booz_xform origin).
        call find_structure(spec, settings, structure)
        call spectral_surface_init(field, spec%xm, spec%xn, spec%bmnc, &
                                   spec%nabla_s_mnc, spec%I, spec%G)
        field%theta_origin = structure%theta_origin
        field%phi_origin = structure%phi_origin
        call make_flock_of_fieldlines(flock, settings%max_n_fieldlines, spec%iota, &
                                      field, settings%M_pol, settings%N_tor, &
                                      spec%nfp, split)
        dr = calc_gradient_scaling_factor_r_eff(flock, spec%psi_edge, &
                                                settings%sign_sqrtg)
        call calc_offset_coefficients(flock, spec%R, dr, standard(1), standard(2))
        standard(3) = calc_finite_boundary_layer_correction(flock, field, spec%R, dr)
        standard(4) = calc_nu_star_crit(flock, spec%R)
        standard(5) = flock%eta_b

        call pack_parameters(spec, params)
        K = size(spec%xm)
        params(2*K + off_iota) = structure%iota_line
        lattice = settings
        lattice%continuum_extrema = .false.
        call evaluate(params, spec, lattice, structure, smooth)
        fine = lattice
        fine%n_gauss = 2*settings%n_gauss
        call evaluate(params, spec, fine, structure, smooth_fine)
        call continuum_shift(name, spec, settings, n_failed)

        rel = abs(smooth - standard)/abs(standard)
        rel_gauss = abs(smooth_fine - smooth)/abs(smooth_fine)
        print "(A,A)", name, ": outputs Lambda_A Lambda_B Lambda_S nu_star_crit eta_b"
        print "(A,5ES14.6)", "  standard pipeline   ", standard
        print "(A,5ES14.6)", "  derivative (lattice)", smooth
        print "(A,5ES10.2)", "  rel difference      ", rel
        print "(A,5ES10.2)", "  rel change n_gauss x2 ", rel_gauss
        if (any(rel > tol)) then
            print *, name, ": derivative primal deviates from the standard pipeline"
            n_failed = n_failed + 1
        end if
        if (any(rel_gauss > tol_gauss)) then
            print *, name, ": Gauss-Legendre rule not converged"
            n_failed = n_failed + 1
        end if
    end subroutine compare

    !> The continuum extrema (global maximum and ridge saddle of B) replace the
    !! discrete max and min over labels; the difference is a label-sampling
    !! error that must shrink when the number of field lines grows.
    subroutine continuum_shift(name, spec, settings, n_failed)
        character(len=*), intent(in) :: name
        type(surface_spectrum_t), intent(in) :: spec
        type(derivative_settings_t), intent(in) :: settings
        integer, intent(inout) :: n_failed

        type(derivative_settings_t) :: cases(2), lattice
        type(derivative_structure_t) :: structure
        real(dp), allocatable :: params(:)
        real(dp) :: v_cont(n_outputs), v_lat(n_outputs), shift(n_outputs, 2)
        integer :: i, K

        K = size(spec%xm)
        cases = settings
        cases(2)%max_n_fieldlines = 4*settings%max_n_fieldlines
        do i = 1, 2
            call find_structure(spec, cases(i), structure)
            call pack_parameters(spec, params)
            params(2*K + off_iota) = structure%iota_line
            lattice = cases(i)
            lattice%continuum_extrema = .false.
            call evaluate(params, spec, cases(i), structure, v_cont)
            call evaluate(params, spec, lattice, structure, v_lat)
            shift(:, i) = abs(v_cont - v_lat)/abs(v_lat)
            print "(A,A,I5,A,5ES10.2)", name, ": lines ", structure%n_lines, &
                "  rel continuum - lattice ", shift(:, i)
        end do
        ! eta_b and nu_star_crit converge with the label spacing
        if (shift(5, 2) > shift(5, 1) .or. shift(4, 2) > shift(4, 1)) then
            print *, name, ": continuum shift does not decrease with more lines"
            n_failed = n_failed + 1
        end if
    end subroutine continuum_shift

end program test_spectrum_derivative_primal
