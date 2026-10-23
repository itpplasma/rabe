program test_derivative_robustness
    !! Regressions for the derivative mode on Landreman-Paul spectra projected
    !! from the reactor-scale low-res references. Usage:
    !! test_derivative_robustness.x mirror|saddle|gauss256
    !!
    !! mirror   LP QA, s = 0.5, 200 lines: lines with two mirror-image maxima on
    !!          one side. A relative spectrum perturbation of 1e-14 must not
    !!          switch the bounding maximum (it moved Lambda_S by 1.4e-4).
    !! saddle   LP QA, s = 0.7, 60 lines: the continuum minimum of the left
    !!          maxima must be found (2D Newton missed it) and nu_star_crit must
    !!          agree with central finite differences.
    !! gauss256 LP QH, s = 0.35, discrete extrema, n_gauss = 256: Gauss nodes
    !!          next to the global maximum give 1 - eta_b*B < 0 by roundoff,
    !!          which must not stop the evaluation. The sign of that roundoff
    !!          depends on the last bits of the spectrum; the surface is one
    !!          where it is negative (-2.4e-15) with the pinned libneo.
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, &
                                   derivative_diagnostics_t, n_outputs, &
                                   n_params, pack_parameters, find_structure, &
                                   evaluate, off_iota, off_G, off_I, i_nu_star_crit
    use spectrum_io, only: spectralize_field

    implicit none

    type(boozer_field_t) :: field
    type(surface_spectrum_t) :: spec
    type(derivative_settings_t) :: settings
    type(derivative_structure_t) :: structure
    type(derivative_diagnostics_t) :: diag
    character(len=16) :: case_name
    real(dp), allocatable :: params(:), jac(:, :), u(:), noise(:), bmnc_0(:)
    real(dp) :: values(n_outputs), values_2(n_outputs), rel_change(n_outputs)
    real(dp) :: s_tor, max_change
    integer :: K, P, n_failed, dir, seed

    call get_command_argument(1, case_name)
    settings%sign_sqrtg = -1
    settings%n_gauss = 64
    n_failed = 0
    select case (trim(case_name))
    case ("mirror")
        s_tor = 0.5_dp
        settings%max_n_fieldlines = 200
        call setup("qa", 16)
        call find_structure(spec, settings, structure)
        print "(A,I4,A,I4)", "  lines ", structure%n_lines, "  split lines ", &
            structure%n_split_lines
        if (structure%n_split_lines == 0) then
            print *, "  no line with more than two maxima: test does not apply"
            n_failed = n_failed + 1
        end if
        call evaluate(params, spec, settings, structure, values)
        max_change = 0.0_dp
        allocate (noise(K))
        bmnc_0 = spec%bmnc
        do seed = 1, 8
            call random_direction(noise, 100 + seed)
            ! structure and evaluation of the perturbed spectrum, as a caller
            ! would do it
            spec%bmnc = bmnc_0*(1.0_dp + 1e-14_dp*noise)
            call find_structure(spec, settings, structure)
            call pack_parameters(spec, params)
            call evaluate(params, spec, settings, structure, values_2)
            rel_change = abs(values_2 - values)/abs(values)
            max_change = max(max_change, maxval(rel_change))
        end do
        print "(A,ES10.2)", "  max relative change under 1e-14 noise ", max_change
        if (max_change > 1e-9_dp) n_failed = n_failed + 1

    case ("saddle")
        s_tor = 0.7_dp
        settings%max_n_fieldlines = 60
        call setup("qa", 8)
        call find_structure(spec, settings, structure)
        allocate (jac(P, n_outputs), u(P))
        call evaluate(params, spec, settings, structure, values, jac, diag)
        print "(A,L2,A,L2,A,ES10.2)", "  continuum extrema ", &
            diag%continuum_extrema, "  crossing ", diag%nu_crit_crossing, &
            "  nu_crit shift ", diag%nu_crit_shift
        if (.not. diag%continuum_extrema) n_failed = n_failed + 1
        do dir = 1, 3
            u = 0.0_dp
            select case (dir)
            case (1)
                call random_direction(u(2:K), 41)
                u(2:K) = u(2:K)*max(abs(spec%bmnc(2:K)), 1e-6_dp*abs(spec%bmnc(1)))
            case (2)
                u(2*K + off_iota) = abs(spec%iota)
            case (3)
                u(2*K + off_G) = spec%G
                u(2*K + off_I) = 0.3_dp*spec%G
            end select
            call fd_check(dir, params, spec, settings, structure, jac, u, n_failed)
        end do

    case ("gauss256")
        s_tor = 0.35_dp
        settings%max_n_fieldlines = 200
        settings%continuum_extrema = .false.
        call setup("qh", 16)
        call find_structure(spec, settings, structure)
        settings%n_gauss = 256
        allocate (jac(P, n_outputs))
        call evaluate(params, spec, settings, structure, values, jac, diag)
        print "(A,ES10.2)", "  min 1 - eta_b*B ", diag%min_lambda_squared
        if (diag%min_lambda_squared >= 0.0_dp) then
            print *, "  no negative roundoff at the nodes: test does not apply"
            n_failed = n_failed + 1
        end if
        if (any(jac /= jac) .or. any(abs(jac) > huge(1.0_dp))) then
            print *, "  Jacobian not finite"
            n_failed = n_failed + 1
        end if
        settings%n_gauss = 128
        call evaluate(params, spec, settings, structure, values_2)
        rel_change = abs(values_2 - values)/abs(values)
        print "(A,7ES10.2)", "  n_gauss 256 vs 128 ", rel_change
        if (maxval(rel_change) > 1e-9_dp) n_failed = n_failed + 1

    case default
        error stop "usage: test_derivative_robustness.x mirror|saddle|gauss256"
    end select

    if (n_failed > 0) then
        print *, "FAILED checks: ", n_failed
        error stop 1
    end if
    print *, "test_derivative_robustness ", trim(case_name), ": all checks passed"

contains

    subroutine setup(equilibrium, n_modes)
        character(len=*), intent(in) :: equilibrium
        integer, intent(in) :: n_modes

        real(dp) :: residual, asymmetry, iota, B_I, B_G

        if (equilibrium == "qh") then
            call field%init_from_vmec( &
                "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc", &
                grid_refinement=6)
            settings%M_pol = -1.0_dp
            settings%N_tor = 4.0_dp
        else
            call field%init_from_vmec( &
                "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc", &
                grid_refinement=6)
            settings%M_pol = 1.0_dp
            settings%N_tor = 0.0_dp
        end if
        call field%fix_to_surface(s_tor)
        call field%get_iota(s_tor, iota)
        call field%get_covariant_components(B_I, B_G)
        spec%iota = iota
        spec%G = B_G
        spec%I = B_I
        spec%psi_edge = field%psi_tor_edge
        spec%R = field%R
        call spectralize_field(field, field%nfp, n_modes, n_modes, spec, &
                               residual, asymmetry)
        call pack_parameters(spec, params)
        K = size(spec%xm)
        P = n_params(K)
    end subroutine setup

    subroutine fd_check(dir, params, template, settings, structure, jac, u, n_failed)
        integer, intent(in) :: dir
        real(dp), intent(in) :: params(:)
        type(surface_spectrum_t), intent(in) :: template
        type(derivative_settings_t), intent(in) :: settings
        type(derivative_structure_t), intent(in) :: structure
        real(dp), intent(in) :: jac(:, :), u(:)
        integer, intent(inout) :: n_failed

        real(dp), parameter :: tol = 1e-5_dp
        integer, parameter :: n_steps = 6
        real(dp) :: steps(n_steps), v_plus(n_outputs), v_minus(n_outputs)
        real(dp) :: analytic(n_outputs), err(n_outputs, n_steps), scale(n_outputs)
        integer :: i, o

        steps = [1e-4_dp, 3e-5_dp, 1e-5_dp, 3e-6_dp, 1e-6_dp, 3e-7_dp]
        analytic = matmul(u, jac)
        scale = matmul(abs(u), abs(jac))
        do i = 1, n_steps
            call evaluate(params + steps(i)*u, template, settings, structure, v_plus)
            call evaluate(params - steps(i)*u, template, settings, structure, v_minus)
            err(:, i) = abs((v_plus - v_minus)/(2.0_dp*steps(i)) - analytic)/ &
                        max(abs(analytic), 1e-6_dp*scale, 1e-300_dp)
        end do
        print "(A,I2,A,7ES10.2)", "  direction ", dir, " best rel err ", &
            minval(err, dim=2)
        print "(A,6ES10.2)", "    nu_star_crit V-curve ", err(i_nu_star_crit, :)
        do o = 1, n_outputs
            if (scale(o) == 0.0_dp) cycle
            if (minval(err(o, :)) > tol) then
                print *, "  FD mismatch direction ", dir, " output ", o, err(o, :)
                n_failed = n_failed + 1
            end if
        end do
    end subroutine fd_check

    subroutine random_direction(v, seed)
        real(dp), intent(out) :: v(:)
        integer, intent(in) :: seed

        integer :: i
        real(dp) :: t

        do i = 1, size(v)
            t = sin(12.9898_dp*real(seed, dp) + 78.233_dp*real(i, dp))*43758.5453_dp
            v(i) = 2.0_dp*(t - floor(t)) - 1.0_dp
        end do
    end subroutine random_direction

end program test_derivative_robustness
