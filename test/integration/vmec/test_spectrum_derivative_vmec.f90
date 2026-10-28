program test_spectrum_derivative_vmec
    !! Derivative mode on a spectrum projected from a real equilibrium
    !! (Landreman-Paul QH or QA reactor-scale low-res references): central
    !! finite differences along random directions of the Boozer spectrum and of
    !! iota, G, I must agree with the analytic Jacobian. Usage:
    !! test_spectrum_derivative_vmec.x qh|qa
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, &
                                   derivative_diagnostics_t, n_outputs, &
                                   n_params, pack_parameters, find_structure, &
                                   evaluate, off_iota, off_G, off_I
    use spectrum_io, only: spectralize_field

    implicit none

    real(dp), parameter :: s_tor = 0.5_dp
    integer, parameter :: mboz = 8, nboz = 8

    type(boozer_field_t) :: field
    type(surface_spectrum_t) :: spec
    type(derivative_settings_t) :: settings
    type(derivative_structure_t) :: structure
    type(derivative_diagnostics_t) :: diag
    character(len=8) :: case_name
    real(dp), allocatable :: params(:), jac(:, :), u(:)
    real(dp) :: values(n_outputs), residual, asymmetry, iota, B_I, B_G
    integer :: K, P, n_failed, dir

    call get_command_argument(1, case_name)
    select case (trim(case_name))
    case ("qh")
        call field%init_from_vmec( &
            "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc", &
            grid_refinement=6)
        settings%M_pol = -1.0_dp
        settings%N_tor = 4.0_dp
    case ("qa")
        call field%init_from_vmec( &
            "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc", &
            grid_refinement=6)
        settings%M_pol = 1.0_dp
        settings%N_tor = 0.0_dp
    case default
        error stop "usage: test_spectrum_derivative_vmec.x qh|qa"
    end select
    settings%sign_sqrtg = -1
    settings%max_n_fieldlines = 60
    settings%n_gauss = 64

    call field%fix_to_surface(s_tor)
    call field%get_iota(s_tor, iota)
    call field%get_covariant_components(B_I, B_G)
    spec%iota = iota
    spec%G = B_G
    spec%I = B_I
    spec%psi_edge = field%psi_tor_edge
    spec%R = field%R
    call spectralize_field(field, field%nfp, mboz, nboz, spec, residual, asymmetry)

    call find_structure(spec, settings, structure)
    ! FD at the input iota with labels frozen at p/q (see the unit test: at the
    ! rational iota itself the S_A, S_B jumps make the primal discontinuous).
    call pack_parameters(spec, params)
    K = size(spec%xm)
    P = n_params(K)
    allocate (jac(P, n_outputs), u(P))
    call evaluate(params, spec, settings, structure, values, jac, diag)
    print "(A,A,I5,A,ES9.2,A,ES9.2)", trim(case_name), ": modes ", K, &
        "  spectral residual ", residual, "  eta_b gap ", diag%eta_b_gap
    print "(A,7ES14.6)", "  values (L_A, L_B, L_S, nu_crit, eta_b, L_alpha, G_well): ", values
    print "(A,2I4,A,2ES10.2)", "  ties eta_b, nu_crit ", diag%eta_b_ties, &
        diag%nu_crit_ties, "  iota kinks ", diag%eta_b_iota_kink, &
        diag%nu_crit_iota_kink

    n_failed = 0
    do dir = 1, 3
        u = 0.0_dp
        select case (dir)
        case (1)
            ! symmetry-breaking modes only, scaled to their size
            call random_direction(u(2:K), 31)
            u(2:K) = u(2:K)*max(abs(spec%bmnc(2:K)), 1e-6_dp*abs(spec%bmnc(1)))
        case (2)
            u(2*K + off_iota) = abs(spec%iota)
        case (3)
            u(2*K + off_G) = spec%G
            u(2*K + off_I) = 0.3_dp*spec%G
            call random_direction(u(K + 1:2*K), 32)
            u(K + 1:2*K) = u(K + 1:2*K)*abs(spec%nabla_s_mnc(1))*1e-2_dp
        end select
        call fd_check(dir, params, spec, settings, structure, jac, u, n_failed)
    end do

    if (n_failed > 0) then
        print *, "FAILED checks: ", n_failed
        error stop 1
    end if
    print *, "test_spectrum_derivative_vmec ", trim(case_name), ": all checks passed"

contains

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
        ! An output whose directional derivative is tiny compared with the sum
        ! of the magnitudes of its contributions is checked in absolute terms.
        scale = matmul(abs(u), abs(jac))
        do i = 1, n_steps
            call evaluate(params + steps(i)*u, template, settings, structure, v_plus)
            call evaluate(params - steps(i)*u, template, settings, structure, v_minus)
            err(:, i) = abs((v_plus - v_minus)/(2.0_dp*steps(i)) - analytic)/ &
                        max(abs(analytic), 1e-6_dp*scale, 1e-300_dp)
        end do
        print "(A,I2,A,7ES10.2)", "  direction ", dir, " best rel err ", &
            minval(err, dim=2)
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

end program test_spectrum_derivative_vmec
