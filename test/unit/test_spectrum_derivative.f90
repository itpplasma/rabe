program test_spectrum_derivative
    !! Finite-difference and adjoint checks of the derivative mode on synthetic
    !! QH and QA spectra. The primal in derivative mode is smooth at frozen
    !! structure, so central differences converge as h^2 until roundoff, and the
    !! best step of a sweep (V curve) must agree with the analytic Jacobian.
    use constants, only: dp
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, &
                                   derivative_diagnostics_t, n_outputs, &
                                   n_params, pack_parameters, find_structure, &
                                   evaluate, spectrum_vjp, spectrum_jvp, &
                                   gauss_legendre_01, off_iota, off_G, off_I, &
                                   off_psi_edge, off_R
    use synthetic_spectra, only: synthetic_qh, synthetic_qa

    implicit none

    type(surface_spectrum_t) :: spec
    type(derivative_settings_t) :: settings
    integer :: n_failed

    n_failed = 0
    call check_gauss_legendre(n_failed)

    call synthetic_qh(spec, settings)
    call check_case("QH", spec, settings, n_failed)
    call synthetic_qa(spec, settings)
    call check_case("QA", spec, settings, n_failed)

    if (n_failed > 0) then
        print *, "FAILED checks: ", n_failed
        error stop 1
    end if
    print *, "test_spectrum_derivative: all checks passed"

contains

    subroutine check_gauss_legendre(n_failed)
        integer, intent(inout) :: n_failed

        real(dp) :: x(32), w(32)
        integer :: p
        real(dp) :: err

        call gauss_legendre_01(32, x, w)
        do p = 0, 63
            err = abs(sum(w*x**p) - 1.0_dp/real(p + 1, dp))
            if (err > 1e-14_dp) then
                print *, "Gauss-Legendre moment ", p, " error ", err
                n_failed = n_failed + 1
            end if
        end do
    end subroutine check_gauss_legendre

    subroutine check_case(name, spec, settings, n_failed)
        character(len=*), intent(in) :: name
        type(surface_spectrum_t), intent(in) :: spec
        type(derivative_settings_t), intent(in) :: settings
        integer, intent(inout) :: n_failed

        type(derivative_structure_t) :: structure
        type(derivative_diagnostics_t) :: diag
        real(dp), allocatable :: params(:), jac(:, :), scale(:), u(:), grad(:)
        real(dp) :: values(n_outputs), vals_jvp(n_outputs), jvp(n_outputs)
        real(dp) :: cot(n_outputs)
        integer :: K, P, o, block
        character(len=16) :: block_name(5)

        K = size(spec%xm)
        P = n_params(K)
        call find_structure(spec, settings, structure)
        call pack_parameters(spec, params)
        ! The FD checks run at the input iota with the labels frozen at p/q,
        ! not at the rational iota itself: there the modes with m*iota_p = 0
        ! mod pi sit on the jumps of S_A and S_B, the primal is discontinuous by
        ! about resonant_weight (PHILOSOPHY item 2), and central differences in
        ! iota diverge like 1/h although the surrogate derivative is defined.
        allocate (jac(P, n_outputs))
        call evaluate(params, spec, settings, structure, values, jac, diag)

        print "(A,A,I4,A,I4,A,I4)", name, ": lines ", structure%n_lines, &
            "  q ", structure%q, "  split lines ", structure%n_split_lines
        print "(A,5ES14.6)", "  values (L_A, L_B, L_S, nu_crit, eta_b): ", values
        print "(A,ES10.2,A,ES10.2,A,2ES10.2)", "  newton shift ", diag%max_newton_shift, &
            "  eta_b gap ", diag%eta_b_gap, "  resonant weight A, B ", &
            diag%resonant_weight_A, diag%resonant_weight_B
        print "(A,2ES10.2)", "  continuum shift eta_b, nu_crit ", diag%eta_b_shift, &
            diag%nu_crit_shift
        if (any(values(1:4) == 0.0_dp)) then
            print *, name, ": an output vanishes, test spectrum not generic"
            n_failed = n_failed + 1
        end if

        ! Parameter scales for the directions: relative for the mode blocks.
        allocate (scale(P), u(P), grad(P))
        scale(1:K) = max(abs(spec%bmnc), 1e-3_dp)
        scale(K + 1:2*K) = max(abs(spec%nabla_s_mnc), 1e-2_dp)
        scale(2*K + 1:P) = abs(params(2*K + 1:P))

        block_name = [character(len=16) :: "bmnc", "nabla_s_mnc", "iota", "G,I", &
                      "psi_edge,R"]
        do block = 1, 6
            u = 0.0_dp
            select case (block)
            case (1)
                call random_direction(u(1:K), 11)
            case (2)
                call random_direction(u(K + 1:2*K), 12)
            case (3)
                u(2*K + off_iota) = 1.0_dp
            case (4)
                u(2*K + off_G) = 1.0_dp
                u(2*K + off_I) = -0.7_dp
            case (5)
                u(2*K + off_psi_edge) = 1.0_dp
                u(2*K + off_R) = 0.5_dp
            case (6)
                call random_direction(u, 13)
            end select
            u = u*scale
            if (block <= 5) then
                call fd_check(name//" "//trim(block_name(block)), params, spec, &
                              settings, structure, jac, u, n_failed)
            else
                call fd_check(name//" random", params, spec, settings, &
                              structure, jac, u, n_failed)
            end if
        end do

        ! Discrete-selection mode (opt-out of continuum extrema): smooth in the
        ! spectrum at fixed rational iota, where tied maxima are the same point.
        block
            type(derivative_settings_t) :: lattice
            real(dp), allocatable :: jac_lattice(:, :)
            real(dp) :: values_lattice(n_outputs)

            lattice = settings
            lattice%continuum_extrema = .false.
            allocate (jac_lattice(P, n_outputs))
            call evaluate(params, spec, lattice, structure, values_lattice, jac_lattice)
            u = 0.0_dp
            call random_direction(u(1:K), 14)
            u = u*scale
            call fd_check(name//" lattice bmnc", params, spec, lattice, structure, &
                          jac_lattice, u, n_failed)
        end block

        ! JVP and VJP entries agree with the dense Jacobian (adjoint duality).
        call random_direction(cot, 21)
        call spectrum_vjp(params, spec, settings, structure, cot, vals_jvp, grad)
        call random_direction(u, 22)
        u = u*scale
        call spectrum_jvp(params, spec, settings, structure, u, vals_jvp, jvp)
        if (abs(dot_product(grad, u) - dot_product(cot, jvp)) > &
            1e-12_dp*max(abs(dot_product(cot, jvp)), tiny(1.0_dp))) then
            print *, name, ": adjoint duality violated ", dot_product(grad, u), &
                dot_product(cot, jvp)
            n_failed = n_failed + 1
        end if
        print "(A,A,ES12.4)", name, ": adjoint duality <v,Ju>-<J^Tv,u> rel ", &
            abs(dot_product(grad, u) - dot_product(cot, jvp))/ &
            abs(dot_product(cot, jvp))
        do o = 1, n_outputs
            if (vals_jvp(o) /= values(o)) then
                print *, name, ": primal not reproducible, output ", o
                n_failed = n_failed + 1
            end if
        end do
    end subroutine check_case

    !> Central differences over a step sweep; the best step must agree with
    !! the analytic directional derivative for every output.
    subroutine fd_check(label, params, template, settings, structure, jac, u, &
                        n_failed)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: params(:)
        type(surface_spectrum_t), intent(in) :: template
        type(derivative_settings_t), intent(in) :: settings
        type(derivative_structure_t), intent(in) :: structure
        real(dp), intent(in) :: jac(:, :), u(:)
        integer, intent(inout) :: n_failed

        real(dp), parameter :: tol = 1e-6_dp
        integer, parameter :: n_steps = 7
        real(dp) :: steps(n_steps), v_plus(n_outputs), v_minus(n_outputs)
        real(dp) :: analytic(n_outputs), fd(n_outputs), err(n_outputs, n_steps)
        real(dp) :: best(n_outputs), denom
        integer :: i, o

        steps = [1e-3_dp, 3e-4_dp, 1e-4_dp, 3e-5_dp, 1e-5_dp, 3e-6_dp, 1e-6_dp]
        analytic = matmul(u, jac)
        do i = 1, n_steps
            call evaluate(params + steps(i)*u, template, settings, structure, v_plus)
            call evaluate(params - steps(i)*u, template, settings, structure, v_minus)
            fd = (v_plus - v_minus)/(2.0_dp*steps(i))
            do o = 1, n_outputs
                denom = max(abs(analytic(o)), 1e-30_dp)
                err(o, i) = abs(fd(o) - analytic(o))/denom
            end do
        end do
        best = minval(err, dim=2)
        print "(2X,A,T26,A,5ES10.2)", label, " best rel err ", best
        do o = 1, n_outputs
            if (abs(analytic(o)) < 1e-14_dp*maxval(abs(jac(:, o)))*maxval(abs(u))) cycle
            if (best(o) > tol) then
                print *, "  FD mismatch ", label, " output ", o, " errors ", err(o, :)
                n_failed = n_failed + 1
            end if
        end do
    end subroutine fd_check

    !> Deterministic pseudo-random direction in [-1, 1].
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

end program test_spectrum_derivative
