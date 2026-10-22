program test_wout_chain
    !! Derivative mode chained with libneo's VMEC -> Boozer tangent and adjoint
    !! on the Landreman-Paul QH or QA wout (s = 0.5): central finite
    !! differences of the composed primal (wout harmonics perturbed in memory
    !! and re-splined, Boozer spectrum recomputed by libneo, offset coefficients
    !! at the frozen structure; fourth-order stencil) against `wout_chain_jvp`, and adjoint duality of
    !! `wout_chain_jvp` and `wout_chain_vjp`. Usage: test_wout_chain.x qh|qa
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, n_outputs, n_params, &
                                   pack_parameters, find_structure, evaluate, &
                                   off_iota
    use spectrum_io, only: spectralize_field
    use wout_chain, only: wout_chain_t, boozer_wout_tangent_t, wout_tangent_alloc, &
                          wout_chain_init, wout_chain_spectrum, wout_chain_jvp, &
                          wout_chain_vjp

    implicit none

    real(dp), parameter :: s_tor = 0.5_dp
    integer, parameter :: mboz = 6, nboz = 6, grid_refinement = 6

    type(boozer_field_t) :: field
    type(wout_chain_t) :: chain
    type(surface_spectrum_t) :: spec
    type(derivative_settings_t) :: settings
    type(derivative_structure_t) :: structure
    type(boozer_wout_tangent_t) :: t, tb
    character(len=8) :: case_name
    character(len=128) :: wout_file
    real(dp), allocatable :: params(:), jac(:, :, :), jac_frozen_iota(:, :, :)
    real(dp) :: values(n_outputs), residual, asymmetry, iota, B_I, B_G
    real(dp) :: dv(n_outputs, 1), cot(n_outputs, 1), lhs, rhs, iota_line, iota_ad
    integer :: n_failed, K

    call get_command_argument(1, case_name)
    select case (trim(case_name))
    case ("qh")
        wout_file = "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
        settings%M_pol = -1.0_dp
        settings%N_tor = 4.0_dp
    case ("qa")
        wout_file = "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
        settings%M_pol = 1.0_dp
        settings%N_tor = 0.0_dp
    case default
        error stop "usage: test_wout_chain.x qh|qa"
    end select
    settings%sign_sqrtg = -1
    settings%max_n_fieldlines = 60
    settings%n_gauss = 64

    ! nabla_s_mnc, psi_edge and R from the rabe projection, as in rabe_deriv.x
    call field%init_from_vmec(trim(wout_file), grid_refinement=grid_refinement)
    call field%fix_to_surface(s_tor)
    call field%get_iota(s_tor, iota)
    call field%get_covariant_components(B_I, B_G)
    spec%psi_edge = field%psi_tor_edge
    spec%R = field%R
    call spectralize_field(field, field%nfp, mboz, nboz, spec, residual, asymmetry)

    call wout_chain_init(chain, trim(wout_file), mboz, nboz, grid_refinement)
    call wout_chain_spectrum(chain, s_tor, spec)
    iota_ad = spec%iota
    print "(A,A,3ES12.4)", trim(case_name), ": iota, G, I libneo AD ", &
        spec%iota, spec%G, spec%I
    print "(A,3ES12.4)", "    iota, G, I rabe field  ", iota, B_G, B_I

    call find_structure(spec, settings, structure)
    iota_line = structure%iota_line
    spec%iota = iota_line
    call pack_parameters(spec, params)
    K = size(spec%xm)
    allocate (jac(n_params(K), n_outputs, 1))
    call evaluate(params, spec, settings, structure, values, jac(:, :, 1))
    print "(A,5ES13.5)", "    values ", values

    n_failed = 0

    ! Geometry direction (rmnc, zmns, lmns), iotaf fixed: the full chain.
    call wout_tangent_alloc(chain%ctx, t)
    call geometry_direction(t)
    call wout_chain_jvp(chain, [s_tor], jac, t, dv)
    call fd_check("rmnc, zmns, lmns", t, dv(:, 1), &
                  abs_scale(jac(:, :, 1), t), .false.)

    ! iotaf direction with the field-line iota held at p/q: the chain without
    ! the iota column (bmnc, G and I still depend on iota through libneo).
    call wout_tangent_alloc(chain%ctx, t)
    call iota_direction(t)
    allocate (jac_frozen_iota, source=jac)
    jac_frozen_iota(2*K + off_iota, :, :) = 0.0_dp
    call wout_chain_jvp(chain, [s_tor], jac_frozen_iota, t, dv)
    call fd_check("iotaf (line iota frozen)", t, dv(:, 1), &
                  abs_scale(jac_frozen_iota(:, :, 1), t), .true.)

    ! Adjoint duality with all inputs, including the iota column.
    call wout_tangent_alloc(chain%ctx, t)
    call geometry_direction(t)
    call iota_direction(t)
    call wout_chain_jvp(chain, [s_tor], jac, t, dv)
    cot(:, 1) = [0.7_dp, -1.3_dp, 2.1_dp, 0.4_dp, -0.9_dp]
    call wout_chain_vjp(chain, [s_tor], jac, cot, tb)
    lhs = sum(cot*dv)
    rhs = sum(tb%rmnc*t%rmnc) + sum(tb%zmns*t%zmns) + sum(tb%lmns*t%lmns) &
          + sum(tb%iotaf*t%iotaf)
    print "(A,ES10.2)", "  adjoint duality rel err ", abs(lhs - rhs)/abs(lhs)
    if (abs(lhs - rhs) > 1e-12_dp*abs(lhs)) n_failed = n_failed + 1

    if (n_failed > 0) then
        print *, "FAILED checks: ", n_failed
        error stop 1
    end if
    print *, "test_wout_chain ", trim(case_name), ": all checks passed"

contains

    subroutine geometry_direction(t)
        type(boozer_wout_tangent_t), intent(inout) :: t
        integer :: i, j
        real(dp) :: s, decay

        ! smooth in s and decaying with m and n like a boundary perturbation
        do j = 0, chain%ctx%kpar
            s = real(j, dp)/real(chain%ctx%kpar, dp)
            do i = 1, chain%ctx%nstrm
                decay = 1.0_dp/(real(1 + chain%ctx%xm(i), dp)**2* &
                                real(1 + abs(chain%ctx%xn(i))/chain%ctx%nper, dp)**2)
                t%rmnc(i, j) = 1e-3_dp*s*cos(0.7_dp*i)*decay
                t%zmns(i, j) = 1e-3_dp*s*sin(1.3_dp*i)*decay
                t%lmns(i, j) = 1e-2_dp*s*cos(2.1_dp*i)*decay
            end do
        end do
    end subroutine geometry_direction

    subroutine iota_direction(t)
        type(boozer_wout_tangent_t), intent(inout) :: t
        integer :: j
        real(dp) :: s

        do j = 0, chain%ctx%kpar
            s = real(j, dp)/real(chain%ctx%kpar, dp)
            t%iotaf(j) = 0.05_dp*(1.0_dp - s**2)
        end do
    end subroutine iota_direction

    !> |J|^T |d spectrum| for outputs whose derivative cancels: checked in
    !! absolute terms relative to the size of their contributions.
    function abs_scale(jac1, t) result(scale)
        real(dp), intent(in) :: jac1(:, :)
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp) :: scale(n_outputs)
        real(dp) :: d(n_outputs, 1)

        call wout_chain_jvp(chain, [s_tor], reshape(abs(jac1), [size(jac1, 1), &
                            n_outputs, 1]), abs_tangent(t), d)
        scale = abs(d(:, 1))
    end function abs_scale

    function abs_tangent(t) result(ta)
        type(boozer_wout_tangent_t), intent(in) :: t
        type(boozer_wout_tangent_t) :: ta

        ta = t
        ta%rmnc = abs(t%rmnc)
        ta%zmns = abs(t%zmns)
        ta%lmns = abs(t%lmns)
        ta%iotaf = abs(t%iotaf)
    end function abs_tangent

    !> Fourth-order central differences of the composed primal along `t`.
    !! libneo's re-splined primal carries a roundoff floor of about 1e-11
    !! relative (bmnc FD error grows like 1e-11/h below h = 3e-4), so
    !! second-order differences cannot resolve the outputs to better than a
    !! few 1e-6; the fourth-order stencil allows larger steps. The tolerance is
    !! that of the other rabe FD tests. The smallest outputs (QA Lambda_A,
    !! Lambda_B of order 1e-8) are noise-limited at about 5e-6.
    subroutine fd_check(label, t, analytic, scale, freeze_iota)
        character(len=*), intent(in) :: label
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(in) :: analytic(n_outputs), scale(n_outputs)
        logical, intent(in) :: freeze_iota

        real(dp), parameter :: tol = 1e-5_dp
        integer, parameter :: n_steps = 2
        real(dp) :: steps(n_steps), fd(n_outputs), err(n_outputs, n_steps)
        integer :: i, o

        steps = [1e-3_dp, 3e-4_dp]
        do i = 1, n_steps
            fd = (8.0_dp*(composed(t, steps(i), freeze_iota) &
                          - composed(t, -steps(i), freeze_iota)) &
                  - (composed(t, 2.0_dp*steps(i), freeze_iota) &
                     - composed(t, -2.0_dp*steps(i), freeze_iota)))/(12.0_dp*steps(i))
            err(:, i) = abs(fd - analytic)/max(abs(analytic), 1e-6_dp*scale, 1e-300_dp)
        end do
        call respline(t, 0.0_dp)
        print "(A,A30,A,5ES10.2)", "  ", label, " best rel err ", minval(err, dim=2)
        do i = 1, n_steps
            print "(A,ES8.1,A,5ES10.2)", "    step ", steps(i), " rel err ", err(:, i)
        end do
        do o = 1, n_outputs
            if (scale(o) == 0.0_dp) cycle
            if (minval(err(o, :)) > tol) then
                print *, "  FD mismatch ", label, " output ", o, err(o, :)
                n_failed = n_failed + 1
            end if
        end do
    end subroutine fd_check

    function composed(t, eps, freeze_iota) result(v)
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(in) :: eps
        logical, intent(in) :: freeze_iota
        real(dp) :: v(n_outputs)
        type(surface_spectrum_t) :: sp
        real(dp), allocatable :: p(:)

        call respline(t, eps)
        sp = spec
        call wout_chain_spectrum(chain, s_tor, sp)
        if (freeze_iota) then
            sp%iota = iota_line
        else
            sp%iota = iota_line + (sp%iota - iota_ad)
        end if
        call pack_parameters(sp, p)
        call evaluate(p, sp, settings, structure, v)
    end function composed

    !> Re-spline libneo's VMEC data with the wout harmonics shifted by eps*t,
    !! as spline_vmec_data does after reading the file (same procedure as
    !! libneo test/source/test_boozer_vmec_ad.f90).
    subroutine respline(t, eps)
        use new_vmec_stuff_mod, only: rmnc, zmns, almns, aiota, kpar
        use vmec_alloc_sub, only: new_allocate_vmec_stuff, new_deallocate_vmec_stuff
        use spline_vmec_sub, only: initialize_vmec_data_and_arrays, &
                                   setup_poloidal_flux_splines, &
                                   setup_angular_grid_and_fourier_synthesis, &
                                   spline_over_phi, spline_over_theta, spline_over_s
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(in) :: eps
        integer :: i

        call new_allocate_vmec_stuff()
        call initialize_vmec_data_and_arrays()
        rmnc = rmnc + eps*chain%ctx%fac_r*t%rmnc
        zmns = zmns + eps*chain%ctx%fac_r*t%zmns
        do i = 1, kpar - 1
            almns(:, i) = almns(:, i) + eps*0.5_dp*(t%lmns(:, i + 1) + t%lmns(:, i))
        end do
        almns(:, kpar) = almns(:, kpar) + eps*(1.5_dp*t%lmns(:, kpar) &
                                               - 0.5_dp*t%lmns(:, kpar - 1))
        aiota = aiota + eps*t%iotaf
        call setup_poloidal_flux_splines()
        call setup_angular_grid_and_fourier_synthesis()
        call new_deallocate_vmec_stuff()
        call spline_over_phi()
        call spline_over_theta()
        call spline_over_s()
    end subroutine respline

end program test_wout_chain
