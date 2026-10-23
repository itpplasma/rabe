program test_wout_chain
    !! Derivative mode chained with libneo's VMEC -> Boozer tangent and adjoint
    !! on the Landreman-Paul QH or QA wout (s = 0.5): central finite
    !! differences of the composed primal (wout data perturbed in memory by
    !! libneo's `boozer_ad_perturb`, every spectrum parameter recomputed by
    !! libneo, offset coefficients at the frozen structure; fourth-order
    !! stencil) against `wout_chain_jvp` for geometry, iotaf, phiedge and
    !! rmajor directions, and adjoint duality of `wout_chain_jvp` and
    !! `wout_chain_vjp`. libneo's |grad s| spectrum, psi_edge and R are also
    !! compared with rabe's own field. Usage: test_wout_chain.x qh|qa
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, n_outputs, n_params, &
                                   pack_parameters, find_structure, evaluate, &
                                   off_iota
    use spectrum_io, only: spectralize_field
    use boozer_vmec_ad, only: boozer_ad_perturb
    use wout_chain, only: wout_chain_t, boozer_wout_tangent_t, wout_tangent_alloc, &
                          wout_chain_init, wout_chain_spectrum, wout_chain_jvp, &
                          wout_chain_vjp

    implicit none

    real(dp), parameter :: s_tor = 0.5_dp
    integer, parameter :: mboz = 6, nboz = 6, grid_refinement = 6

    type(boozer_field_t) :: field
    type(wout_chain_t) :: chain
    type(surface_spectrum_t) :: spec, spec_field
    type(derivative_settings_t) :: settings
    type(derivative_structure_t) :: structure
    type(boozer_wout_tangent_t) :: t, tb
    character(len=8) :: case_name
    character(len=128) :: wout_file
    real(dp), allocatable :: params(:), jac(:, :, :), jac_frozen_iota(:, :, :)
    real(dp) :: values(n_outputs), residual, asymmetry, iota, B_I, B_G
    real(dp) :: dv(n_outputs, 1), cot(n_outputs, 1), lhs, rhs, iota_line, iota_ad
    real(dp) :: scale_phi(n_outputs)
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

    n_failed = 0

    ! rabe's own field and projection, to check the spectrum parameters
    ! that libneo now provides (rabe_deriv.x uses libneo's).
    call field%init_from_vmec(trim(wout_file), grid_refinement=grid_refinement)
    call field%fix_to_surface(s_tor)
    call field%get_iota(s_tor, iota)
    call field%get_covariant_components(B_I, B_G)
    call spectralize_field(field, field%nfp, mboz, nboz, spec_field, residual, &
                           asymmetry)

    call wout_chain_init(chain, trim(wout_file), mboz, nboz, grid_refinement)
    call wout_chain_spectrum(chain, s_tor, spec)
    iota_ad = spec%iota
    print "(A,A,3ES12.4)", trim(case_name), ": iota, G, I libneo AD ", &
        spec%iota, spec%G, spec%I
    print "(A,3ES12.4)", "    iota, G, I rabe field  ", iota, B_G, B_I
    call compare_with_field()

    call find_structure(spec, settings, structure)
    iota_line = structure%iota_line
    spec%iota = iota_line
    call pack_parameters(spec, params)
    K = size(spec%xm)
    allocate (jac(n_params(K), n_outputs, 1))
    call evaluate(params, spec, settings, structure, values, jac(:, :, 1))
    print "(A,5ES13.5)", "    values ", values

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

    ! Toroidal flux (whole phi profile scaled): B, G, I and psi_edge scale
    ! together, which leaves Lambda_A, Lambda_B, Lambda_S and nu_star_crit
    ! invariant and scales eta_b like 1/B. The invariant outputs are a
    ! cancellation of contributions of size `scale`, so the analytic tangent
    ! must vanish relative to `scale`, and the FD is compared in absolute
    ! terms relative to `scale` as well.
    call wout_tangent_alloc(chain%ctx, t)
    t%phiedge = chain%ctx%phiedge
    call wout_chain_jvp(chain, [s_tor], jac, t, dv)
    scale_phi = abs_scale(jac(:, :, 1), t)
    print "(A,4ES10.2)", "    phiedge: |analytic|/scale of the invariant outputs ", &
        abs(dv(1:4, 1))/scale_phi(1:4)
    if (any(abs(dv(1:4, 1)) > 1e-10_dp*scale_phi(1:4))) then
        print *, "  phiedge: invariant outputs have a nonzero tangent"
        n_failed = n_failed + 1
    end if
    call fd_check("phiedge", t, dv(:, 1), scale_phi, .false., &
                  absolute=[.true., .true., .true., .true., .false.])
    ! Rmajor_p, read from the file and passed through to R.
    call wout_tangent_alloc(chain%ctx, t)
    t%rmajor = 1e-2_dp*chain%ctx%rmajor
    call wout_chain_jvp(chain, [s_tor], jac, t, dv)
    call fd_check("rmajor", t, dv(:, 1), abs_scale(jac(:, :, 1), t), .false.)

    ! Adjoint duality with all inputs, including the iota column.
    call wout_tangent_alloc(chain%ctx, t)
    call geometry_direction(t)
    call iota_direction(t)
    t%phiedge = 3e-3_dp*chain%ctx%phiedge
    t%rmajor = -2e-3_dp*chain%ctx%rmajor
    call wout_chain_jvp(chain, [s_tor], jac, t, dv)
    cot(:, 1) = [0.7_dp, -1.3_dp, 2.1_dp, 0.4_dp, -0.9_dp]
    call wout_chain_vjp(chain, [s_tor], jac, cot, tb)
    lhs = sum(cot*dv)
    rhs = sum(tb%rmnc*t%rmnc) + sum(tb%zmns*t%zmns) + sum(tb%lmns*t%lmns) &
          + sum(tb%iotaf*t%iotaf) + tb%phiedge*t%phiedge + tb%rmajor*t%rmajor
    print "(A,ES10.2)", "  adjoint duality rel err ", abs(lhs - rhs)/abs(lhs)
    if (abs(lhs - rhs) > 1e-12_dp*abs(lhs)) n_failed = n_failed + 1

    if (n_failed > 0) then
        print *, "FAILED checks: ", n_failed
        error stop 1
    end if
    print *, "test_wout_chain ", trim(case_name), ": all checks passed"

contains

    !> libneo's |grad s| spectrum (1/m), psi_edge and R against rabe's field
    !! (projection of compute_nabla_s, psi_tor_edge, R): a check of units,
    !! sign and mode convention. The spectra come from different
    !! discretizations (libneo's transform vs rabe's projection of its Boozer
    !! splines). The leading modes agree to five digits, and the largest
    !! difference sits in small high-|n| modes, 2e-6 (QH) and 1.3e-5 (QA) of
    !! the (0,0) mode, like the bmnc difference of the two transforms; a unit
    !! or convention error would be of order one.
    subroutine compare_with_field()
        real(dp) :: d_grads, d_psi, d_R

        d_grads = maxval(abs(spec%nabla_s_mnc - spec_field%nabla_s_mnc)) &
                  /abs(spec_field%nabla_s_mnc(1))
        d_psi = abs(spec%psi_edge - field%psi_tor_edge)/abs(field%psi_tor_edge)
        d_R = abs(spec%R - field%R)
        print "(A,3ES10.2)", "    |grad s| spectrum, psi_edge, R vs rabe field ", &
            d_grads, d_psi, d_R
        if (d_grads > 1e-4_dp .or. d_psi > 1e-12_dp .or. d_R > 0.0_dp) then
            print *, "  libneo spectrum parameters differ from rabe's field"
            n_failed = n_failed + 1
        end if
    end subroutine compare_with_field

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
        ta%phiedge = abs(t%phiedge)
        ta%rmajor = abs(t%rmajor)
    end function abs_tangent

    !> Fourth-order central differences of the composed primal along `t`.
    !! libneo's primal carries a roundoff floor of about 1e-11
    !! relative (bmnc FD error grows like 1e-11/h below h = 3e-4), so
    !! second-order differences cannot resolve the outputs to better than a
    !! few 1e-6; the fourth-order stencil allows larger steps. The tolerance is
    !! that of the other rabe FD tests. The smallest outputs (QA Lambda_A,
    !! Lambda_B of order 1e-8) are noise-limited at about 5e-6.
    subroutine fd_check(label, t, analytic, scale, freeze_iota, absolute)
        character(len=*), intent(in) :: label
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(in) :: analytic(n_outputs), scale(n_outputs)
        logical, intent(in) :: freeze_iota
        logical, intent(in), optional :: absolute(n_outputs)
            !! error relative to `scale` instead of the tangent (invariants)

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
            if (present(absolute)) then
                where (absolute) err(:, i) = abs(fd - analytic)/max(scale, 1e-300_dp)
            end if
        end do
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
        type(wout_chain_t) :: shifted
        real(dp), allocatable :: p(:)

        ! libneo's transform reads the wout data held by its context, not the
        ! VMEC splines, so the shift goes through boozer_ad_perturb on a copy.
        shifted = chain
        call boozer_ad_perturb(shifted%ctx, t, eps)
        sp = spec
        call wout_chain_spectrum(shifted, s_tor, sp)
        if (freeze_iota) then
            sp%iota = iota_line
        else
            sp%iota = iota_line + (sp%iota - iota_ad)
        end if
        call pack_parameters(sp, p)
        call evaluate(p, sp, settings, structure, v)
    end function composed

end program test_wout_chain
