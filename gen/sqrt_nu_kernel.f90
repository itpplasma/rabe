program gen_sqrt_nu_kernel
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, rat, pi_expr, i_expr, &
        operator(+), operator(-), operator(*), operator(/), operator(**), &
        exp, cos, sqrt, log, erfc, abs
    use fortsym_subs, only: subs
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: action, x, a, eta_b, pi, two, profile, phi, roots(3)
    type(expr_t) :: mode, omega, tau_b, tau_d, radius, major_r, nu, ib, bpp, ve
    type(expr_t) :: a0, a1, a2, b1, b2, b3, c1, c2, c3, residual, gamma, ii
    type(expr_t) :: prefactor, amplitude, covariant, speed, bfield, eta, bounce
    type(expr_t) :: iref, l_c, da, de, sa, sb, average_b2, lb, dr, norm, g
    integer :: order
    type(kernel_spec_t) :: spec
    integer :: unit

    call arena%init()
    engine = make_symengine_engine(arena)
    action = sym(arena, 'action')
    x = sym(arena, 'x')
    a = sym(arena, 'a')
    eta_b = sym(arena, 'eta_b')
    pi = pi_expr(arena)
    two = num(arena, 2)
    call verify_identities()
    mode = sym(arena, 'mode')
    omega = sym(arena, 'omega')
    tau_b = sym(arena, 'tau_b')
    tau_d = sym(arena, 'tau_d')
    open(newunit=unit, file='src/generated/sqrt_nu_kernel.f90', status='replace')
    write(unit, '(a)') 'module sqrt_nu_kernel'
    write(unit, '(a)') '    implicit none'
    write(unit, '(a)') '    private'
    write(unit, '(a)') '    public :: wkb_pointwise, log_pointwise, amplitude_pointwise'
    write(unit, '(a)') '    public :: bounce_time_pointwise, e77_boundary_pointwise'
    write(unit, '(a)') 'contains'
    profile = exp(-action)*cos(action)
    roots = [sqrt(abs(mode*omega*tau_b/tau_d)/2), profile, diff(profile, action)]
    spec%name = str('wkb_pointwise')
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str('gen/sqrt_nu_kernel.f90')
    allocate(spec%args(5), spec%outputs(3))
    spec%args = [str('mode'), str('omega'), str('tau_b'), str('tau_d'), str('action')]
    spec%outputs = [str('root'), str('profile'), str('dprofile')]
    call emit(roots)
    deallocate(spec%args, spec%outputs)
    phi = x*sqrt(-log(x)) + sqrt(pi)/2*erfc(sqrt(-log(x)))
    call assert_zero(diff(phi, x) - sqrt(-log(x)), 'logarithmic action derivative')
    roots = [phi, exp(-a*phi/sqrt(two))* &
        cos(a*phi/sqrt(two)), &
        -eta_b/(sqrt(2*a**2*log(a)))]
    spec%name = str('log_pointwise')
    allocate(spec%args(3), spec%outputs(3))
    spec%args = [str('x'), str('a'), str('eta_b')]
    spec%outputs = [str('phi'), str('integrand'), str('asymptotic')]
    call emit(roots)
    deallocate(spec%args, spec%outputs)
    amplitude = sqrt(abs(pi*major_r*mode*ve*eta_b**2/(radius*nu*ib)))* &
        abs(2*eta_b*bpp)**rat(arena, -1_int64, 4_int64)
    spec%name = str('amplitude_pointwise')
    allocate(spec%args(8), spec%outputs(1))
    spec%args = [str('major_r'), str('mode'), str('ve'), str('eta_b'), &
        str('radius'), str('nu'), str('ib'), str('bpp')]
    spec%outputs = [str('amplitude')]
    block
        type(expr_t) :: amplitude_root(1)
        amplitude_root = [amplitude]
        call emit(amplitude_root)
    end block
    deallocate(spec%args, spec%outputs)
    bfield = sym(arena, 'b')
    eta = sym(arena, 'eta')
    bounce = 2*covariant/(speed*bfield*sqrt(1 - eta*bfield))
    roots = [bounce, diff(bounce, bfield), diff(bounce, eta)]
    spec%name = str('bounce_time_pointwise')
    allocate(spec%args(4), spec%outputs(3))
    spec%args = [str('b'), str('eta'), str('covariant'), str('speed')]
    spec%outputs = [str('value'), str('dvalue_db'), str('dvalue_deta')]
    call emit(roots)
    deallocate(spec%args, spec%outputs)
    iref = sym(arena, 'iref')
    l_c = sym(arena, 'l_c')
    da = sym(arena, 'delta_aspect_mode')
    de = sym(arena, 'delta_eta_mode')
    sa = sym(arena, 'sa')
    sb = sym(arena, 'sb')
    average_b2 = sym(arena, 'average_b2')
    lb = sym(arena, 'average_lambda_b')
    dr = sym(arena, 'dr_dAtheta')
    norm = sym(arena, 'average_normalization')
    g = average_b2/lb*(sqrt(eta_b*iref*l_c)*da*sa + l_c*de*sb/2)
    roots = [g, pi*dr/norm, diff(g, l_c)]
    spec%name = str('e77_boundary_pointwise')
    allocate(spec%args(11), spec%outputs(3))
    spec%args = [str('eta_b'), str('iref'), str('l_c'), str('delta_aspect_mode'), &
        str('delta_eta_mode'), str('sa'), str('sb'), str('average_b2'), &
        str('average_lambda_b'), str('dr_dAtheta'), str('average_normalization')]
    spec%outputs = [str('g'), str('normalization'), str('dg_dlc')]
    call emit(roots)
    write(unit, '(a)') 'end module sqrt_nu_kernel'
    close(unit)
contains
    subroutine verify_identities()
        a0 = sym(arena, 'a0')
        a1 = sym(arena, 'a1')
        a2 = sym(arena, 'a2')
        b1 = sym(arena, 'b1')
        b2 = sym(arena, 'b2')
        b3 = sym(arena, 'b3')
        ii = i_expr(arena)
        c1 = ii*a0/b1
        c2 = ((ii*a0 - 2*b2)*c1 + ii*a1)/(4*b1)
        c3 = ((ii*a0 - 6*b2)*c2 + (ii*a1 - 3*b3)*c1 + ii*a2)/(9*b1)
        gamma = 1 + c1*x + c2*x**2 + c3*x**3
        residual = ii*(a0 + a1*x + a2*x**2)*gamma - &
            diff((b1*x + b2*x**2 + b3*x**3)*diff(gamma, x), x)
        do order = 0, 2
            call assert_zero(subs(residual, x, num(arena, 0)), 'series recurrence')
            residual = diff(residual, x)
        end do
        major_r = sym(arena, 'major_r')
        mode = sym(arena, 'mode')
        ve = sym(arena, 've')
        radius = sym(arena, 'radius')
        nu = sym(arena, 'nu')
        ib = sym(arena, 'ib')
        bpp = sym(arena, 'bpp')
        covariant = sym(arena, 'covariant')
        speed = sym(arena, 'speed')
        ! All magnitudes here are positive. Compare fourth powers to avoid
        ! branch assumptions, then select the positive physical A_m root.
        prefactor = (mode*ve*speed/radius/(4*nu*eta_b*covariant*ib/(pi*major_r)))**2 &
            *eta_b**4*(4*covariant/speed)**2*eta_b/(2*bpp)
        amplitude = (pi*major_r*mode*ve*eta_b**2/(radius*nu*ib))**2/(2*eta_b*bpp)
        call assert_zero(prefactor - amplitude, 'A_m fourth powers')
        print *, 'FortSym exact identities: recurrence orders 0..2 and A_m passed'
    end subroutine verify_identities

    subroutine assert_zero(expression, label)
        type(expr_t), intent(in) :: expression
        character(*), intent(in) :: label
        type(engine_result_t) :: verdict
        verdict = engine%zero_test(expression)
        if (verdict%verdict /= VERDICT_TRUE) then
            print *, 'Exact symbolic identity failed or unknown: ', label
            error stop 1
        end if
    end subroutine assert_zero

    subroutine emit(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        type(engine_result_t) :: simplified
        integer :: k
        do k = 1, size(expressions)
            simplified = engine%simplify(expressions(k))
            if (simplified%ok) expressions(k) = simplified%value
        end do
        write(unit, '(a)') chars(emit_kernel(expressions, spec))
    end subroutine emit
end program gen_sqrt_nu_kernel
