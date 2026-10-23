!> Symbolic checks of the identities behind Lambda_alpha and Gamma_well
!! (src/derivative/alpha_criteria.md), with fortsym
!! (https://github.com/lazy-fortran/fortsym). Not part of the rabe build:
!!
!!     FS=<fortsym checkout>; gfortran -I$FS/build/include check_alpha_criteria.f90 \
!!         -L$FS/build/lib -lfortsym -Wl,-rpath,$FS/build/lib -o check && ./check
!!
!! check_identity refuses probe evidence, so every PASS is a decision of the
!! SymEngine zero test.
program check_alpha_criteria
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-)
    use fortsym_parse, only: parse_expr
    use fortsym_diff, only: diff
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: eng
    type(suite_t) :: s

    call arena%init()
    eng = make_symengine_engine(arena)
    call suite_begin(s, "rabe alpha criteria")

    ! dJ/dalpha integrand: d/dtheta [lambda/B] = B_t d/dB [sqrt(1 - eta B)/B]
    ! = -B_t (1 + lambda^2)/(2 lambda B^2), lambda^2 = 1 - eta B
    call check_identity(s, eng, "d/dB (lambda/B) = -(1 + lambda^2)/(2 lambda B^2)", &
        diff(p("sqrt(1 - eta*B)/B"), sym(arena, "B")) &
        + p("(2 - eta*B)/(2*B**2*sqrt(1 - eta*B))"))
    ! the part from d lambda/d theta alone is the bounce-time weight
    call check_identity(s, eng, "d/deta J integrand = -1/(2 lambda)", &
        diff(p("sqrt(1 - eta*B)/B"), sym(arena, "eta")) + p("1/(2*sqrt(1 - eta*B))"))
    ! tangent of h(lambda) = (1 + lambda^2)/(2 lambda) used in the Fortran sweep
    call check_identity(s, eng, "h'(lambda) = (lambda^2 - 1)/(2 lambda^2)", &
        diff(p("(1 + lam**2)/(2*lam)"), sym(arena, "lam")) - p("(lam**2 - 1)/(2*lam**2)"))
    ! lambda^2 next to a turning point: B(a + d) - B(a) per mode,
    ! cos(y + 2c) - cos(y) = -2 sin(y + c) sin(c) with c = (m iota - n) d/2
    call check_identity(s, eng, "cosine difference", &
        p("cos(y + 2*c) - cos(y) + 2*sin(y + c)*sin(c)"))
    ! label variation: delta_1 = sqrt(J_ref/J_1) - 1, J_ref = (3/S)^2,
    ! S = sum_j J_j^(-1/2); the coded tangent of delta_1 with respect to J_1
    call check_identity(s, eng, "d delta_1/d J_1", &
        diff(p("sqrt((3/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2)))**2/J1)"), sym(arena, "J1")) &
        - p("1/2*sqrt((3/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2)))**2/J1)*" // &
            "((-2*9/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2))**3*(-1/2*J1**(-3/2)))/" // &
            "((3/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2)))**2) - 1/J1)"))
    call check_identity(s, eng, "d delta_1/d J_2", &
        diff(p("sqrt((3/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2)))**2/J1)"), sym(arena, "J2")) &
        - p("1/2*sqrt((3/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2)))**2/J1)*" // &
            "((-2*9/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2))**3*(-1/2*J2**(-3/2)))/" // &
            "((3/(J1**(-1/2) + J2**(-1/2) + J3**(-1/2)))**2))"))
    ! turning point tangent: 1 - eta B(phi*(p), p) = 0 differentiated in p
    call check_identity(s, eng, "implicit turning point", &
        p("-(B*deta + eta*dBp)/(eta*Bphi)*eta*Bphi + B*deta + eta*dBp"))

    call suite_end(s, "check_alpha_criteria.json")

contains

    function p(text) result(e)
        character(*), intent(in) :: text
        type(expr_t) :: e
        logical :: ok
        character(:), allocatable :: msg

        e = parse_expr(arena, text, ok, msg)
        if (.not. ok) then
            print *, "parse error: ", text, " ", msg
            error stop 1
        end if
    end function p

end program check_alpha_criteria
