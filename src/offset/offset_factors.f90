module offset_factors
    use constants, only: dp, pi
    use fit_functions, only: S_A_fit => S_A, S_B_fit => S_B

    implicit none

    private

    public :: S_A, S_B, set_offset_factor_mode
    public :: base_wiener_hopf_residual

    !> Discretization of the universal Wiener-Hopf offset-factor problem.
    !! Reproduces the Python oracle tools/offset_factor_solver.py on the
    !! documented check grid (points=121, extent=7.0).
    type :: offset_grid_t
        integer :: points
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: kernel(:, :)
        real(dp), allocatable :: phi(:)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: phi_a(:)
        real(dp), allocatable :: phi_b(:)
    end type offset_grid_t

    !> Grid matching the documented oracle reference table. The factors are not
    !! converged in the truncation extent, so these values are load-bearing.
    integer, parameter :: oracle_points = 121
    real(dp), parameter :: oracle_extent = 7.0_dp
    integer, parameter :: phi_source_order = 96

    logical :: use_exact = .false.
    logical :: base_ready = .false.
    type(offset_grid_t) :: base

contains

    !> Select between the analytic fit (default) and the exact Wiener-Hopf
    !! factors. Mirrors set_unsafe_mode: a namelist flag applied once in main.
    subroutine set_offset_factor_mode(exact)
        logical, intent(in) :: exact

        use_exact = exact
    end subroutine set_offset_factor_mode

    !> Universal offset factor S_A evaluated at each angle. Branches on the
    !! module mode flag; the fit path is the elemental fit_functions routine.
    function S_A(angle) result(out)
        real(dp), intent(in) :: angle(:)
        real(dp) :: out(size(angle))

        if (use_exact) then
            call ensure_base()
            out = exact_factors(angle, want_b=.false.)
        else
            out = S_A_fit(angle)
        end if
    end function S_A

    !> Universal offset factor S_B evaluated at each angle.
    function S_B(angle) result(out)
        real(dp), intent(in) :: angle(:)
        real(dp) :: out(size(angle))

        if (use_exact) then
            call ensure_base()
            out = exact_factors(angle, want_b=.true.)
        else
            out = S_B_fit(angle)
        end if
    end function S_B

    !> Build the angle-independent base solution once and cache it.
    subroutine ensure_base()
        if (base_ready) return
        call solve_base(oracle_points, oracle_extent, base)
        base_ready = .true.
    end subroutine ensure_base

    !> Evaluate the exact factor at each angle. want_b selects S_B, else S_A.
    function exact_factors(angle, want_b) result(out)
        real(dp), intent(in) :: angle(:)
        logical, intent(in) :: want_b
        real(dp) :: out(size(angle))

        integer :: k

        do k = 1, size(angle)
            out(k) = sample_factor(base, angle(k), want_b)
        end do
    end function exact_factors

    !> Residual max|alpha - (K (sign alpha) + phi)| of the base Wiener-Hopf
    !! solve on the cached grid. Used by the unit test to assert convergence.
    function base_wiener_hopf_residual() result(residual)
        real(dp) :: residual

        real(dp), allocatable :: rhs(:), signed(:)
        integer :: i

        call ensure_base()
        allocate (signed(base%points))
        do i = 1, base%points
            signed(i) = sign_of(base%x(i))*base%alpha(i)
        end do
        rhs = matmul(base%kernel, signed) + base%phi
        residual = maxval(abs(base%alpha - rhs))
    end function base_wiener_hopf_residual

    !> Uniform grid on [-extent, extent] with trapezoidal weights.
    subroutine trapezoid_grid(points, extent, x, weights)
        integer, intent(in) :: points
        real(dp), intent(in) :: extent
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable, intent(out) :: weights(:)

        integer :: i
        real(dp) :: h

        allocate (x(points))
        allocate (weights(points))
        h = 2.0_dp*extent/real(points - 1, dp)
        do i = 1, points
            x(i) = -extent + real(i - 1, dp)*h
        end do
        weights = h
        weights(1) = 0.5_dp*h
        weights(points) = 0.5_dp*h
    end subroutine trapezoid_grid

    !> Right-hand side phi_i = 2/sqrt(pi) int_0^1 exp(-x_i^2/t^2) dt, evaluated
    !! with Gauss-Legendre quadrature mapped from [-1,1] to [0,1].
    subroutine phi_source(x, phi)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: phi(:)

        real(dp) :: nodes(phi_source_order), gw(phi_source_order)
        real(dp) :: t(phi_source_order), wq(phi_source_order)
        integer :: i, k
        real(dp) :: acc

        call gauss_legendre(phi_source_order, nodes, gw)
        t = 0.5_dp*(nodes + 1.0_dp)
        wq = 0.5_dp*gw

        allocate (phi(size(x)))
        do i = 1, size(x)
            acc = 0.0_dp
            do k = 1, phi_source_order
                acc = acc + exp(-(x(i)*x(i))/(t(k)*t(k)))*wq(k)
            end do
            phi(i) = 2.0_dp/sqrt(pi)*acc
        end do
    end subroutine phi_source

    !> Column-scaled Gaussian kernel K_ij = exp(-(x_i-x_j)^2) w_j / sqrt(pi).
    subroutine gaussian_kernel(x, weights, kernel)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: weights(:)
        real(dp), allocatable, intent(out) :: kernel(:, :)

        integer :: i, j, n
        real(dp) :: delta

        n = size(x)
        allocate (kernel(n, n))
        do j = 1, n
            do i = 1, n
                delta = x(i) - x(j)
                kernel(i, j) = exp(-(delta*delta))*weights(j)/sqrt(pi)
            end do
        end do
    end subroutine gaussian_kernel

    !> Solve the base Wiener-Hopf equation and derive the angle-independent
    !! source vectors phi_a and phi_b.
    subroutine solve_base(points, extent, grid)
        integer, intent(in) :: points
        real(dp), intent(in) :: extent
        type(offset_grid_t), intent(out) :: grid

        real(dp), allocatable :: matrix(:, :), phi(:), alpha(:)
        integer, allocatable :: ipiv(:)
        real(dp) :: sign_j, alpha0
        integer :: i, j, info

        grid%points = points
        call trapezoid_grid(points, extent, grid%x, grid%weights)
        call gaussian_kernel(grid%x, grid%weights, grid%kernel)
        call phi_source(grid%x, phi)

        allocate (matrix(points, points))
        allocate (alpha(points))
        allocate (ipiv(points))
        do j = 1, points
            sign_j = sign_of(grid%x(j))
            do i = 1, points
                matrix(i, j) = -grid%kernel(i, j)*sign_j
            end do
            matrix(j, j) = matrix(j, j) + 1.0_dp
        end do

        alpha = phi
        call dgesv(points, 1, matrix, points, ipiv, alpha, points, info)
        if (info /= 0) then
            print *, "offset_factors: base dgesv failed, info = ", info
            error stop
        end if

        grid%phi = phi
        grid%alpha = alpha

        alpha0 = interp_zero(grid%x, alpha)

        allocate (grid%phi_b(points))
        grid%phi_b = 2.0_dp*alpha0/sqrt(pi)*exp(-(grid%x*grid%x))

        call aspect_source(grid%x, grid%weights, alpha, grid%phi_a)
    end subroutine solve_base

    !> Source vector phi_a from the aspect-ratio kernel acting on the base alpha.
    subroutine aspect_source(x, weights, alpha, phi_a)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: alpha(:)
        real(dp), allocatable, intent(out) :: phi_a(:)

        integer :: i, j, n
        real(dp) :: delta, k2, src, acc, sign_j

        n = size(x)
        allocate (phi_a(n))
        do i = 1, n
            acc = 0.0_dp
            do j = 1, n
                delta = x(i) - x(j)
                k2 = exp(-(delta*delta))*(1.0_dp - 2.0_dp*delta*delta)
                sign_j = sign_of(x(j))
                src = k2*weights(j)*sign_j/sqrt(pi)
                acc = acc + src*alpha(j)
            end do
            phi_a(i) = acc - 2.0_dp/sqrt(pi)*exp(-(x(i)*x(i)))
        end do
    end subroutine aspect_source

    !> Solve the per-angle 2N x 2N complex block system and return the top block.
    subroutine solve_offset_factor(grid, angle, source, alpha)
        type(offset_grid_t), intent(in) :: grid
        real(dp), intent(in) :: angle
        complex(dp), intent(in) :: source(:)
        complex(dp), allocatable, intent(out) :: alpha(:)

        complex(dp), allocatable :: left(:, :), rhs(:)
        integer, allocatable :: ipiv(:)
        complex(dp) :: phase
        real(dp) :: pos, neg
        integer :: n, two_n, i, j, info

        n = grid%points
        two_n = 2*n
        phase = cmplx(cos(2.0_dp*angle), sin(2.0_dp*angle), dp)

        allocate (left(two_n, two_n))
        allocate (rhs(two_n))
        allocate (ipiv(two_n))
        left = (0.0_dp, 0.0_dp)

        do j = 1, n
            pos = positive_of(grid%x(j))
            neg = 1.0_dp - pos
            do i = 1, n
                left(i, j) = -phase*grid%kernel(i, j)*pos
                left(i, n + j) = -grid%kernel(i, j)*neg
                left(n + i, j) = -grid%kernel(i, j)*neg
                left(n + i, n + j) = -conjg(phase)*grid%kernel(i, j)*pos
            end do
            left(j, j) = left(j, j) + (1.0_dp, 0.0_dp)
            left(n + j, n + j) = left(n + j, n + j) + (1.0_dp, 0.0_dp)
        end do

        rhs(1:n) = source
        rhs(n + 1:two_n) = conjg(source)

        call zgesv(two_n, 1, left, two_n, ipiv, rhs, two_n, info)
        if (info /= 0) then
            print *, "offset_factors: angle zgesv failed, info = ", info
            error stop
        end if

        allocate (alpha(n))
        alpha = rhs(1:n)
    end subroutine solve_offset_factor

    !> Sample one offset factor at a single angle. want_b selects S_B, else S_A.
    function sample_factor(grid, angle, want_b) result(value)
        type(offset_grid_t), intent(in) :: grid
        real(dp), intent(in) :: angle
        logical, intent(in) :: want_b
        real(dp) :: value

        complex(dp), allocatable :: source(:), alpha(:)
        complex(dp) :: phase_b

        allocate (source(grid%points))
        if (want_b) then
            phase_b = cmplx(cos(angle), sin(angle), dp)
            source = phase_b*cmplx(grid%phi_b, 0.0_dp, dp)/cmplx(0.0_dp, 2.0_dp, dp)
        else
            source = cmplx(grid%phi_a, 0.0_dp, dp)/cmplx(0.0_dp, 2.0_dp, dp)
        end if

        call solve_offset_factor(grid, angle, source, alpha)
        value = real(alpha(1), dp)
    end function sample_factor

    !> Sign with sign(0) = 0 (matches numpy.sign on the grid node at x = 0).
    pure function sign_of(value) result(s)
        real(dp), intent(in) :: value
        real(dp) :: s

        if (value > 0.0_dp) then
            s = 1.0_dp
        else if (value < 0.0_dp) then
            s = -1.0_dp
        else
            s = 0.0_dp
        end if
    end function sign_of

    !> 1.0 for strictly positive x (x = 0 counts as negative side), else 0.0.
    pure function positive_of(value) result(p)
        real(dp), intent(in) :: value
        real(dp) :: p

        if (value > 0.0_dp) then
            p = 1.0_dp
        else
            p = 0.0_dp
        end if
    end function positive_of

    !> Linear interpolation of alpha at x = 0 (returns the node value on the
    !! odd grid that has an exact x = 0 node).
    pure function interp_zero(x, alpha) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha(:)
        real(dp) :: value

        integer :: i, n
        real(dp) :: frac

        n = size(x)
        if (x(1) >= 0.0_dp) then
            value = alpha(1)
            return
        end if
        if (x(n) <= 0.0_dp) then
            value = alpha(n)
            return
        end if
        do i = 1, n - 1
            if (x(i) <= 0.0_dp .and. x(i + 1) >= 0.0_dp) then
                frac = (0.0_dp - x(i))/(x(i + 1) - x(i))
                value = alpha(i) + frac*(alpha(i + 1) - alpha(i))
                return
            end if
        end do
        value = alpha(n)
    end function interp_zero

    !> Gauss-Legendre nodes and weights on [-1,1] via Newton iteration on the
    !! Legendre polynomial. Deterministic and self-contained.
    subroutine gauss_legendre(order, nodes, weights)
        integer, intent(in) :: order
        real(dp), intent(out) :: nodes(order)
        real(dp), intent(out) :: weights(order)

        integer :: i, iter, mid
        real(dp) :: x0, p, dp_val, dx

        mid = (order + 1)/2
        do i = 1, mid
            x0 = cos(pi*(real(i, dp) - 0.25_dp)/(real(order, dp) + 0.5_dp))
            do iter = 1, 100
                call legendre(order, x0, p, dp_val)
                dx = -p/dp_val
                x0 = x0 + dx
                if (abs(dx) < 1.0e-15_dp) exit
            end do
            call legendre(order, x0, p, dp_val)
            nodes(i) = -x0
            nodes(order + 1 - i) = x0
            weights(i) = 2.0_dp/((1.0_dp - x0*x0)*dp_val*dp_val)
            weights(order + 1 - i) = weights(i)
        end do
    end subroutine gauss_legendre

    !> Legendre polynomial P_n(x) and its derivative via the recurrence.
    pure subroutine legendre(n, x, p, dp_val)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p
        real(dp), intent(out) :: dp_val

        integer :: k
        real(dp) :: p_prev, p_curr, p_next

        p_prev = 1.0_dp
        p_curr = x
        if (n == 0) then
            p = 1.0_dp
            dp_val = 0.0_dp
            return
        end if
        do k = 2, n
            p_next = (real(2*k - 1, dp)*x*p_curr - real(k - 1, dp)*p_prev)/real(k, dp)
            p_prev = p_curr
            p_curr = p_next
        end do
        p = p_curr
        dp_val = real(n, dp)*(x*p_curr - p_prev)/(x*x - 1.0_dp)
    end subroutine legendre

end module offset_factors
