module spectrum_derivative
    !! Exact derivatives of the rabe offset coefficients per flux surface with
    !! respect to the Boozer spectrum (derivative mode, opt-in).
    !!
    !! The surface is given by `surface_spectrum_t` (bmnc, nabla_s_mnc, iota, G, I,
    !! psi_edge, R), and all parameters are packed into one flat vector
    !!
    !!     params = [bmnc(1:K), nabla_s_mnc(1:K), iota, G, I, psi_edge, R]
    !!
    !! with K modes. The outputs are `[Lambda_A, Lambda_B, Lambda_S, nu_star_crit,
    !! eta_b]`, and `jacobian(i_param, i_output) = d output / d param`.
    !!
    !! The primal in derivative mode follows the standard rabe pipeline, with the
    !! changes that make it a smooth function of the parameters:
    !!
    !! * the field is the exact Fourier sum (`spectral_surface_t`, no spline),
    !! * the discrete structure is frozen by `find_structure`: the rational
    !!   approximation p/q of iota, the field-line labels, and which local maxima
    !!   bound each line (taken from the standard bracket search),
    !! * the maxima are polished by Newton on dB/dphi = 0 along the line and are
    !!   not nudged inward, and their derivatives follow from the implicit
    !!   function theorem, while B_max follows from the envelope theorem,
    !! * the four line integrals use a fixed Gauss-Legendre rule in the square-root
    !!   substitution phi = a + (b-a) u^2/2 (left half) and phi = b - (b-a) u^2/2
    !!   (right half), which regularises the moving endpoints, and the tangent is
    !!   the exact derivative of that discrete rule including node motion.
    !!
    !! The iota entry of `params` is the field-line rotational transform
    !! (`iota_line`), which the standard pipeline sets to the rational
    !! approximation of the input iota. Its derivative is taken at frozen labels
    !! and frozen p/q (surrogate, see src/derivative/README.md), and S_A, S_B are
    !! differentiated piecewise (S_A' = 0.3, S_B' = 0) except at the modes that sit
    !! exactly on their jumps, which the primal sets to zero and which therefore
    !! get zero slope.
    use constants, only: dp, pi, machine_eps
    use spectral_surface, only: spectral_surface_t, spectral_surface_init

    implicit none
    private

    integer, parameter, public :: n_outputs = 5
    integer, parameter, public :: i_Lambda_A = 1, i_Lambda_B = 2, &
                                  i_Lambda_S = 3, i_nu_star_crit = 4, &
                                  i_eta_b = 5
    character(len=*), parameter, public :: output_names = &
                                           "Lambda_A Lambda_B Lambda_S nu_star_crit eta_b"
    character(len=*), parameter, public :: parameter_layout = &
                                           "bmnc(1:K) nabla_s_mnc(1:K) iota G I psi_edge R"

    !> Offsets of the scalar parameters after the two mode blocks.
    integer, parameter, public :: off_iota = 1, off_G = 2, off_I = 3, &
                                  off_psi_edge = 4, off_R = 5

    real(dp), parameter :: S_A_slope = 0.3_dp
    real(dp), parameter :: fit_eps = 1e-15_dp
        !! must match `eps` in fit_functions
    real(dp), parameter :: c_Lambda_S = 1.5_dp*0.855_dp
    real(dp), parameter :: resonance_tol = 1e-9_dp
    real(dp), parameter :: tie_tol = 1e-10_dp
    !> Floor of the relative symmetry violation used to identify mirror maxima
    !! on one side of a line (`pick_maximum_on_each_side`), well above the
    !! roundoff of B (a few machine_eps times sum |bmnc|/B_00) and well below
    !! physical differences between distinct maxima.
    real(dp), parameter :: pick_roundoff_floor = 5e-13_dp
    !> 1 - eta_b*B below zero by less than this many machine_eps times
    !! eta_b*sum |bmnc| is roundoff of the Fourier sum, not a field point above
    !! the global maximum.
    real(dp), parameter :: L2_roundoff_factor = 64.0_dp
    !> Line maxima within this relative distance of the lattice extremum start
    !! the Newton search for the continuum extremum. It must exceed the
    !! sampling deficit of the lattice, O(B''/B dxi^2).
    real(dp), parameter :: multistart_window = 1e-2_dp
        !! relative tolerance for ties of the discrete max (eta_b) and min
        !! (nu_star_crit) over field-line maxima

    type, public :: surface_spectrum_t
        real(dp) :: nfp = 1.0_dp
        real(dp), allocatable :: xm(:)
            !! poloidal mode numbers
        real(dp), allocatable :: xn(:)
            !! toroidal mode numbers including nfp, B ~ cos(xm*theta - xn*phi)
        real(dp), allocatable :: bmnc(:)
            !! B cosine coefficients in T
        real(dp), allocatable :: nabla_s_mnc(:)
            !! |grad s| cosine coefficients in 1/m
        real(dp) :: iota = 0.0_dp
        real(dp) :: G = 0.0_dp
            !! covariant toroidal component B_phi in T*m
        real(dp) :: I = 0.0_dp
            !! covariant poloidal component B_theta in T*m
        real(dp) :: psi_edge = 0.0_dp
            !! toroidal flux at the edge in Wb (sign as in rabe)
        real(dp) :: R = 0.0_dp
            !! major radius in m
    end type surface_spectrum_t

    type, public :: derivative_settings_t
        real(dp) :: M_pol = 1.0_dp
        real(dp) :: N_tor = 0.0_dp
        integer :: max_n_fieldlines = 200
        integer :: n_gauss = 64
        integer :: sign_sqrtg = 1
        logical :: continuum_extrema = .true.
            !! eta_b and nu_star_crit from the global maximum and the ridge
            !! saddle of B (continuum limit of the discrete max and min over
            !! labels, smooth in iota) instead of the discrete selections
    end type derivative_settings_t

    !> Discrete choices frozen for the derivative mode.
    type, public :: derivative_structure_t
        integer :: p = 0, q = 0, n_lines = 0
        real(dp) :: iota_input = 0.0_dp
        real(dp) :: iota_line = 0.0_dp
            !! rational approximation used as field-line iota
        real(dp), allocatable :: xi(:)
        real(dp), allocatable :: phi_max(:, :)
            !! (2, n_lines) bracket estimates of the bounding maxima
        integer :: n_split_lines = 0
            !! lines on which more than two maxima were found
        real(dp) :: symmetry_violation = 0.0_dp
        real(dp) :: theta_origin = 0.0_dp, phi_origin = 0.0_dp
            !! angular origin of the field-line labels in the spectrum's
            !! angles; half a period when the given origin is a maximum (QI)
    end type derivative_structure_t

    !> Margins and residuals that tell whether the frozen structure is valid.
    type, public :: derivative_diagnostics_t
        integer :: argmax_line = 0, argmax_side = 0, argmin_line = 0
        real(dp) :: eta_b_gap = 0.0_dp
            !! relative gap between global B_max and the next distinct maximum
        real(dp) :: nu_crit_gap = 0.0_dp
            !! relative gap between min B_max(1) and the next distinct one
        integer :: eta_b_ties = 0, nu_crit_ties = 0
            !! sizes of the tie sets of the discrete max and min
        logical :: continuum_extrema = .false.
            !! whether eta_b and nu_star_crit use the continuum extrema
        logical :: nu_crit_crossing = .false.
            !! the continuum min of B_max(1) is a crossing of two maxima on
            !! one line (kink of the ridge), not a saddle of B
        real(dp) :: eta_b_shift = 0.0_dp, nu_crit_shift = 0.0_dp
            !! relative change of B_max and min B_max(1) from the lattice values
        real(dp) :: eta_b_iota_kink = 0.0_dp, nu_crit_iota_kink = 0.0_dp
            !! spread of d B_max / d iota over the tie set divided by B, i.e.
            !! the jump between one-sided iota derivatives (0 if smooth)
        real(dp) :: max_newton_shift = 0.0_dp
        real(dp) :: min_lambda_squared = 0.0_dp
        real(dp) :: resonant_weight_A = 0.0_dp
        real(dp) :: resonant_weight_B = 0.0_dp
            !! share of |rd_sin*delta_cos| on modes with m*iota_p = 0 mod pi
    end type derivative_diagnostics_t

    !> Maximum of B along one line, followed in the label by `ridge_saddle`.
    type :: branch_point_t
        real(dp) :: xi = 0.0_dp, phi = 0.0_dp, B = 0.0_dp
        real(dp) :: f = 0.0_dp
            !! dF/dxi = (N_tor - iota*M_pol)/nfp*B_t
        real(dp) :: slope = 0.0_dp
            !! dphi/dxi of the maximum
        real(dp) :: g_iota = 0.0_dp
            !! dF/d iota at fixed label, B_t*(phi - phi_0)
        real(dp), allocatable :: cos_k(:)
            !! dF/d bmnc at fixed label
        logical :: ok = .false.
    end type branch_point_t

    public :: n_params, pack_parameters, unpack_parameters
    public :: find_structure, evaluate, spectrum_jvp, spectrum_vjp
    public :: gauss_legendre_01, label_dft, lambda_off

contains

    pure integer function n_params(n_modes)
        integer, intent(in) :: n_modes

        n_params = 2*n_modes + 5
    end function n_params

    subroutine pack_parameters(spec, params)
        type(surface_spectrum_t), intent(in) :: spec
        real(dp), allocatable, intent(out) :: params(:)

        integer :: K

        K = size(spec%xm)
        allocate (params(n_params(K)))
        params(1:K) = spec%bmnc
        params(K + 1:2*K) = spec%nabla_s_mnc
        params(2*K + off_iota) = spec%iota
        params(2*K + off_G) = spec%G
        params(2*K + off_I) = spec%I
        params(2*K + off_psi_edge) = spec%psi_edge
        params(2*K + off_R) = spec%R
    end subroutine pack_parameters

    subroutine unpack_parameters(params, template, spec)
        real(dp), intent(in) :: params(:)
        type(surface_spectrum_t), intent(in) :: template
            !! provides xm, xn, nfp
        type(surface_spectrum_t), intent(out) :: spec

        integer :: K

        K = size(template%xm)
        if (size(params) /= n_params(K)) &
            error stop "unpack_parameters: wrong parameter vector length"
        spec%nfp = template%nfp
        spec%xm = template%xm
        spec%xn = template%xn
        spec%bmnc = params(1:K)
        spec%nabla_s_mnc = params(K + 1:2*K)
        spec%iota = params(2*K + off_iota)
        spec%G = params(2*K + off_G)
        spec%I = params(2*K + off_I)
        spec%psi_edge = params(2*K + off_psi_edge)
        spec%R = params(2*K + off_R)
    end subroutine unpack_parameters

    subroutine make_field(spec, field, structure)
        type(surface_spectrum_t), intent(in) :: spec
        type(spectral_surface_t), intent(out) :: field
        type(derivative_structure_t), intent(in), optional :: structure

        call spectral_surface_init(field, spec%xm, spec%xn, spec%bmnc, &
                                   spec%nabla_s_mnc, spec%I, spec%G)
        if (present(structure)) then
            field%theta_origin = structure%theta_origin
            field%phi_origin = structure%phi_origin
        end if
    end subroutine make_field

    !> Freeze the discrete structure (labels, p/q, bounding maxima) with the same
    !! algorithms as `make_flock_of_fieldlines`.
    subroutine find_structure(spec, settings, structure)
        use diophantine, only: rational_approx
        use fieldline_labels, only: calc_iota_p, get_labels
        use field_checks, only: suspect_omnigenous_origin_not_minimum
        type(surface_spectrum_t), intent(in) :: spec
        type(derivative_settings_t), intent(in) :: settings
        type(derivative_structure_t), intent(out) :: structure

        type(spectral_surface_t) :: field
        real(dp) :: M_pol, N_tor, nfp, iota_p
        integer :: j
        logical :: split

        M_pol = settings%M_pol
        N_tor = settings%N_tor
        nfp = spec%nfp
        call check_settings(M_pol, N_tor, nfp)
        call make_field(spec, field)

        structure%iota_input = spec%iota
        iota_p = calc_iota_p(spec%iota, M_pol, N_tor, nfp)
        call rational_approx(iota_p/(2.0_dp*pi), settings%max_n_fieldlines, &
                             structure%p, structure%q)
        call get_labels(settings%max_n_fieldlines, spec%iota, M_pol, N_tor, nfp, &
                        structure%xi, structure%iota_line)
        structure%n_lines = size(structure%xi)

        structure%symmetry_violation = symmetry_violation(field, &
                                                          structure%iota_line, nfp)
        if (suspect_omnigenous_origin_not_minimum(field, N_tor=N_tor, &
                                                  M_pol=M_pol, &
                                                  retol=structure%symmetry_violation)) then
            ! rabe puts the well centre (a minimum of the ideal omnigenous
            ! field) at the origin. A quasi-isodynamic spectrum in booz_xform
            ! convention has its maximum there, so move the origin to the
            ! other stellarator-symmetric point, half a period in the
            ! helical angle (the half-period shift of the QI .bc files).
            if (N_tor /= 0.0_dp) then
                field%phi_origin = pi/nfp
            else
                field%theta_origin = pi
            end if
            if (suspect_omnigenous_origin_not_minimum(field, N_tor=N_tor, &
                                                      M_pol=M_pol, &
                                                      retol=structure%symmetry_violation)) then
                print *, "warning in find_structure: origin of the field is "// &
                    "suspected not to be a minimum, also half a period away"
                field%phi_origin = 0.0_dp
                field%theta_origin = 0.0_dp
            else
                print *, "find_structure: origin is not a minimum, labels "// &
                    "measured from theta, phi = ", field%theta_origin, field%phi_origin
            end if
            structure%theta_origin = field%theta_origin
            structure%phi_origin = field%phi_origin
        end if

        allocate (structure%phi_max(2, structure%n_lines))
        structure%n_split_lines = 0
        do j = 1, structure%n_lines
            call line_bounding_maxima(field, structure%xi(j), structure%iota_line, &
                                      M_pol, N_tor, nfp, &
                                      structure%symmetry_violation, &
                                      structure%phi_max(:, j), split)
            if (split) structure%n_split_lines = structure%n_split_lines + 1
        end do
    end subroutine find_structure

    subroutine check_settings(M_pol, N_tor, nfp)
        real(dp), intent(in) :: M_pol, N_tor, nfp

        if (abs(M_pol - nint(M_pol)) > 1e-15_dp .or. &
            abs(N_tor - nint(N_tor)) > 1e-15_dp) &
            error stop "spectrum_derivative: M_pol and N_tor must be integers"
        if (nint(N_tor) /= 0 .and. abs(N_tor - nfp) > 1e-15_dp*nfp) &
            error stop "spectrum_derivative: nonzero N_tor must equal nfp"
        if (nint(N_tor) == 0 .and. nint(M_pol) /= 1) &
            error stop "spectrum_derivative: M_pol must be 1 if N_tor = 0"
    end subroutine check_settings

    function symmetry_violation(field, iota, nfp) result(violation)
        use utils, only: linspace
        type(spectral_surface_t), intent(in) :: field
        real(dp), intent(in) :: iota, nfp
        real(dp) :: violation

        integer, parameter :: n_points = 1000
        integer :: idx
        real(dp) :: phi(n_points), B_ref, B_sym

        call linspace(0.0_dp, 2.0_dp*pi/nfp, n_points, phi)
        violation = 0.0_dp
        do idx = 1, n_points
            call field%compute_B_mod(iota*phi(idx), phi(idx), B_ref)
            call field%compute_B_mod(-iota*phi(idx), -phi(idx), B_sym)
            violation = max(violation, abs(B_ref - B_sym)/B_ref)
        end do
    end function symmetry_violation

    !> Gauss-Legendre nodes and weights on [0, 1].
    subroutine gauss_legendre_01(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        integer :: i, j, iter
        real(dp) :: z, z_old, p0, p1, p2, dp_dz

        do i = 1, n
            z = cos(pi*(real(i, dp) - 0.25_dp)/(real(n, dp) + 0.5_dp))
            do iter = 1, 100
                p0 = 1.0_dp
                p1 = z
                do j = 2, n
                    p2 = p1
                    p1 = ((2.0_dp*j - 1.0_dp)*z*p2 - (j - 1.0_dp)*p0)/real(j, dp)
                    p0 = p2
                end do
                ! p1 = P_n(z), p0 = P_{n-1}(z)
                if (n == 1) p0 = 1.0_dp
                dp_dz = real(n, dp)*(z*p1 - p0)/(z*z - 1.0_dp)
                z_old = z
                z = z_old - p1/dp_dz
                if (abs(z - z_old) <= 4.0_dp*machine_eps) exit
            end do
            x(n + 1 - i) = 0.5_dp*(1.0_dp + z)
            w(n + 1 - i) = 1.0_dp/((1.0_dp - z*z)*dp_dz*dp_dz)
        end do
    end subroutine gauss_legendre_01

    !> Discrete Fourier transform over the field-line label, as `real_ft`, with
    !! a forward tangent. `dshift` is the tangent of a common shift of all
    !! labels (x_j -> x_j + shift).
    subroutine label_dft(x, f, f_cos, f_sin, df, dshift, df_cos, df_sin)
        real(dp), intent(in) :: x(:), f(:)
        real(dp), intent(out) :: f_cos(:), f_sin(:)
            !! size(f)/2 + 1
        real(dp), intent(in), optional :: df(:, :)
            !! (np, N)
        real(dp), intent(in), optional :: dshift(:)
        real(dp), intent(out), optional :: df_cos(:, :), df_sin(:, :)
            !! (np, N/2 + 1)

        integer :: N, k, j
        real(dp) :: E_re, E_im, c, s, scale
        logical :: tangent
        real(dp), allocatable :: dE_re(:), dE_im(:)

        N = size(f)
        tangent = present(df)
        if (tangent) then
            allocate (dE_re(size(df, 1)), dE_im(size(df, 1)))
        end if
        do k = 0, N/2
            E_re = 0.0_dp
            E_im = 0.0_dp
            if (tangent) then
                dE_re = 0.0_dp
                dE_im = 0.0_dp
            end if
            do j = 1, N
                c = cos(real(k, dp)*x(j))
                s = sin(real(k, dp)*x(j))
                E_re = E_re + f(j)*c
                E_im = E_im - f(j)*s
                if (tangent) then
                    dE_re = dE_re + df(:, j)*c
                    dE_im = dE_im - df(:, j)*s
                end if
            end do
            if (tangent .and. present(dshift)) then
                ! dE/dshift = -i k E
                dE_re = dE_re + real(k, dp)*E_im*dshift
                dE_im = dE_im - real(k, dp)*E_re*dshift
            end if
            if (k == 0 .or. (mod(N, 2) == 0 .and. k == N/2)) then
                scale = 1.0_dp/real(N, dp)
                f_cos(k + 1) = E_re*scale
                f_sin(k + 1) = 0.0_dp
                if (tangent) then
                    df_cos(:, k + 1) = dE_re*scale
                    df_sin(:, k + 1) = 0.0_dp
                end if
            else
                scale = 2.0_dp/real(N, dp)
                f_cos(k + 1) = E_re*scale
                f_sin(k + 1) = -E_im*scale
                if (tangent) then
                    df_cos(:, k + 1) = dE_re*scale
                    df_sin(:, k + 1) = -dE_im*scale
                end if
            end if
        end do
    end subroutine label_dft

    !> Newton polish of a maximum of B along the line theta = theta_0 +
    !! iota*(phi - phi_0), starting from a bracket estimate.
    !! Without `ok` a non-concave start is a fatal error; with `ok` it, a move
    !! by more than `max_move` or no convergence returns ok = .false.
    subroutine polish_maximum(field, theta_0, phi_0, iota, phi, cos_k, sin_k, ok)
        type(spectral_surface_t), intent(in) :: field
        real(dp), intent(in) :: theta_0, phi_0, iota
        real(dp), intent(inout) :: phi
        real(dp), intent(out) :: cos_k(:), sin_k(:)
        logical, intent(out), optional :: ok

        integer, parameter :: max_iter = 50
        real(dp), parameter :: max_move = 0.5_dp, newton_roundoff = 1e-10_dp
        integer :: iter
        real(dp) :: B, B_t, B_p, B_tt, B_tp, B_pp, g, dg, step, phi_start

        phi_start = phi
        if (present(ok)) ok = .false.
        do iter = 1, max_iter
            call field%basis(theta_0 + iota*(phi - phi_0), phi, cos_k, sin_k)
            call field%B_derivatives(cos_k, sin_k, B, B_t, B_p, B_tt, B_tp, B_pp)
            g = B_p + iota*B_t
            dg = B_pp + 2.0_dp*iota*B_tp + iota**2*B_tt
            if (dg >= 0.0_dp) then
                if (present(ok)) return
                print *, "polish_maximum: d2B/dphi2 >= 0 at phi = ", phi
                error stop
            end if
            step = -g/dg
            phi = phi + step
            if (present(ok)) then
                if (abs(phi - phi_start) > max_move) return
            end if
            if (abs(step) <= 4.0_dp*machine_eps*(1.0_dp + abs(phi))) then
                if (present(ok)) ok = .true.
                exit
            end if
        end do
        ! On flat maxima (small dg) the roundoff of g/dg can stay above the
        ! tolerance; a last step at roundoff level is converged.
        if (present(ok) .and. iter > max_iter) &
            ok = abs(step) <= newton_roundoff*(1.0_dp + abs(phi))
        call field%basis(theta_0 + iota*(phi - phi_0), phi, cos_k, sin_k)
    end subroutine polish_maximum

    !> Lowest point of the ridge of left line maxima near the label xi_c, the
    !! continuum limit of the discrete min over labels of B_max(1).
    !!
    !! F(xi) is B at the left maximum that bounds the line of label xi, picked
    !! as in `find_structure` (bracket search, `line_bounding_maxima`) and
    !! polished by Newton. By the envelope theorem F'(xi) = c*B_t with
    !! c = (N_tor - iota*M_pol)/nfp. F is the larger of the local-maximum
    !! branches on that side of the line, so its minimum is either
    !!
    !! * a smooth minimum on one branch, F' = 0: there B_t = 0 and
    !!   B_p + iota*B_t = 0, a saddle of B (mountain pass), or
    !! * a crossing of two branches, F_a = F_b with F_a' < 0 < F_b': two equal
    !!   maxima on one line, e.g. mirror images at a symmetry label.
    !!
    !! F' is sampled on [xi_c - dxi, xi_c + dxi], the sign change from - to +
    !! closest to xi_c is bracketed and bisected to roundoff (the first steps
    !! with the bracket search, then the larger of the two endpoint branches,
    !! and the final endpoints are checked against the bracket search). Bisection needs
    !! neither the ridge curvature, which is tiny on the almost flat ridges of
    !! quasi-axisymmetric fields, nor smoothness across branch crossings. If
    !! the bisection ends on two branches that are not equal there, the
    !! bounding maximum disappears (a fold) and F jumps, so there is no
    !! continuum minimum and ok = .false.
    !!
    !! Outputs the minimum B and its derivatives with respect to bmnc and to
    !! the field-line iota: at a smooth minimum the envelope cos_k and
    !! B_t*(phi - phi_0) (= 0 there), at a crossing the implicit-function
    !! derivative (f_a*g_b - f_b*g_a)/(f_a - f_b) of F_a(xi*) = F_b(xi*).
    subroutine ridge_saddle(field, M_pol, N_tor, nfp, iota, xi_c, dxi, &
                            violation, B, dB_dbmnc, dB_diota, crossing, ok)
        type(spectral_surface_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp, iota, xi_c, dxi
        real(dp), intent(in) :: violation
        real(dp), intent(out) :: B
        real(dp), intent(out) :: dB_dbmnc(:)
        real(dp), intent(out) :: dB_diota
        logical, intent(out) :: crossing, ok

        integer, parameter :: n_sub = 2, max_iter = 200
        integer, parameter :: n_bracket_search = 10
            !! bisection steps with the full bracket search, before the
            !! bracket (dxi/2**11) is followed by continuation of its ends
        integer, parameter :: n_pts = 2*n_sub + 1
        real(dp), parameter :: same_point_tol = 1e-9_dp
        real(dp), parameter :: crossing_tol = 1e-10_dp
        type(branch_point_t) :: pts(n_pts), a, b_, mid, pa, pb
        real(dp) :: c, xm, step
        integer :: i, ia, iter

        ok = .false.
        crossing = .false.
        B = 0.0_dp
        dB_dbmnc = 0.0_dp
        dB_diota = 0.0_dp
        c = (N_tor - iota*M_pol)/nfp

        do i = 1, n_pts
            call bounding_point(xi_c + real(i - n_sub - 1, dp)*dxi/real(n_sub, dp), &
                                pts(i))
            if (.not. pts(i)%ok) return
        end do
        ! sign change - to + of F' closest to xi_c
        ia = 0
        do i = 1, n_pts - 1
            if (pts(i)%f <= 0.0_dp .and. pts(i + 1)%f >= 0.0_dp) then
                if (ia == 0) then
                    ia = i
                elseif (abs(real(i - n_sub, dp) - 0.5_dp) < &
                        abs(real(ia - n_sub, dp) - 0.5_dp)) then
                    ia = i
                end if
            end if
        end do
        if (ia == 0) return
        a = pts(ia)
        b_ = pts(ia + 1)

        do iter = 1, max_iter
            if (.not. (a%f < 0.0_dp)) then
                b_ = a
                exit
            end if
            if (.not. (b_%f > 0.0_dp)) then
                a = b_
                exit
            end if
            xm = 0.5_dp*(a%xi + b_%xi)
            if (.not. (xm > a%xi .and. xm < b_%xi)) exit
            if (iter <= n_bracket_search) then
                call bounding_point(xm, mid)
                if (.not. mid%ok) return
            else
                ! the larger of the two endpoint branches continued to xm
                call line_point(xm, a%phi + a%slope*(xm - a%xi), pa)
                call line_point(xm, b_%phi + b_%slope*(xm - b_%xi), pb)
                if (pa%ok .and. pb%ok) then
                    mid = pa
                    if (pb%B > pa%B) mid = pb
                elseif (pa%ok) then
                    mid = pa
                elseif (pb%ok) then
                    mid = pb
                else
                    return
                end if
            end if
            if (mid%f < 0.0_dp) then
                a = mid
            else
                b_ = mid
            end if
        end do

        ! the endpoints must still be the bounding maxima of their lines
        call bounding_point(a%xi, pa)
        call bounding_point(b_%xi, pb)
        if (.not. (pa%ok .and. pb%ok)) return
        if (abs(pa%B - a%B) > crossing_tol*a%B .or. &
            abs(pb%B - b_%B) > crossing_tol*b_%B) return
        if (abs(a%phi - b_%phi) <= same_point_tol*(1.0_dp + abs(a%phi))) then
            ! smooth minimum: envelope derivatives at the point with smaller F'
            if (abs(a%f) <= abs(b_%f)) then
                mid = a
            else
                mid = b_
            end if
            B = mid%B
            dB_dbmnc = mid%cos_k
            dB_diota = mid%g_iota
            ok = .true.
            return
        end if
        ! two branches: both continued to the midpoint must be equal there
        xm = 0.5_dp*(a%xi + b_%xi)
        call line_point(xm, a%phi + a%slope*(xm - a%xi), pa)
        call line_point(xm, b_%phi + b_%slope*(xm - b_%xi), pb)
        if (.not. (pa%ok .and. pb%ok)) return
        if (abs(pa%B - pb%B) > crossing_tol*pa%B) return
        if (.not. (pa%f < 0.0_dp .and. pb%f > 0.0_dp)) return
        ! Newton on F_a - F_b = 0 removes the offset of the bisection end,
        ! which follows the tie tolerance of the bracket-search pick
        do iter = 1, 5
            step = -(pa%B - pb%B)/(pa%f - pb%f)
            xm = xm + step
            call line_point(xm, pa%phi + pa%slope*step, pa)
            call line_point(xm, pb%phi + pb%slope*step, pb)
            if (.not. (pa%ok .and. pb%ok)) return
            if (abs(step) <= 4.0_dp*machine_eps*(1.0_dp + abs(xm))) exit
        end do
        if (.not. (pa%f < 0.0_dp .and. pb%f > 0.0_dp)) return
        crossing = .true.
        B = max(pa%B, pb%B)
        dB_dbmnc = (pa%f*pb%cos_k - pb%f*pa%cos_k)/(pa%f - pb%f)
        dB_diota = (pa%f*pb%g_iota - pb%f*pa%g_iota)/(pa%f - pb%f)
        ok = .true.

    contains

        !> The bounding left maximum of the line xi.
        subroutine bounding_point(xi, p)
            real(dp), intent(in) :: xi
            type(branch_point_t), intent(out) :: p

            real(dp) :: phi_bracket(2)
            logical :: split

            call line_bounding_maxima(field, xi, iota, M_pol, N_tor, nfp, &
                                      violation, phi_bracket, split)
            call line_point(xi, phi_bracket(1), p)
        end subroutine bounding_point

        !> Maximum of B along the line xi from phi_guess.
        subroutine line_point(xi, phi_guess, p)
            real(dp), intent(in) :: xi, phi_guess
            type(branch_point_t), intent(out) :: p

            real(dp) :: theta_0, phi_0, dg, B_t, B_p, B_tt, B_tp, B_pp
            real(dp) :: sin_k(size(dB_dbmnc))

            allocate (p%cos_k(size(dB_dbmnc)))
            theta_0 = N_tor*xi/nfp
            phi_0 = M_pol*xi/nfp
            p%xi = xi
            p%phi = phi_guess
            call polish_maximum(field, theta_0, phi_0, iota, p%phi, p%cos_k, &
                                sin_k, p%ok)
            if (.not. p%ok) return
            call field%B_derivatives(p%cos_k, sin_k, p%B, B_t, B_p, B_tt, B_tp, B_pp)
            dg = B_pp + 2.0_dp*iota*B_tp + iota**2*B_tt
            p%f = c*B_t
            ! d theta/d xi = c at fixed phi, so dg/dxi = (B_tp + iota*B_tt)*c
            p%slope = -(B_tp + iota*B_tt)*c/dg
            p%g_iota = B_t*(p%phi - phi_0)
        end subroutine line_point
    end subroutine ridge_saddle

    !> Bracket estimates of the two maxima that bound the line of label xi, as
    !! in `make_flock_of_fieldlines`. Mirror maxima are identified within
    !! max(violation, pick_roundoff_floor).
    subroutine line_bounding_maxima(field, xi, iota_line, M_pol, N_tor, nfp, &
                                    violation, phi_max, split)
        use fieldline_mod, only: fieldline_t
        use make_fieldline, only: maxima_t, find_maxima_along_fieldline, &
                                  pick_maximum_on_each_side
        type(spectral_surface_t), intent(in) :: field
        real(dp), intent(in) :: xi, iota_line, M_pol, N_tor, nfp, violation
        real(dp), intent(out) :: phi_max(2)
        logical, intent(out) :: split

        type(fieldline_t) :: line
        type(maxima_t) :: maxima
        real(dp) :: interval(2), phi_max_error(2)

        line%xi_0 = xi
        line%iota = iota_line
        line%theta_0 = N_tor*xi/nfp
        line%phi_0 = M_pol*xi/nfp
        interval = [-1.5_dp*pi, 1.5_dp*pi]/abs(N_tor - iota_line*M_pol) &
                   + line%phi_0
        call find_maxima_along_fieldline(field, line, interval, maxima)
        split = .false.
        if (maxima%n < 2) then
            print *, "line_bounding_maxima: fewer than two maxima on line xi = ", xi
            error stop
        elseif (maxima%n > 2) then
            ! Maxima equal within the symmetry violation are mirror images
            ! and the one closest to phi_0 bounds the well. The exact Fourier
            ! sum of a symmetric spectrum has zero violation, so without a
            ! roundoff floor the choice between two mirror maxima would follow
            ! their last-bit difference (LP QA, s = 0.5: a 1e-14 spectrum
            ! perturbation moved Lambda_S by 1.4e-4).
            call pick_maximum_on_each_side(maxima, line%phi_0, &
                                           max(violation, pick_roundoff_floor), &
                                           phi_max, phi_max_error)
            split = .true.
        else
            phi_max = maxima%phi(1:2)
        end if
    end subroutine line_bounding_maxima

    !> 2D Newton for a critical point of B(theta, phi): a maximum (kind = +1)
    !! or a saddle (kind = 0), started near a field-line maximum.
    subroutine critical_point(field, theta, phi, kind, B, cos_k, ok)
        type(spectral_surface_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi
        integer, intent(in) :: kind
        real(dp), intent(out) :: B
        real(dp), intent(out) :: cos_k(:)
        logical, intent(out) :: ok

        integer, parameter :: max_iter = 60
        real(dp), parameter :: max_move = 0.5_dp
        real(dp) :: th, ph, B_t, B_p, B_tt, B_tp, B_pp, det, d_th, d_ph
        real(dp) :: sin_k(size(cos_k))
        integer :: iter

        th = theta
        ph = phi
        ok = .false.
        do iter = 1, max_iter
            call field%basis(th, ph, cos_k, sin_k)
            call field%B_derivatives(cos_k, sin_k, B, B_t, B_p, B_tt, B_tp, B_pp)
            det = B_tt*B_pp - B_tp**2
            if (abs(det) <= tiny(1.0_dp)) return
            d_th = -(B_pp*B_t - B_tp*B_p)/det
            d_ph = -(-B_tp*B_t + B_tt*B_p)/det
            th = th + d_th
            ph = ph + d_ph
            if (abs(th - theta) + abs(ph - phi) > max_move) return
            if (abs(d_th) + abs(d_ph) <= 8.0_dp*machine_eps*(1.0_dp + abs(th) + abs(ph))) then
                ok = .true.
                exit
            end if
        end do
        if (.not. ok) return
        call field%basis(th, ph, cos_k, sin_k)
        call field%B_derivatives(cos_k, sin_k, B, B_t, B_p, B_tt, B_tp, B_pp)
        det = B_tt*B_pp - B_tp**2
        if (kind == 1) then
            ok = det > 0.0_dp .and. B_tt < 0.0_dp
        else
            ok = det < 0.0_dp
        end if
    end subroutine critical_point

    !> Smooth primal of the frozen-structure pipeline and, if requested, its
    !! dense forward Jacobian `jacobian(n_params, n_outputs)`.
    subroutine evaluate(params, template, settings, structure, values, &
                        jacobian, diagnostics)
        use fieldline_labels, only: calc_iota_p
        use fit_functions, only: S_A, S_B
        real(dp), intent(in) :: params(:)
        type(surface_spectrum_t), intent(in) :: template
        type(derivative_settings_t), intent(in) :: settings
        type(derivative_structure_t), intent(in) :: structure
        real(dp), intent(out) :: values(n_outputs)
        real(dp), intent(out), optional :: jacobian(:, :)
        type(derivative_diagnostics_t), intent(out), optional :: diagnostics

        type(surface_spectrum_t) :: spec
        type(spectral_surface_t) :: field
        type(derivative_diagnostics_t) :: diag
        integer :: K, np, n, j, side, i, half, ng, i_iota, i_G, i_I, i_psi, i_R
        integer :: n_modes_label, kk
        real(dp) :: M_pol, N_tor, nfp, iota, G, I_cov, psi_edge, R
        real(dp) :: theta_0, phi_0, phi, theta, B, B_t, B_p, B_tt, B_tp, B_pp
        real(dp) :: s_val, s_t, s_p, dg, L, u, h, wu, cA, cB, L2, lam
        real(dp) :: f(4), acc(4), eta_b, B_glob, c0, dxi
        real(dp) :: Gt, S_sum, I_ref, iota_p, diota_p
        real(dp) :: norm, B_sq, lam_b, nab, raw_A, raw_B, dev_A, dev_B, dr
        real(dp) :: Lambda_A, Lambda_B, Lambda_S, hf, I_hat, nu_crit, B1_min
        real(dp) :: delta_max, x_mode, w_res, w_all, slope, factor_A, factor_S
        real(dp) :: tie_lo, tie_hi, B_crit, B_best, dB_diota, dB_diota_best
        real(dp) :: L2_roundoff
        logical :: ok, crossing
        real(dp), allocatable :: cos_c(:), cos_best(:)
        real(dp), allocatable :: d_tie(:)
        real(dp), allocatable :: cos_k(:), sin_k(:), xg(:), wg(:)
        real(dp), allocatable :: phi_m(:, :), B_m(:, :), dphi_m(:, :, :), dB_m(:, :, :)
        real(dp), allocatable :: F_int(:, :), dF_int(:, :, :), dacc(:, :)
        real(dp), allocatable :: d_eta(:), dphi(:), dtheta(:), dB(:), dBt(:), ds(:)
        real(dp), allocatable :: dlam(:), df(:, :), dg_p(:)
        real(dp), allocatable :: I_j(:), dI_j(:, :), d_S(:), dI_ref(:)
        real(dp), allocatable :: delta_A(:), d_delta_A(:, :), delta_eta(:), d_delta_eta(:, :)
        real(dp), allocatable :: rd(:), d_rd(:, :), x_lab(:), x_shift(:), d_shift(:)
        real(dp), allocatable :: rd_c(:), rd_s(:), dA_c(:), dA_s(:), de_c(:), de_s(:)
        real(dp), allocatable :: d_rd_c(:, :), d_rd_s(:, :), d_dA_c(:, :), d_dA_s(:, :)
        real(dp), allocatable :: d_de_c(:, :), d_de_s(:, :)
        real(dp), allocatable :: d_norm(:), d_Bsq(:), d_lamb(:), d_nab(:)
        real(dp), allocatable :: d_rawA(:), d_rawB(:), d_devA(:), d_devB(:), d_dr(:)
        real(dp), allocatable :: d_LA(:), d_LB(:), d_LS(:), d_nu(:), d_hf(:), d_Ihat(:)
        real(dp), allocatable :: d_delta_max(:), e_iota(:), d_Gt(:)

        K = size(template%xm)
        if (size(params) /= n_params(K)) error stop "evaluate: wrong params length"
        np = 0
        if (present(jacobian)) then
            np = n_params(K)
            if (size(jacobian, 1) /= np .or. size(jacobian, 2) /= n_outputs) &
                error stop "evaluate: jacobian must have shape (n_params, n_outputs)"
        end if
        i_iota = 2*K + off_iota
        i_G = 2*K + off_G
        i_I = 2*K + off_I
        i_psi = 2*K + off_psi_edge
        i_R = 2*K + off_R

        call unpack_parameters(params, template, spec)
        call make_field(spec, field, structure)
        M_pol = settings%M_pol
        N_tor = settings%N_tor
        nfp = template%nfp
        iota = spec%iota
        G = spec%G
        I_cov = spec%I
        psi_edge = spec%psi_edge
        R = spec%R
        n = structure%n_lines

        allocate (cos_k(K), sin_k(K))
        allocate (e_iota(np))
        e_iota = 0.0_dp
        if (np > 0) e_iota(i_iota) = 1.0_dp

        ! ---- maxima: Newton polish, implicit function theorem, envelope ----
        allocate (phi_m(2, n), B_m(2, n), dphi_m(np, 2, n), dB_m(np, 2, n))
        allocate (dg_p(np), dtheta(np))
        diag%max_newton_shift = 0.0_dp
        do j = 1, n
            theta_0 = N_tor*structure%xi(j)/nfp
            phi_0 = M_pol*structure%xi(j)/nfp
            do side = 1, 2
                phi = structure%phi_max(side, j)
                call polish_maximum(field, theta_0, phi_0, iota, phi, cos_k, sin_k)
                diag%max_newton_shift = max(diag%max_newton_shift, &
                                            abs(phi - structure%phi_max(side, j)))
                call field%B_derivatives(cos_k, sin_k, B, B_t, B_p, B_tt, B_tp, B_pp)
                phi_m(side, j) = phi
                B_m(side, j) = B
                if (np > 0) then
                    dg = B_pp + 2.0_dp*iota*B_tp + iota**2*B_tt
                    dg_p = 0.0_dp
                    dg_p(1:K) = (template%xn - iota*template%xm)*sin_k
                    dg_p(i_iota) = B_t + (B_tp + iota*B_tt)*(phi - phi_0)
                    dphi_m(:, side, j) = -dg_p/dg
                    dtheta = iota*dphi_m(:, side, j) + (phi - phi_0)*e_iota
                    dB_m(:, side, j) = B_t*dtheta + B_p*dphi_m(:, side, j)
                    dB_m(1:K, side, j) = dB_m(1:K, side, j) + cos_k
                end if
            end do
            if (phi_m(1, j) >= phi_m(2, j)) error stop "evaluate: maxima crossed"
        end do

        ! ---- eta_b from the global maximum ----
        ! Maxima equal to the global one within tie_tol are mirror images or,
        ! at rational iota, the same physical point reached from two labels.
        ! Their derivatives agree except in the iota (surrogate) direction,
        ! where the discrete max has a kink; the mean over the tie set is the
        ! central (Clarke) derivative and the spread is recorded.
        c0 = 1.0_dp - 2.0_dp*machine_eps
        B_glob = maxval(B_m)
        allocate (d_eta(np), d_tie(np))
        d_tie = 0.0_dp
        diag%eta_b_ties = 0
        tie_lo = huge(1.0_dp)
        tie_hi = -huge(1.0_dp)
        do j = 1, n
            do side = 1, 2
                if (B_m(side, j) >= B_glob*(1.0_dp - tie_tol)) then
                    diag%eta_b_ties = diag%eta_b_ties + 1
                    if (B_m(side, j) == B_glob) then
                        diag%argmax_line = j
                        diag%argmax_side = side
                    end if
                    if (np > 0) then
                        d_tie = d_tie + dB_m(:, side, j)
                        tie_lo = min(tie_lo, dB_m(i_iota, side, j))
                        tie_hi = max(tie_hi, dB_m(i_iota, side, j))
                    end if
                end if
            end do
        end do
        diag%eta_b_gap = relative_gap_below(reshape(B_m, [2*n]), B_glob)
        if (np > 0) then
            d_tie = d_tie/real(diag%eta_b_ties, dp)
            diag%eta_b_iota_kink = (tie_hi - tie_lo)/B_glob
        end if
        diag%continuum_extrema = .false.
        allocate (cos_c(K), cos_best(K))
        if (settings%continuum_extrema) then
            ! Continuum limit of the discrete max over labels: the global maximum
            ! of B on the surface, polished by 2D Newton from every line maximum
            ! near the lattice maximum. A single start from the lattice argmax
            ! can stop at a lower local maximum when the lattice samples that one
            ! at its peak (a line through a symmetry point) and the higher one
            ! off its peak; the primal then jumps when the argmax line switches
            ! (LP QA, s = 0.3: 1.3e-6 relative in eta_b).
            B_best = -huge(1.0_dp)
            do j = 1, n
                do side = 1, 2
                    if (B_m(side, j) < B_glob*(1.0_dp - multistart_window)) cycle
                    theta_0 = N_tor*structure%xi(j)/nfp
                    phi_0 = M_pol*structure%xi(j)/nfp
                    phi = phi_m(side, j)
                    call critical_point(field, theta_0 + iota*(phi - phi_0), phi, &
                                        +1, B_crit, cos_c, ok)
                    if (ok .and. B_crit > B_best .and. &
                        B_crit <= B_glob*(1.0_dp + multistart_window)) then
                        B_best = B_crit
                        cos_best = cos_c
                    end if
                end do
            end do
            if (B_best >= B_glob*(1.0_dp - 1e-13_dp)) then
                diag%continuum_extrema = .true.
                diag%eta_b_shift = (B_best - B_glob)/B_glob
                B_glob = B_best
                if (np > 0) then
                    d_tie = 0.0_dp
                    d_tie(1:K) = cos_best
                end if
            end if
        end if
        eta_b = c0/B_glob
        if (np > 0) d_eta = -c0/B_glob**2*d_tie

        ! ---- line integrals by fixed Gauss-Legendre in sqrt substitution ----
        ng = settings%n_gauss
        ! roundoff bound of eta_b*B for the Fourier sum of B
        L2_roundoff = L2_roundoff_factor*machine_eps*eta_b*sum(abs(spec%bmnc))
        allocate (xg(ng), wg(ng))
        call gauss_legendre_01(ng, xg, wg)
        allocate (F_int(4, n), dF_int(np, 4, n), dacc(np, 4))
        allocate (dphi(np), dB(np), dBt(np), ds(np), dlam(np), df(np, 4))
        diag%min_lambda_squared = huge(1.0_dp)
        do j = 1, n
            theta_0 = N_tor*structure%xi(j)/nfp
            phi_0 = M_pol*structure%xi(j)/nfp
            L = phi_m(2, j) - phi_m(1, j)
            acc = 0.0_dp
            if (np > 0) dacc = 0.0_dp
            do half = 1, 2
                do i = 1, ng
                    u = xg(i)
                    h = 0.5_dp*L*u*u
                    if (half == 1) then
                        phi = phi_m(1, j) + h
                        cA = 1.0_dp - 0.5_dp*u*u
                        cB = 0.5_dp*u*u
                    else
                        phi = phi_m(2, j) - h
                        cA = 0.5_dp*u*u
                        cB = 1.0_dp - 0.5_dp*u*u
                    end if
                    theta = theta_0 + iota*(phi - phi_0)
                    call field%basis(theta, phi, cos_k, sin_k)
                    call field%B_derivatives(cos_k, sin_k, B, B_t, B_p, &
                                             B_tt, B_tp, B_pp)
                    call field%nabla_s_derivatives(cos_k, sin_k, s_val, s_t, s_p)
                    L2 = 1.0_dp - eta_b*B
                    diag%min_lambda_squared = min(diag%min_lambda_squared, L2)
                    if (L2 < -L2_roundoff) then
                        print *, "evaluate: 1 - eta_b*B < 0 on line ", j, &
                            " phi = ", phi, " value ", L2, " roundoff ", &
                            L2_roundoff
                        error stop
                    end if
                    ! B <= B_glob, so 1 - eta_b*B >= 1 - c0 holds exactly; a
                    ! smaller value is roundoff of B at nodes next to a maximum
                    ! (n_gauss = 256 on LP QH gave -5.6e-16)
                    L2 = max(L2, 1.0_dp - c0)
                    lam = sqrt(L2)
                    f(1) = 0.5_dp*lam*(3.0_dp + L2)*B_t/B**3
                    f(2) = lam/B**2
                    f(3) = 1.0_dp/B**2
                    f(4) = s_val/B**2
                    wu = wg(i)*u
                    acc = acc + wu*f
                    if (np > 0) then
                        dphi = cA*dphi_m(:, 1, j) + cB*dphi_m(:, 2, j)
                        dtheta = iota*dphi + (phi - phi_0)*e_iota
                        dB = B_t*dtheta + B_p*dphi
                        dB(1:K) = dB(1:K) + cos_k
                        dBt = B_tt*dtheta + B_tp*dphi
                        dBt(1:K) = dBt(1:K) - template%xm*sin_k
                        ds = s_t*dtheta + s_p*dphi
                        ds(K + 1:2*K) = ds(K + 1:2*K) + cos_k
                        dlam = -(eta_b*dB + B*d_eta)/(2.0_dp*lam)
                        df(:, 1) = 0.5_dp*(3.0_dp*(1.0_dp + L2)*B_t/B**3*dlam &
                                           + lam*(3.0_dp + L2)* &
                                           (dBt/B**3 - 3.0_dp*B_t/B**4*dB))
                        df(:, 2) = dlam/B**2 - 2.0_dp*lam/B**3*dB
                        df(:, 3) = -2.0_dp/B**3*dB
                        df(:, 4) = ds/B**2 - 2.0_dp*s_val/B**3*dB
                        dacc = dacc + wu*df
                    end if
                end do
            end do
            F_int(:, j) = L*acc
            if (np > 0) then
                do kk = 1, 4
                    dF_int(:, kk, j) = L*dacc(:, kk) + &
                                       acc(kk)*(dphi_m(:, 2, j) - dphi_m(:, 1, j))
                end do
            end if
        end do

        ! ---- per-line reductions ----
        allocate (I_j(n), dI_j(np, n), d_S(np), dI_ref(np), d_Gt(np))
        Gt = G + iota*I_cov
        if (Gt <= machine_eps) error stop "evaluate: G + iota*I must be positive"
        d_Gt = 0.0_dp
        if (np > 0) then
            d_Gt(i_G) = 1.0_dp
            d_Gt(i_I) = iota
            d_Gt(i_iota) = I_cov
        end if
        I_j = F_int(2, :)*Gt
        do j = 1, n
            dI_j(:, j) = dF_int(:, 2, j)*Gt + F_int(2, j)*d_Gt
        end do
        S_sum = sum(1.0_dp/sqrt(I_j))
        I_ref = (real(n, dp)/S_sum)**2
        d_S = 0.0_dp
        do j = 1, n
            d_S = d_S - 0.5_dp*I_j(j)**(-1.5_dp)*dI_j(:, j)
        end do
        dI_ref = -2.0_dp*real(n, dp)**2/S_sum**3*d_S

        allocate (delta_A(n), d_delta_A(np, n), delta_eta(n), d_delta_eta(np, n))
        allocate (rd(n), d_rd(np, n))
        do j = 1, n
            delta_A(j) = sqrt(I_ref/I_j(j)) - 1.0_dp
            d_delta_A(:, j) = 0.5_dp*sqrt(I_ref/I_j(j))* &
                              (dI_ref/I_ref - dI_j(:, j)/I_j(j))
            delta_eta(j) = 1.0_dp/B_m(1, j) - eta_b
            d_delta_eta(:, j) = -dB_m(:, 1, j)/B_m(1, j)**2 - d_eta
            rd(j) = F_int(1, j)
            d_rd(:, j) = dF_int(:, 1, j)
        end do

        ! ---- Fourier transform over labels ----
        iota_p = calc_iota_p(iota, M_pol, N_tor, nfp)
        diota_p = -sign(pi, iota*M_pol - N_tor)*nfp/(iota*M_pol - N_tor)**2
        n_modes_label = n/2 + 1
        allocate (x_lab(n), x_shift(n), d_shift(np))
        x_lab = structure%xi
        x_shift = structure%xi - iota_p
        d_shift = -diota_p*e_iota
        allocate (rd_c(n_modes_label), rd_s(n_modes_label), dA_c(n_modes_label), &
                  dA_s(n_modes_label), de_c(n_modes_label), de_s(n_modes_label))
        allocate (d_rd_c(np, n_modes_label), d_rd_s(np, n_modes_label), &
                  d_dA_c(np, n_modes_label), d_dA_s(np, n_modes_label), &
                  d_de_c(np, n_modes_label), d_de_s(np, n_modes_label))
        if (np > 0) then
            call label_dft(x_lab, rd, rd_c, rd_s, d_rd, df_cos=d_rd_c, df_sin=d_rd_s)
            call label_dft(x_lab, delta_A, dA_c, dA_s, d_delta_A, &
                           df_cos=d_dA_c, df_sin=d_dA_s)
            call label_dft(x_shift, delta_eta, de_c, de_s, d_delta_eta, d_shift, &
                           d_de_c, d_de_s)
        else
            call label_dft(x_lab, rd, rd_c, rd_s)
            call label_dft(x_lab, delta_A, dA_c, dA_s)
            call label_dft(x_shift, delta_eta, de_c, de_s)
        end if

        ! ---- surface averages ----
        dxi = 2.0_dp*pi/real(n, dp)
        allocate (d_norm(np), d_Bsq(np), d_lamb(np), d_nab(np))
        norm = sum(F_int(3, :))*dxi
        B_sq = sum(phi_m(2, :) - phi_m(1, :))*dxi/norm
        lam_b = sum(F_int(2, :))*dxi/norm
        nab = sum(F_int(4, :))*dxi/norm
        d_norm = sum(dF_int(:, 3, :), dim=2)*dxi
        d_Bsq = sum(dphi_m(:, 2, :) - dphi_m(:, 1, :), dim=2)*dxi/norm - B_sq*d_norm/norm
        d_lamb = sum(dF_int(:, 2, :), dim=2)*dxi/norm - lam_b*d_norm/norm
        d_nab = sum(dF_int(:, 4, :), dim=2)*dxi/norm - nab*d_norm/norm

        ! ---- deviations with frozen-rational surrogate for S_A, S_B ----
        allocate (d_rawA(np), d_rawB(np))
        raw_A = 0.0_dp
        raw_B = 0.0_dp
        d_rawA = 0.0_dp
        d_rawB = 0.0_dp
        w_res = 0.0_dp
        w_all = 0.0_dp
        do kk = 1, n_modes_label
            x_mode = iota_p*real(kk - 1, dp)
            raw_A = raw_A + rd_s(kk)*dA_c(kk)*S_A(x_mode)
            slope = S_A_slope
            if (abs(mod(x_mode, pi)) < fit_eps) slope = 0.0_dp
            d_rawA = d_rawA + (d_rd_s(:, kk)*dA_c(kk) + rd_s(kk)*d_dA_c(:, kk))* &
                     S_A(x_mode) + rd_s(kk)*dA_c(kk)*slope*real(kk - 1, dp)* &
                     diota_p*e_iota
            raw_B = raw_B + rd_s(kk)*de_c(kk)*S_B(x_mode)
            d_rawB = d_rawB + (d_rd_s(:, kk)*de_c(kk) + rd_s(kk)*d_de_c(:, kk))* &
                     S_B(x_mode)
            w_all = w_all + abs(rd_s(kk)*dA_c(kk))
            if (abs(sin(x_mode)) < resonance_tol) &
                w_res = w_res + abs(rd_s(kk)*dA_c(kk))
        end do
        diag%resonant_weight_A = safe_ratio(w_res, w_all)
        w_res = 0.0_dp
        w_all = 0.0_dp
        do kk = 1, n_modes_label
            x_mode = iota_p*real(kk - 1, dp)
            w_all = w_all + abs(rd_s(kk)*de_c(kk))
            if (abs(sin(x_mode)) < resonance_tol) &
                w_res = w_res + abs(rd_s(kk)*de_c(kk))
        end do
        diag%resonant_weight_B = safe_ratio(w_res, w_all)
        raw_A = pi*raw_A
        raw_B = pi*raw_B
        d_rawA = pi*d_rawA
        d_rawB = pi*d_rawB

        allocate (d_devA(np), d_devB(np), d_dr(np))
        factor_A = B_sq/lam_b*sqrt(eta_b)*sqrt(I_ref)/norm
        dev_A = raw_A*factor_A
        d_devA = d_rawA*factor_A + raw_A*factor_A* &
                 (d_Bsq/B_sq - d_lamb/lam_b + 0.5_dp*d_eta/eta_b &
                  + 0.5_dp*dI_ref/I_ref - d_norm/norm)
        dev_B = raw_B*B_sq/lam_b*0.5_dp/norm
        d_devB = d_rawB*B_sq/lam_b*0.5_dp/norm + &
                 dev_B*(d_Bsq/B_sq - d_lamb/lam_b - d_norm/norm)

        ! ---- coefficients ----
        dr = real(settings%sign_sqrtg, dp)/(nab*psi_edge)
        d_dr = -dr*d_nab/nab
        if (np > 0) d_dr(i_psi) = d_dr(i_psi) - dr/psi_edge

        allocate (d_LA(np), d_LB(np), d_LS(np), d_nu(np), d_hf(np), d_Ihat(np))
        allocate (d_delta_max(np))
        Lambda_A = dev_A*dr*sqrt(0.5_dp*R*pi)
        d_LA = (d_devA*dr + dev_A*d_dr)*sqrt(0.5_dp*R*pi)
        if (np > 0) d_LA(i_R) = d_LA(i_R) + 0.5_dp*Lambda_A/R
        Lambda_B = dev_B*0.5_dp*R*pi*dr
        d_LB = (d_devB*dr + dev_B*d_dr)*0.5_dp*R*pi
        if (np > 0) d_LB(i_R) = d_LB(i_R) + Lambda_B/R

        hf = (G*M_pol + I_cov*N_tor)/(M_pol*iota - N_tor)
        d_hf = 0.0_dp
        if (np > 0) then
            d_hf(i_G) = M_pol/(M_pol*iota - N_tor)
            d_hf(i_I) = N_tor/(M_pol*iota - N_tor)
            d_hf(i_iota) = -M_pol*hf/(M_pol*iota - N_tor)
        end if
        I_hat = I_ref/eta_b/(2.0_dp*pi*R)
        d_Ihat = I_hat*(dI_ref/I_ref - d_eta/eta_b)
        if (np > 0) d_Ihat(i_R) = d_Ihat(i_R) - I_hat/R
        factor_S = -c_Lambda_S/sqrt(pi)*sqrt(I_hat)*B_sq*eta_b**2/lam_b
        Lambda_S = hf*dr*factor_S
        d_LS = (d_hf*dr + hf*d_dr)*factor_S + Lambda_S* &
               (0.5_dp*d_Ihat/I_hat + d_Bsq/B_sq + 2.0_dp*d_eta/eta_b - d_lamb/lam_b)

        B1_min = minval(B_m(1, :))
        diag%argmin_line = minloc(B_m(1, :), dim=1)
        diag%nu_crit_gap = relative_gap_above(B_m(1, :), B1_min)
        delta_max = eta_b - 1.0_dp/B1_min
        ! tie set of the discrete min, treated as for eta_b
        d_tie = 0.0_dp
        diag%nu_crit_ties = 0
        tie_lo = huge(1.0_dp)
        tie_hi = -huge(1.0_dp)
        do j = 1, n
            if (B_m(1, j) <= B1_min*(1.0_dp + tie_tol)) then
                diag%nu_crit_ties = diag%nu_crit_ties + 1
                if (np > 0) then
                    d_tie = d_tie + dB_m(:, 1, j)
                    tie_lo = min(tie_lo, dB_m(i_iota, 1, j))
                    tie_hi = max(tie_hi, dB_m(i_iota, 1, j))
                end if
            end if
        end do
        if (np > 0) then
            d_tie = d_tie/real(diag%nu_crit_ties, dp)
            diag%nu_crit_iota_kink = (tie_hi - tie_lo)/B1_min
        end if
        if (diag%continuum_extrema) then
            ! Continuum limit of the discrete min over labels of the left
            ! maxima: the lowest saddle of B on the maximum ridge (mountain
            ! pass), by a bracketed 1D search along the ridge from every
            ! lattice local minimum near the global one (`ridge_saddle`). 2D
            ! Newton on grad B = 0 was used before and converged to nearby
            ! maxima of B on the flat ridge of LP QA, s = 0.7.
            ! Saddles further below than the window are not on the ridge and
            ! are rejected.
            B_best = huge(1.0_dp)
            dB_diota_best = 0.0_dp
            dxi = 2.0_dp*pi/real(n, dp)
            do j = 1, n
                if (B_m(1, j) > B1_min*(1.0_dp + multistart_window)) cycle
                ! lattice local minima only (cyclic in the label)
                if (B_m(1, j) > B_m(1, modulo(j - 2, n) + 1) .or. &
                    B_m(1, j) > B_m(1, modulo(j, n) + 1)) cycle
                call ridge_saddle(field, M_pol, N_tor, nfp, iota, &
                                  structure%xi(j), dxi, &
                                  structure%symmetry_violation, B_crit, cos_c, &
                                  dB_diota, crossing, ok)
                if (ok .and. B_crit < B_best .and. &
                    B_crit >= B1_min*(1.0_dp - multistart_window)) then
                    B_best = B_crit
                    cos_best = cos_c
                    dB_diota_best = dB_diota
                    diag%nu_crit_crossing = crossing
                end if
            end do
            if (B_best <= B1_min*(1.0_dp + 1e-13_dp)) then
                diag%nu_crit_shift = (B_best - B1_min)/B1_min
                B1_min = B_best
                delta_max = eta_b - 1.0_dp/B1_min
                if (np > 0) then
                    d_tie = 0.0_dp
                    d_tie(1:K) = cos_best
                    d_tie(i_iota) = dB_diota_best
                end if
            else
                diag%continuum_extrema = .false.
                print *, "warning: saddle of B not found, nu_star_crit from "// &
                    "the discrete min over labels"
            end if
        end if
        d_delta_max = d_eta + d_tie/B1_min**2
        nu_crit = 0.125_dp*(delta_max/eta_b)**2/I_hat
        d_nu = 0.25_dp*delta_max/(eta_b**2*I_hat)*d_delta_max &
               - nu_crit*(2.0_dp*d_eta/eta_b + d_Ihat/I_hat)

        values(i_Lambda_A) = Lambda_A
        values(i_Lambda_B) = Lambda_B
        values(i_Lambda_S) = Lambda_S
        values(i_nu_star_crit) = nu_crit
        values(i_eta_b) = eta_b
        if (np > 0) then
            jacobian(:, i_Lambda_A) = d_LA
            jacobian(:, i_Lambda_B) = d_LB
            jacobian(:, i_Lambda_S) = d_LS
            jacobian(:, i_nu_star_crit) = d_nu
            jacobian(:, i_eta_b) = d_eta
        end if
        if (present(diagnostics)) diagnostics = diag
    end subroutine evaluate

    pure real(dp) function safe_ratio(a, b)
        real(dp), intent(in) :: a, b

        if (b > 0.0_dp) then
            safe_ratio = a/b
        else
            safe_ratio = 0.0_dp
        end if
    end function safe_ratio

    !> Relative gap between `top` and the largest value that is not equal to
    !! it within 1e-10 (equal values are symmetric copies, which keep the
    !! maximum differentiable).
    pure real(dp) function relative_gap_below(v, top)
        real(dp), intent(in) :: v(:), top

        real(dp) :: tol

        tol = 1e-10_dp*abs(top)
        if (any(v < top - tol)) then
            relative_gap_below = (top - maxval(v, mask=v < top - tol))/abs(top)
        else
            relative_gap_below = 0.0_dp
        end if
    end function relative_gap_below

    pure real(dp) function relative_gap_above(v, bottom)
        real(dp), intent(in) :: v(:), bottom

        real(dp) :: tol

        tol = 1e-10_dp*abs(bottom)
        if (any(v > bottom + tol)) then
            relative_gap_above = (minval(v, mask=v > bottom + tol) - bottom)/abs(bottom)
        else
            relative_gap_above = 0.0_dp
        end if
    end function relative_gap_above

    !> Jacobian-vector product: out_tangent(o) = sum_p J(p, o)*tangent(p),
    !! with J(param, output) = d output / d param.
    subroutine spectrum_jvp(params, template, settings, structure, tangent, &
                            values, out_tangent)
        real(dp), intent(in) :: params(:)
        type(surface_spectrum_t), intent(in) :: template
        type(derivative_settings_t), intent(in) :: settings
        type(derivative_structure_t), intent(in) :: structure
        real(dp), intent(in) :: tangent(:)
        real(dp), intent(out) :: values(n_outputs), out_tangent(n_outputs)

        real(dp), allocatable :: jac(:, :)

        allocate (jac(size(params), n_outputs))
        call evaluate(params, template, settings, structure, values, jac)
        out_tangent = matmul(tangent, jac)
    end subroutine spectrum_jvp

    !> Vector-Jacobian product with respect to the spectrum parameters:
    !! gradient(p) = sum_o cotangent(o) * d output_o / d param_p.
    subroutine spectrum_vjp(params, template, settings, structure, cotangent, &
                            values, gradient)
        real(dp), intent(in) :: params(:)
        type(surface_spectrum_t), intent(in) :: template
        type(derivative_settings_t), intent(in) :: settings
        type(derivative_structure_t), intent(in) :: structure
        real(dp), intent(in) :: cotangent(n_outputs)
        real(dp), intent(out) :: values(n_outputs)
        real(dp), intent(out) :: gradient(:)

        real(dp), allocatable :: jac(:, :)

        allocate (jac(size(params), n_outputs))
        call evaluate(params, template, settings, structure, values, jac)
        gradient = matmul(jac, cotangent)
    end subroutine spectrum_vjp

    !> Above-critical asymptote lambda_off = Lambda_A/sqrt(nu) + Lambda_B/nu and
    !! its gradient with respect to the parameters, given values and jacobian.
    subroutine lambda_off(values, jacobian, nu_star, value, gradient)
        real(dp), intent(in) :: values(n_outputs), jacobian(:, :), nu_star
        real(dp), intent(out) :: value, gradient(:)

        value = values(i_Lambda_A)/sqrt(nu_star) + values(i_Lambda_B)/nu_star
        gradient = jacobian(:, i_Lambda_A)/sqrt(nu_star) + &
                   jacobian(:, i_Lambda_B)/nu_star
    end subroutine lambda_off

end module spectrum_derivative
