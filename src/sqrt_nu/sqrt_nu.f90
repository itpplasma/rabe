! Leading, decaying WKB solution of nearomnigen.tex, Eqs. (apprsol_s, FEdef).
! This module is deliberately independent of the finite-precession BVP solver.
module sqrt_nu_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use sqrt_nu_kernel, only: wkb_pointwise, log_pointwise, amplitude_pointwise
    implicit none
    private
    public :: sqrt_nu_result_t, evaluate_sqrt_nu_modes
    public :: fe_log_model, fe_log_asymptotic, compute_A_m

    type :: sqrt_nu_result_t
        real(dp), allocatable :: flux_modes(:), lambda_modes(:), fe_modes(:), action(:)
        logical, allocatable :: valid(:)
        real(dp) :: lambda_off = 0.0_dp
        integer :: status = 1
        logical :: zero_precession = .false.
    end type sqrt_nu_result_t
contains
    ! eta increases from the transition boundary to the deeply trapped end.
    ! tau_b [s], tau_d=tau_b*Detaeta [eta**2], omega_theta [s**-1].
    ! dIv_deta is the signed eta derivative. Optional iv_deep supplies its
    ! deeply trapped endpoint value, defaulting to zero for a single well.
    ! normalization includes pi*(dr/dpsi)/integral(dtheta0*dphi/B**2).
    ! valid is only a strength/resolution screen, not a complete physics domain
    ! certificate (near omnigeneity, mild rotation and slowly varying coefficients
    ! must be established by the caller). The action>=10 threshold is diagnostic.
    subroutine evaluate_sqrt_nu_modes(eta, tau_b, tau_d, omega_theta, dIv_deta, &
                                     g_modes, normalization, result, &
                                     magnetic_drift_included, iv_deep)
        real(dp), intent(in) :: eta(:), tau_b(:), tau_d(:), omega_theta(:)
        real(dp), intent(in) :: dIv_deta(:, :), g_modes(:), normalization
        type(sqrt_nu_result_t), intent(out) :: result
        logical, intent(in), optional :: magnetic_drift_included
        real(dp), intent(in), optional :: iv_deep(:)
        integer :: n, nm, j, m
        real(dp) :: nan, root, root_previous, profile, profile_previous
        real(dp) :: derivative, accumulated, step, width, flux, fe, max_step
        logical :: zero_precession

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        result%lambda_off = nan
        if (present(magnetic_drift_included)) then
            if (magnetic_drift_included) then
                result%status = 2
                return
            end if
        end if
        n = size(eta)
        nm = size(g_modes)
        if (n < 2 .or. nm < 1) return
        if (size(tau_b) /= n .or. size(tau_d) /= n) return
        if (size(omega_theta) /= n) return
        if (size(dIv_deta, 1) /= n .or. size(dIv_deta, 2) /= nm) return
        if (present(iv_deep)) then
            if (size(iv_deep) /= nm) return
            if (.not. all(ieee_is_finite(iv_deep))) return
        end if
        if (.not. all(ieee_is_finite(eta))) return
        if (.not. all(ieee_is_finite(tau_b))) return
        if (.not. all(ieee_is_finite(tau_d))) return
        if (.not. all(ieee_is_finite(omega_theta))) return
        if (.not. all(ieee_is_finite(dIv_deta))) return
        if (.not. all(ieee_is_finite(g_modes))) return
        if (.not. ieee_is_finite(normalization)) return
        if (any(eta(2:) <= eta(:n - 1))) return
        if (any(tau_b <= 0.0_dp) .or. any(tau_d <= 0.0_dp)) return
        zero_precession = maxval(abs(omega_theta)) <= 0.0_dp
        result%zero_precession = zero_precession
        if (.not. zero_precession) then
            if (minval(abs(omega_theta)) <= 0.0_dp) then
                result%status = 3
                return
            end if
            if (minval(omega_theta) < 0.0_dp .and. maxval(omega_theta) > 0.0_dp) then
                result%status = 3
                return
            end if
        end if
        allocate(result%flux_modes(nm), result%lambda_modes(nm), &
                 result%fe_modes(nm), result%action(nm), result%valid(nm))
        !$omp parallel do default(none) &
        !$omp shared(n,nm,eta,tau_b,tau_d,omega_theta,dIv_deta,g_modes,normalization) &
        !$omp shared(result,zero_precession,iv_deep) &
        !$omp private(m,j,root,root_previous,profile,profile_previous,derivative) &
        !$omp private(accumulated,step,width,flux,fe,max_step)
        do m = 1, nm
            accumulated = 0.0_dp
            flux = 0.0_dp
            fe = 0.0_dp
            max_step = 0.0_dp
            if (zero_precession) then
                ! Exact E77 branch, deliberately separate from the WKB approximation.
                do j = 2, n
                    width = eta(j) - eta(j - 1)
                    flux = flux - 0.5_dp*width*(dIv_deta(j, m) + dIv_deta(j - 1, m))
                end do
                if (present(iv_deep)) flux = flux + iv_deep(m)
                result%flux_modes(m) = flux
                result%fe_modes(m) = eta(1) - eta(n)
                result%action(m) = 0.0_dp
                result%lambda_modes(m) = normalization*g_modes(m)*flux
                result%valid(m) = .false.
                cycle
            end if
            call wkb_pointwise(real(m, dp), omega_theta(1), tau_b(1), tau_d(1), &
                               accumulated, root_previous, profile_previous, derivative)
            do j = 2, n
                call wkb_pointwise(real(m, dp), omega_theta(j), tau_b(j), tau_d(j), &
                                   accumulated, root, profile, derivative)
                width = eta(j) - eta(j - 1)
                step = 0.5_dp*width*(root_previous + root)
                accumulated = accumulated + step
                max_step = max(max_step, step)
                call wkb_pointwise(real(m, dp), omega_theta(j), tau_b(j), tau_d(j), &
                                   accumulated, root, profile, derivative)
                flux = flux - 0.5_dp*width*(profile*dIv_deta(j, m) + &
                                             profile_previous*dIv_deta(j - 1, m))
                fe = fe - 0.5_dp*width*(profile + profile_previous)
                root_previous = root
                profile_previous = profile
            end do
            if (present(iv_deep)) flux = flux + profile_previous*iv_deep(m)
            result%flux_modes(m) = flux
            result%fe_modes(m) = fe
            result%action(m) = sqrt(2.0_dp)*accumulated
            result%lambda_modes(m) = normalization*g_modes(m)*flux
            result%valid(m) = .not. zero_precession
            if (result%action(m) < 10.0_dp) result%valid(m) = .false.
            if (max_step > 0.25_dp) result%valid(m) = .false.
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(result%lambda_modes))) return
        result%lambda_off = sum(result%lambda_modes)
        result%status = 0
    end subroutine evaluate_sqrt_nu_modes

    ! FEdef with tau_b proportional to -log(x), x=(eta/eta_b)-1, 0<x<1/2.
    ! The logarithmic variable removes the endpoint singularity. Simpson
    ! quadrature is over t=-log(x), with a tail < exp(-tmax)*eta_b.
    function fe_log_model(a, eta_b, n_intervals) result(fe)
        real(dp), intent(in) :: a, eta_b
        integer, intent(in), optional :: n_intervals
        real(dp) :: fe, t, tmax, h, x, phi, integrand, asymptotic, term, total
        integer :: n, j, weight
        fe = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(eta_b)) return
        if (a <= 1.0_dp .or. eta_b <= 0.0_dp) return
        n = 4096
        if (present(n_intervals)) n = n_intervals
        if (n < 2 .or. mod(n, 2) /= 0) return
        tmax = max(50.0_dp, log(a) + 40.0_dp)
        h = (tmax - log(2.0_dp))/real(n, dp)
        total = 0.0_dp
        do j = 0, n
            t = log(2.0_dp) + real(j, dp)*h
            x = exp(-t)
            call log_pointwise(x, a, eta_b, phi, integrand, asymptotic)
            weight = 2
            if (mod(j, 2) == 1) weight = 4
            if (j == 0 .or. j == n) weight = 1
            term = real(weight, dp)*x*integrand
            total = total + term
        end do
        fe = -eta_b*h*total/3.0_dp
    end function fe_log_model

    function fe_log_asymptotic(a, eta_b) result(fe)
        real(dp), intent(in) :: a, eta_b
        real(dp) :: fe, phi, integrand
        fe = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(eta_b)) return
        if (a <= 1.0_dp .or. eta_b <= 0.0_dp) return
        call log_pointwise(0.5_dp, a, eta_b, phi, integrand, fe)
    end function fe_log_asymptotic
    ! Logarithmic-model amplitude, nearomnigen.tex Eq. (largeparam).
    ! All lengths, field and pitch units must be mutually consistent. The sign
    ! of v_E* is immaterial here; the physical Omega convention belongs to input.
    function compute_A_m(major_r, m, v_e_star, eta_b, radius, nu_star, ib, bpp) &
        result(a)
        real(dp), intent(in) :: major_r, v_e_star, eta_b, radius, nu_star, ib, bpp
        integer, intent(in) :: m
        real(dp) :: a
        a = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. all(ieee_is_finite([major_r, v_e_star, eta_b, radius, &
                                     nu_star, ib, bpp]))) return
        if (min(major_r, eta_b, radius, nu_star, ib) <= 0.0_dp) return
        if (m < 1 .or. abs(bpp) <= 0.0_dp) return
        call amplitude_pointwise(major_r, real(m, dp), v_e_star, eta_b, &
                                 radius, nu_star, ib, bpp, a)
    end function compute_A_m
end module sqrt_nu_mod
