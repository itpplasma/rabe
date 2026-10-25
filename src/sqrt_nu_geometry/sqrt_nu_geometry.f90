module sqrt_nu_geometry_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constants, only: pi
    use field_base, only: field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use fieldline_labels, only: fieldline_modes_t, fourier_transform_over_label
    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients, only: calc_gradient_scaling_factor_r_eff
    use fit_functions, only: S_A, S_B
    use fourier, only: real_ft
    use integrate_substituted, only: integrate_1d_substituted
    use find_extrema, only: find_local_minima
    use fieldline_integrands, only: local_radial_drift
    use sqrt_nu_kernel, only: bounce_time_pointwise, e77_boundary_pointwise
    implicit none
    private
    public :: sqrt_nu_geometry_t, compute_sqrt_nu_geometry

    type :: sqrt_nu_geometry_t
        real(dp), allocatable :: eta(:), tau_b(:), tau_d_per_nu(:)
        ! Per-label coefficients for the B2 construction: tau_b_line(eta,line)
        ! and tau_d_per_nu_line(line); their label means are tau_b, tau_d_per_nu.
        real(dp), allocatable :: tau_b_line(:, :), tau_d_per_nu_line(:)
        real(dp), allocatable :: dIv_deta(:, :), iv_deep(:), &
                                 iv_transition(:), iv_eta_b(:)
        real(dp), allocatable :: g_modes(:, :)
        real(dp) :: normalization, eta_b, eta_c, eta_t, endpoint_cutoff
    end type sqrt_nu_geometry_t
contains
    subroutine compute_sqrt_nu_geometry(field, flock, R, psi_tor_edge, sign_sqrtg, &
                                        speed, nu_star, n_eta, n_mode, cutoff, data)
        class(field_t), intent(in) :: field
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: R, psi_tor_edge, speed, nu_star(:), cutoff
        integer, intent(in) :: sign_sqrtg, n_eta, n_mode
        type(sqrt_nu_geometry_t), intent(out) :: data
        type(surface_average_t) :: average
        type(fieldline_modes_t) :: modes
        real(dp), allocatable :: minima(:), b_min(:), phi_min(:), values(:)
        real(dp), allocatable :: iv(:, :), line_iv(:), cosine(:), sine(:), labels(:)
        real(dp) :: cov_theta, cov_phi, cov, dr_dAtheta, I_j, l_c, dg_dlc
        real(dp) :: pitch, lo, hi, value, width, grid_step
        integer :: iline, j, m, inu, nline, imin

        if (abs(flock%M_pol) > 1.e-12_dp .or. &
            abs(flock%N_tor - flock%nfp) > 1.e-12_dp) &
            error stop 'sqrt_nu_geometry: only QI helicity M=0,N=nfp is supported'
        if (n_eta < 3 .or. n_mode < 1) error stop 'sqrt_nu_geometry: invalid grid'
        if (.not. ieee_is_finite(cutoff)) &
            error stop 'sqrt_nu_geometry: nonfinite cutoff'
        if (.not. ieee_is_finite(speed) .or. .not. ieee_is_finite(R)) &
            error stop 'sqrt_nu_geometry: nonfinite speed or R'
        if (.not. all(ieee_is_finite(nu_star))) &
            error stop 'sqrt_nu_geometry: nonfinite collisionality'
        if (cutoff <= 0.0_dp .or. cutoff >= 0.1_dp) &
            error stop 'sqrt_nu_geometry: cutoff must be between 0 and 0.1'
        if (speed <= 0.0_dp .or. R <= 0.0_dp .or. any(nu_star <= 0.0_dp)) &
            error stop 'sqrt_nu_geometry: positive speed,R,nu_star required'
        nline = size(flock%fieldlines)
        if (n_mode > (nline - 1)/2) error stop 'sqrt_nu_geometry: unresolved harmonics'
        allocate (b_min(nline), phi_min(nline), line_iv(nline), labels(nline))
        labels = flock%fieldlines%xi_0
        allocate (cosine(nline/2 + 1), sine(nline/2 + 1), iv(n_eta, n_mode))
        call field%get_covariant_components(cov_theta, cov_phi)
        cov = cov_phi + flock%iota*cov_theta
        if (cov <= 0.0_dp) error stop 'sqrt_nu_geometry: nonpositive covariant factor'
        call calc_surface_averages(flock, average)
        call fourier_transform_over_label(flock, modes)
        dr_dAtheta = calc_gradient_scaling_factor_r_eff(flock, psi_tor_edge, sign_sqrtg)
        data%eta_b = flock%eta_b
        data%eta_c = 1.0_dp/minval(flock%fieldlines%B_max(1))
        do iline = 1, nline
            call find_local_minima(b_along_line, &
                                   flock%fieldlines(iline)%phi_max, minima)
            if (size(minima) /= 1) &
                error stop 'sqrt_nu_geometry: exactly one minimum per well is required'
            allocate (values(size(minima)))
            call b_along_line(minima, values)
            imin = minloc(values, dim=1)
            phi_min(iline) = minima(imin)
            b_min(iline) = values(imin)
            deallocate (values, minima)
        end do
        data%eta_t = 1.0_dp/maxval(b_min)
        if (data%eta_t <= data%eta_c) &
            error stop 'sqrt_nu_geometry: no common trapped well'
        width = data%eta_t - data%eta_c
        data%endpoint_cutoff = cutoff*width
        allocate (data%eta(n_eta), data%tau_b(n_eta), data%tau_d_per_nu(n_eta))
        allocate (data%tau_b_line(n_eta, nline), data%tau_d_per_nu_line(nline))
        allocate (data%dIv_deta(n_eta, n_mode), data%iv_deep(n_mode))
        allocate (data%iv_transition(n_mode), data%iv_eta_b(n_mode))
        allocate (data%g_modes(n_mode, size(nu_star)))
        grid_step = (1.0_dp - 2.0_dp*cutoff)*width/real(n_eta - 1, dp)
        data%eta = [(data%eta_c + cutoff*width + real(j - 1, dp)*grid_step, j=1, n_eta)]
        ! FEdef freezes the regular diffusion product at the boundary. It does
        ! not assert the unconfirmed extrapolation 2*eta*I_j(eta)/l_c away from it.
        I_j = cov*sum(flock%fieldlines%integral_lambda_b_over_B_squared)/real(nline, dp)
        data%tau_d_per_nu = 4.0_dp*data%eta_b*I_j/(pi*R)
        data%tau_d_per_nu_line = 4.0_dp*data%eta_b*cov* &
                                flock%fieldlines%integral_lambda_b_over_B_squared/(pi*R)
        data%tau_b = 0.0_dp
        do j = 1, n_eta
            pitch = data%eta(j)
            do iline = 1, nline
                call trapped_interval(lo, hi)
                call integrate_1d_substituted(bounce_integrand, lo, hi, value)
                data%tau_b_line(j, iline) = value
                data%tau_b(j) = data%tau_b(j) + value/real(nline, dp)
                call integrate_1d_substituted(iv_integrand, lo, hi, line_iv(iline))
            end do
            call real_ft(labels, line_iv, cosine, sine)
            iv(j, :) = sine(2:n_mode + 1)
        end do
        data%iv_transition = iv(1, :)
        data%iv_deep = iv(n_eta, :)
        data%iv_eta_b = modes%radial_drift%sin_coeffs(2:n_mode + 1)
        data%dIv_deta(1, :) = &
            (-3.0_dp*iv(1, :) + 4.0_dp*iv(2, :) - iv(3, :))/(2*grid_step)
        data%dIv_deta(n_eta, :) = (3.0_dp*iv(n_eta, :) - 4.0_dp*iv(n_eta - 1, :) + &
                                   iv(n_eta - 2, :))/(2*grid_step)
        do j = 2, n_eta - 1
            data%dIv_deta(j, :) = (iv(j + 1, :) - iv(j - 1, :))/(2*grid_step)
        end do
        do inu = 1, size(nu_star)
            l_c = pi*R/(2.0_dp*nu_star(inu))
            do m = 1, n_mode
                call e77_boundary_pointwise(data%eta_b, flock%I_ref, l_c, &
                                           modes%delta_aspect_ratio%cos_coeffs(m + 1), &
                                            modes%delta_eta%cos_coeffs(m + 1), &
                                            S_A(flock%iota_p*real(m, dp)), &
                                     S_B(flock%iota_p*real(m, dp)), average%B_squared, &
                                  average%lambda_b, dr_dAtheta, average%normalization, &
                                       data%g_modes(m, inu), data%normalization, dg_dlc)
            end do
        end do
        if (.not. all(ieee_is_finite(data%tau_b))) &
            error stop 'sqrt_nu_geometry: nonfinite bounce integral'
        if (.not. all(ieee_is_finite(data%tau_b_line)) .or. &
            .not. all(ieee_is_finite(data%tau_d_per_nu_line))) &
            error stop 'sqrt_nu_geometry: nonfinite per-line coefficient'
        if (any(data%tau_b_line <= 0.0_dp) .or. &
            any(data%tau_d_per_nu_line <= 0.0_dp)) &
            error stop 'sqrt_nu_geometry: nonpositive per-line coefficient'
    contains
        subroutine b_along_line(phi, b)
            real(dp), intent(in) :: phi(:)
            real(dp), intent(out) :: b(:)
            integer :: k
            do k = 1, size(phi)
                call field%compute_B_mod(flock%fieldlines(iline)%get_theta(phi(k)), &
                                         phi(k), b(k))
            end do
        end subroutine b_along_line

        subroutine trapped_interval(left, right)
            real(dp), intent(out) :: left, right
            real(dp) :: outside, inside, middle, b
            integer :: side, iteration
            do side = 1, 2
                outside = flock%fieldlines(iline)%phi_max(side)
                inside = phi_min(iline)
                do iteration = 1, 48
                    middle = 0.5_dp*(outside + inside)
                    call field%compute_B_mod( &
                        flock%fieldlines(iline)%get_theta(middle), middle, b)
                    if (pitch*b >= 1.0_dp) then
                        outside = middle
                    else
                        inside = middle
                    end if
                end do
                ! Move inside by a bounded roundoff guard. The coordinate
                ! omission is explicit and much smaller than the pitch cutoff.
                inside = inside + sign(1.e-12_dp* &
                                       (flock%fieldlines(iline)%phi_max(2) - &
                           flock%fieldlines(iline)%phi_max(1)), phi_min(iline) - inside)
                if (side == 1) left = inside
                if (side == 2) right = inside
            end do
        end subroutine trapped_interval

        function bounce_integrand(phi) result(integrand)
            real(dp), intent(in) :: phi
            real(dp) :: integrand, b, db, deta
            call field%compute_B_mod(flock%fieldlines(iline)%get_theta(phi), phi, b)
            if (pitch*b >= 1.0_dp) &
                error stop 'sqrt_nu_geometry: invalid turning bracket'
            call bounce_time_pointwise(b, pitch, cov, speed, integrand, db, deta)
        end function bounce_integrand

        function iv_integrand(phi) result(integrand)
            real(dp), intent(in) :: phi
            real(dp) :: integrand
            integrand = local_radial_drift(field, &
                                     flock%fieldlines(iline)%get_theta(phi), phi, pitch)
        end function iv_integrand
    end subroutine compute_sqrt_nu_geometry
end module sqrt_nu_geometry_mod
