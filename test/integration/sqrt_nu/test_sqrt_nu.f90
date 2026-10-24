program test_sqrt_nu
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sqrt_nu_mod, only: sqrt_nu_result_t, evaluate_sqrt_nu_modes, &
        fe_log_model, compute_A_m
    use sqrt_nu_kernel, only: bounce_time_pointwise
    implicit none
    integer, parameter :: n = 8193
    real(dp) :: eta(n), tau(n), diffusion(n), omega(n), div(n, 2), g(2)
    real(dp) :: exact, error_coarse, error_fine, a, tau_d_boundary, tau_log, omega_model
    real(dp) :: bounce, dbounce_db, dbounce_deta, shifted, ignored1, ignored2
    complex(dp) :: z
    type(sqrt_nu_result_t) :: result, reverse, coarse
    character(32) :: arg
    integer :: j, m

    do j = 1, n
        eta(j) = real(j - 1, dp)/real(n - 1, dp)
    end do
    tau = 1.0_dp
    diffusion = 1.0_dp
    omega = 400.0_dp
    div(:, 1) = -1.0_dp
    div(:, 2) = -2.0_dp
    g = [2.0_dp, 3.0_dp]
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, result)
    if (result%status /= 0) error stop 'constant-coefficient status'
    if (.not. all(result%valid)) error stop 'constant strength/resolution screen'
    do m = 1, 2
        z = sqrt(cmplx(0.0_dp, real(m, dp)*omega(1), dp))
        exact = real(m, dp)*real((1.0_dp - exp(-z))/z, dp)
        if (abs(result%flux_modes(m) - exact) > 1.0e-7_dp) &
            error stop 'constant-coefficient exact solution'
    end do
    exact = sum(0.5_dp*g*result%flux_modes)
    if (abs(result%lambda_off - exact) > 1.0e-14_dp) error stop 'mode normalization'
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, -omega, div, g, 0.5_dp, reverse)
    if (maxval(abs(result%flux_modes - reverse%flux_modes)) > 1.0e-14_dp) &
        error stop 'precession reversal symmetry'
    call evaluate_sqrt_nu_modes(eta(::8), tau(::8), diffusion(::8), omega(::8), &
                               div(::8, :), g, 0.5_dp, coarse)
    z = sqrt(cmplx(0.0_dp, omega(1), dp))
    exact = real((1.0_dp - exp(-z))/z, dp)
    error_coarse = abs(coarse%flux_modes(1) - exact)
    error_fine = abs(result%flux_modes(1) - exact)
    if (error_fine > error_coarse/50.0_dp) error stop 'eta convergence'
    omega = 0.0_dp
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, result)
    if (result%status /= 0) error stop 'E77 status'
    if (.not. result%zero_precession) error stop 'E77 dispatch flag'
    if (any(result%valid)) error stop 'zero precession is outside WKB'
    if (maxval(abs(result%flux_modes - [1.0_dp, 2.0_dp])) > 1.0e-13_dp) &
        error stop 'E77 integrated drift'
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, &
                               result, iv_deep=[0.25_dp, -0.5_dp])
    if (maxval(abs(result%flux_modes - [1.25_dp, 1.5_dp])) > 1.0e-13_dp) &
        error stop 'E77 nonzero deep endpoint'
    omega = 400.0_dp
    omega(n) = -400.0_dp
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, result)
    if (result%status /= 3) error stop 'precession resonance rejected'
    omega = 400.0_dp
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, &
                               result, magnetic_drift_included=.true.)
    if (result%status /= 2) error stop 'magnetic drift rejected'
    omega = 4.0_dp
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, &
                               result, iv_deep=[0.25_dp, -0.5_dp])
    z = sqrt(cmplx(0.0_dp, omega(1), dp))
    exact = real((1.0_dp - exp(-z))/z + 0.25_dp*exp(-z), dp)
    if (abs(result%flux_modes(1) - exact) > 1.0e-8_dp) &
        error stop 'finite precession nonzero deep endpoint'
    omega = 0.01_dp
    call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, 0.5_dp, result)
    if (any(result%valid)) error stop 'weak precession screen'
    ! Recover A_m from FEdef's primitive bounce-time and diffusion prefactors,
    ! independently of the reduced expression used in the generated kernel.
    tau_d_boundary = 4.0_dp*0.001_dp*0.5_dp*3.0_dp*0.3_dp/(acos(-1.0_dp)*10.0_dp)
    tau_log = 4.0_dp*3.0_dp/4.0_dp*sqrt(0.5_dp/(2.0_dp*2.0_dp))
    omega_model = 0.01_dp*4.0_dp/1.0_dp
    exact = 0.5_dp*sqrt(2.0_dp*omega_model*tau_log/tau_d_boundary)
    a = compute_A_m(10.0_dp, 2, 0.01_dp, 0.5_dp, 1.0_dp, 0.001_dp, 0.3_dp, -2.0_dp)
    if (abs(a - exact) > 1.0e-13_dp*exact) error stop 'A_m primitive normalization'
    exact = compute_A_m(10.0_dp, 2, -0.01_dp, 0.5_dp, 1.0_dp, 0.004_dp, &
                        0.3_dp, -2.0_dp)
    if (abs(2.0_dp*exact - a) > 1.0e-13_dp*a) error stop 'A_m sign and nu scaling'
    exact = fe_log_model(100.0_dp, 1.0_dp, 8192)
    if (abs(fe_log_model(100.0_dp, 1.0_dp) - exact) > 1.0e-12_dp) &
        error stop 'log quadrature convergence'
    call bounce_time_pointwise(2.0_dp, 0.1_dp, 3.0_dp, 4.0_dp, &
                               bounce, dbounce_db, dbounce_deta)
    if (abs(bounce*4.0_dp*2.0_dp*sqrt(0.8_dp) - 6.0_dp) > 1.0e-13_dp) &
        error stop 'bounce time normalization'
    call bounce_time_pointwise(2.0_dp + 1.0e-6_dp, 0.1_dp, 3.0_dp, 4.0_dp, &
                               shifted, ignored1, ignored2)
    if (abs((shifted - bounce)/1.0e-6_dp - dbounce_db) > 1.0e-6_dp) &
        error stop 'bounce derivative in field'
    call bounce_time_pointwise(2.0_dp, 0.1_dp + 1.0e-6_dp, 3.0_dp, 4.0_dp, &
                               shifted, ignored1, ignored2)
    if (abs((shifted - bounce)/1.0e-6_dp - dbounce_deta) > 1.0e-5_dp) &
        error stop 'bounce derivative in pitch'
    call get_command_argument(1, arg)
    if (trim(arg) == '--reference-output') then
        do j = 1, 3
            a = 10.0_dp**j
            write(*, '(a,2es25.16)') 'LOG ', a, fe_log_model(a, 1.0_dp)
        end do
        do j = 1, 2
            a = 10.0_dp + real(j - 1, dp)*20.0_dp
            ! tau*D=1-eta gives the regular Kelvin/Bessel equation at eta=1.
            eta = eta*(1.0_dp - 1.0e-9_dp)/eta(n)
            diffusion = 1.0_dp - eta
            omega = a*a
            call evaluate_sqrt_nu_modes(eta, tau, diffusion, omega, div, g, &
                                       0.5_dp, result)
            write(*, '(a,2es25.16)') 'KELVIN ', a, result%flux_modes(1)
        end do
    end if
    print *, 'sqrt(nu) exact, E77, symmetry, convergence and domain checks passed'
end program test_sqrt_nu
