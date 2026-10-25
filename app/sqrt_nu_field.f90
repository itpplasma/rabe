program sqrt_nu_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines
    use sqrt_nu_geometry_mod, only: sqrt_nu_geometry_t, compute_sqrt_nu_geometry
    use error_handling, only: set_unsafe_mode
    implicit none
    integer :: n_surface = 0, n_eta = 0, n_mode = 0, n_nu = 0, n_vE = 0
    integer :: max_n_fieldlines = 32, sign_sqrtg = 1, split_maxima
    real(dp) :: M_pol = 0.0_dp, N_tor = -1.0_dp
    real(dp) :: speed = 1.0_dp, endpoint_cutoff = 1.e-4_dp
    real(dp) :: reference_origin_theta = 0.0_dp, reference_origin_phi = 0.0_dp
    logical :: magnetic_drift_included = .false.
    character(len=1024) :: input_file = 'sqrt_nu_field.in', output_file = 'sqrt_nu.in'
    character(len=1024) :: field_file = ''
    character(len=20) :: field_type = 'vmec_nc'
    character(len=128) :: eta_units = 'T-1', iv_units = 'T-2'
    real(dp), allocatable :: s_tor(:), nu_star(:), v_E_star(:), radius_dkes(:)
    real(dp), allocatable :: eta(:, :), tau_b_seconds(:, :), tau_d_per_nu(:, :)
    real(dp), allocatable :: omega_theta_per_vE(:, :)
    real(dp), allocatable :: dIv_deta(:, :, :), g_modes(:, :, :)
    real(dp), allocatable :: normalization(:), iv_deep(:, :)
    real(dp), allocatable :: iv_transition(:, :), iv_eta_b(:, :)
    character(len=1024) :: geometry_source
    character(len=128) :: normalization_convention = &
                          'native RABE effective-radius gradient'
    character(len=128) :: precession_radius_convention = &
                          'supplied DKES radius, not r_eff'
    real(dp) :: particle_speed, endpoint_cutoff_fraction
    real(dp) :: iota, nan_value
    integer :: unit, ios, is
    type(boozer_field_t) :: field
    type(flock_of_fieldlines_t) :: flock
    type(sqrt_nu_geometry_t) :: data
    namelist /sqrt_nu_sizes/ n_surface, n_eta, n_mode, n_nu, n_vE
    namelist /sqrt_nu_field_config/ field_file, field_type, s_tor, nu_star, v_E_star, &
        radius_dkes, speed, M_pol, N_tor, max_n_fieldlines, sign_sqrtg, &
        endpoint_cutoff, reference_origin_theta, reference_origin_phi, &
        magnetic_drift_included
    namelist /sqrt_nu_config/ s_tor, nu_star, v_E_star, eta, &
        tau_b_seconds, tau_d_per_nu, &
        omega_theta_per_vE, dIv_deta, g_modes, normalization, magnetic_drift_included, &
        eta_units, iv_units, iv_deep, iv_transition, iv_eta_b, &
        geometry_source, normalization_convention, precession_radius_convention, &
        particle_speed, radius_dkes, endpoint_cutoff_fraction, &
        reference_origin_theta, reference_origin_phi

    if (command_argument_count() > 0) call get_command_argument(1, input_file)
    if (command_argument_count() > 1) call get_command_argument(2, output_file)
    open (newunit=unit, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'sqrt_nu_field: cannot open input'
    read (unit, nml=sqrt_nu_sizes, iostat=ios)
    if (ios /= 0) error stop 'sqrt_nu_field: invalid sizes'
    if (min(n_surface, n_mode, n_nu, n_vE) < 1 .or. n_eta < 3) &
        error stop 'sqrt_nu_field: invalid sizes'
    allocate (s_tor(n_surface), nu_star(n_nu), v_E_star(n_vE), radius_dkes(n_surface))
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    s_tor = nan_value
    nu_star = nan_value
    v_E_star = nan_value
    radius_dkes = nan_value
    read (unit, nml=sqrt_nu_field_config, iostat=ios)
    close (unit)
    if (ios /= 0) error stop 'sqrt_nu_field: invalid field config'
    if (magnetic_drift_included) error stop 'sqrt_nu_field: magnetic drift unsupported'
    if (.not. all(ieee_is_finite(radius_dkes))) &
        error stop 'sqrt_nu_field: missing DKES radii'
    if (any(radius_dkes <= 0.0_dp)) &
        error stop 'sqrt_nu_field: DKES radii must be positive'
    if (.not. all(ieee_is_finite(s_tor))) error stop 'sqrt_nu_field: missing surfaces'
    if (.not. all(ieee_is_finite(nu_star))) error stop 'sqrt_nu_field: missing nu_star'
    if (.not. all(ieee_is_finite(v_E_star))) &
        error stop 'sqrt_nu_field: missing v_E_star'
    call set_unsafe_mode(.false.)
    select case (trim(field_type))
    case ('vmec_nc')
        call field%init_from_vmec(trim(field_file), grid_refinement=6)
    case ('booz_xform')
        call field%init_from_boozmn(trim(field_file))
    case ('chartmap')
        call field%init_from_chartmap(trim(field_file))
    case default
        error stop 'sqrt_nu_field: unsupported field_type'
    end select
    if (N_tor < 0.0_dp) N_tor = field%nfp
    allocate (eta(n_eta, n_surface), tau_b_seconds(n_eta, n_surface))
    allocate (tau_d_per_nu(n_eta, n_surface), omega_theta_per_vE(n_eta, n_surface))
    allocate (dIv_deta(n_eta, n_mode, n_surface), g_modes(n_mode, n_nu, n_surface))
    allocate (normalization(n_surface), iv_deep(n_mode, n_surface))
    allocate (iv_transition(n_mode, n_surface), iv_eta_b(n_mode, n_surface))
    do is = 1, n_surface
        call field%fix_to_surface(s_tor(is))
        call field%get_iota(s_tor(is), iota)
        call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, field, &
                                      M_pol, N_tor, field%nfp, split_maxima, &
                        reference_origin=[reference_origin_theta, reference_origin_phi])
        if (split_maxima /= 0) &
            error stop 'sqrt_nu_field: split maxima outside supported well'
        call compute_sqrt_nu_geometry(field, flock, field%R, field%psi_tor_edge, &
                       sign_sqrtg, speed, nu_star, n_eta, n_mode, endpoint_cutoff, data)
        eta(:, is) = data%eta
        tau_b_seconds(:, is) = data%tau_b
        tau_d_per_nu(:, is) = data%tau_d_per_nu
        omega_theta_per_vE(:, is) = -speed/radius_dkes(is)
        dIv_deta(:, :, is) = data%dIv_deta
        g_modes(:, :, is) = data%g_modes
        normalization(is) = data%normalization
        iv_deep(:, is) = data%iv_deep
        iv_transition(:, is) = data%iv_transition
        iv_eta_b(:, is) = data%iv_eta_b
        print *, 's_tor, eta_b, eta_c, eta_t, cutoff: ', s_tor(is), data%eta_b, &
            data%eta_c, data%eta_t, data%endpoint_cutoff
    end do
    geometry_source = trim(field_file)
    particle_speed = speed
    endpoint_cutoff_fraction = endpoint_cutoff
    open (newunit=unit, file=trim(output_file), status='replace', action='write')
    write (unit, nml=sqrt_nu_sizes)
    write (unit, nml=sqrt_nu_config)
    close (unit)
end program sqrt_nu_field
