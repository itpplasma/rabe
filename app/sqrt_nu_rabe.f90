program sqrt_nu_rabe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use sqrt_nu_mod, only: sqrt_nu_result_t, evaluate_sqrt_nu_modes
    use netcdf_mod, only: netcdf_t
    use git_version, only: git_hash
    implicit none

    integer :: n_surface = 0, n_eta = 0, n_mode = 0, n_nu = 0, n_vE = 0
    logical :: magnetic_drift_included = .false.
    character(len=1024) :: input_file = 'sqrt_nu.in', output_file = 'sqrt_nu.nc'
    character(len=128) :: eta_units = 'T-1', iv_units = 'caller supplied'
    character(len=1024) :: geometry_source = 'caller supplied tables'
    character(len=128) :: normalization_convention = 'caller supplied'
    character(len=128) :: precession_radius_convention = 'caller supplied'
    real(dp) :: particle_speed, endpoint_cutoff_fraction
    real(dp) :: reference_origin_theta = 0.0_dp, reference_origin_phi = 0.0_dp
    real(dp), allocatable :: radius_dkes(:)
    real(dp), allocatable :: s_tor(:), nu_star(:), v_E_star(:), normalization(:)
    real(dp), allocatable :: eta(:, :), tau_b_seconds(:, :), tau_d_per_nu(:, :)
    real(dp), allocatable :: omega_theta_per_vE(:, :), dIv_deta(:, :, :)
    real(dp), allocatable :: g_modes(:, :, :)
    real(dp), allocatable :: iv_deep(:, :), iv_transition(:, :), iv_eta_b(:, :)
    real(dp), allocatable :: e77_anchor(:), transition_out(:), boundary_out(:)
    real(dp), allocatable :: lambda_off(:), flux(:), lambda_mode(:), fe(:), action(:)
    real(dp), allocatable :: g_out(:)
    integer, allocatable :: surface_index(:), nu_index(:), vE_index(:), status(:)
    integer, allocatable :: zero_precession(:), case_index(:), mode_index(:), valid(:)
    integer :: unit, ios, is, inu, ie, im, icase, sample, ncase, nsample
    real(dp) :: nan_value
    type(sqrt_nu_result_t) :: result
    type(netcdf_t) :: nc
    namelist /sqrt_nu_sizes/ n_surface, n_eta, n_mode, n_nu, n_vE
    namelist /sqrt_nu_config/ s_tor, nu_star, v_E_star, eta, tau_b_seconds, &
        tau_d_per_nu, omega_theta_per_vE, dIv_deta, g_modes, normalization, &
        magnetic_drift_included, eta_units, iv_units, iv_deep, &
        iv_transition, iv_eta_b, &
        geometry_source, normalization_convention, precession_radius_convention, &
        particle_speed, radius_dkes, endpoint_cutoff_fraction, &
        reference_origin_theta, reference_origin_phi

    if (command_argument_count() > 0) call get_command_argument(1, input_file)
    if (command_argument_count() > 1) call get_command_argument(2, output_file)
    open (newunit=unit, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'sqrt_nu_rabe: cannot open input'
    read (unit, nml=sqrt_nu_sizes, iostat=ios)
    if (ios /= 0) error stop 'sqrt_nu_rabe: invalid sqrt_nu_sizes'
    if (min(n_surface, n_mode, n_nu, n_vE) < 1 .or. n_eta < 2) &
        error stop 'sqrt_nu_rabe: positive sizes and n_eta >= 2 required'
    allocate (s_tor(n_surface), nu_star(n_nu), v_E_star(n_vE))
    allocate (eta(n_eta, n_surface), tau_b_seconds(n_eta, n_surface))
    allocate (tau_d_per_nu(n_eta, n_surface), omega_theta_per_vE(n_eta, n_surface))
    allocate (dIv_deta(n_eta, n_mode, n_surface), g_modes(n_mode, n_nu, n_surface))
    allocate (normalization(n_surface))
    allocate (iv_deep(n_mode, n_surface), iv_transition(n_mode, n_surface))
    allocate (iv_eta_b(n_mode, n_surface), radius_dkes(n_surface))
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    s_tor = nan_value
    nu_star = nan_value
    v_E_star = nan_value
    eta = nan_value
    tau_b_seconds = nan_value
    tau_d_per_nu = nan_value
    omega_theta_per_vE = nan_value
    dIv_deta = nan_value
    g_modes = nan_value
    normalization = nan_value
    iv_deep = 0.0_dp
    iv_transition = nan_value
    iv_eta_b = nan_value
    radius_dkes = nan_value
    particle_speed = nan_value
    endpoint_cutoff_fraction = nan_value
    read (unit, nml=sqrt_nu_config, iostat=ios)
    close (unit)
    if (ios /= 0) error stop 'sqrt_nu_rabe: invalid sqrt_nu_config'
    if (magnetic_drift_included) &
        error stop 'sqrt_nu_rabe: magnetic precession is outside the SD2.2 model'
    if (.not. all(ieee_is_finite(s_tor))) error stop 'sqrt_nu_rabe: invalid s_tor'
    if (.not. all(ieee_is_finite(nu_star))) error stop 'sqrt_nu_rabe: invalid nu_star'
    if (any(nu_star <= 0.0_dp)) error stop 'sqrt_nu_rabe: nu_star must be positive'
    if (.not. all(ieee_is_finite(v_E_star))) error stop 'sqrt_nu_rabe: invalid v_E_star'

    ncase = n_surface*n_nu*n_vE
    nsample = ncase*n_mode
    allocate (lambda_off(ncase), surface_index(ncase), nu_index(ncase), vE_index(ncase))
    allocate (status(ncase), zero_precession(ncase))
    allocate (flux(nsample), lambda_mode(nsample), fe(nsample), action(nsample))
    allocate (g_out(nsample), case_index(nsample), mode_index(nsample), valid(nsample))
    allocate (e77_anchor(ncase), transition_out(nsample), boundary_out(nsample))
    icase = 0
    sample = 0
    do is = 1, n_surface
        do inu = 1, n_nu
            do ie = 1, n_vE
                icase = icase + 1
                call evaluate_sqrt_nu_modes(eta(:, is), tau_b_seconds(:, is), &
                                            nu_star(inu)*tau_d_per_nu(:, is), &
                           v_E_star(ie)*omega_theta_per_vE(:, is), dIv_deta(:, :, is), &
                                       g_modes(:, inu, is), normalization(is), result, &
                                            magnetic_drift_included, iv_deep(:, is))
                surface_index(icase) = is
                nu_index(icase) = inu
                vE_index(icase) = ie
                status(icase) = result%status
                e77_anchor(icase) = nan_value
                if (all(ieee_is_finite(iv_eta_b(:, is)))) &
                    e77_anchor(icase) = normalization(is)* &
                                        sum(g_modes(:, inu, is)*iv_eta_b(:, is))
                zero_precession(icase) = merge(1, 0, result%zero_precession)
                lambda_off(icase) = nan_value
                if (result%status == 0) lambda_off(icase) = result%lambda_off
                do im = 1, n_mode
                    sample = sample + 1
                    case_index(sample) = icase
                    mode_index(sample) = im
                    g_out(sample) = g_modes(im, inu, is)
                    transition_out(sample) = iv_transition(im, is)
                    boundary_out(sample) = iv_eta_b(im, is)
                    flux(sample) = nan_value
                    lambda_mode(sample) = nan_value
                    fe(sample) = nan_value
                    action(sample) = nan_value
                    valid(sample) = 0
                    if (result%status /= 0) cycle
                    flux(sample) = result%flux_modes(im)
                    lambda_mode(sample) = result%lambda_modes(im)
                    fe(sample) = result%fe_modes(im)
                    action(sample) = result%action(im)
                    valid(sample) = merge(1, 0, result%valid(im))
                end do
            end do
        end do
    end do

    call nc%create(trim(output_file))
    call nc%add_global_attribute('title', &
                                 'SD2.2 tabulated strong-precession WKB offset')
    call nc%add_global_attribute('git_hash', git_hash)
    call nc%add_global_attribute('geometry_source', trim(geometry_source))
    call nc%add_global_attribute('normalization_radius', trim(normalization_convention))
    call nc%add_global_attribute('precession_radius', &
                                 trim(precession_radius_convention))
    call nc%add_global_attribute('definition', &
                                 'lambda_off=sum_m '// &
                       'normalization*g_m*(gamma_s(deep)*Iv_deep-integral gamma_s*dIv)')
    call nc%add_global_attribute('scope', &
                   'Tabulated geometry and normalization; no finite-Er benchmark claim')
    call nc%add_global_attribute('precession', &
                 'Electric drift only; Omega_theta = v_E_star*omega_theta_per_vE [s-1]')
    call nc%add_global_attribute('Omega_hat_convention', &
                                 'Omega_hat=Omega_theta/v=-v_E_star/r for positive '// &
                                 'v_E_star=c*Er/(v*Bref)')
    call nc%add_global_attribute('magnetic_drift_included', 'false; true is rejected')
    call nc%add_global_attribute('validity', &
                    'valid_wkb tests numerical action and resolution, not the full '// &
                                 'physics domain')
    call nc%def_dim('surface', n_surface)
    call nc%def_dim('collisionality', n_nu)
    call nc%def_dim('electric_field', n_vE)
    call nc%def_dim('case', ncase)
    call nc%def_dim('mode_sample', nsample)
    call real_variable('s_tor', 'surface', '1', 'normalized toroidal flux')
    call real_variable('radius_dkes', 'surface', 'm', 'Supplied precession radius')
    call nc%add_real('particle_speed')
    call nc%add_attr('particle_speed', 'units', 'm s-1')
    call nc%add_real('reference_origin_theta')
    call nc%add_attr('reference_origin_theta', 'units', 'rad')
    call nc%add_real('reference_origin_phi')
    call nc%add_attr('reference_origin_phi', 'units', 'rad')
    call nc%add_real('endpoint_cutoff_fraction')
    call nc%add_attr('endpoint_cutoff_fraction', 'units', '1')
    call real_variable('nu_star', 'collisionality', '1', 'pi*R/(2*mean_free_path)')
    call real_variable('v_E_star', 'electric_field', '1', 'signed c*Er/(v*Bref)')
    call real_variable('lambda_off', 'case', '1', 'sum of lambda_mode over harmonics')
    call real_variable('lambda_E77_anchor', 'case', '1', &
                  'Independent E77 anchor from supplied I_V(eta_b), NaN if unavailable')
    call real_variable('Iv_transition_mode', 'mode_sample', trim(iv_units), &
                       'I_V at first pitch grid point, NaN if unavailable')
    call real_variable('Iv_eta_b_mode', 'mode_sample', trim(iv_units), &
                 'I_V at E77 eta_b; compare with Iv_transition_mode for truncation '// &
                       'sensitivity')
    call real_variable('flux_mode', 'mode_sample', trim(iv_units), 'I_V^m-I_E^m')
    call real_variable('lambda_mode', 'mode_sample', '1', 'per-harmonic offset')
    call real_variable('F_E_mode', 'mode_sample', trim(eta_units), &
                       'negative integral of the WKB sine distribution')
    call real_variable('wkb_action', 'mode_sample', '1', &
                       'integral sqrt(abs(m*Omega_theta*tau_b/tau_d)) deta; not '// &
                       'log-model A_m')
    call real_variable('g_mode', 'mode_sample', 'caller supplied', &
                       'E77 boundary-layer harmonic for this collisionality')
    call integer_variable('surface_index', 'case', 'one-based index into s_tor')
    call integer_variable('nu_index', 'case', 'one-based index into nu_star')
    call integer_variable('vE_index', 'case', 'one-based index into v_E_star')
    call integer_variable('status', 'case', &
                          '0 success; 1 invalid input; 2 magnetic drift; 3 zero or '// &
                          'sign-changing Omega')
    call integer_variable('zero_precession', 'case', '1 exact E77 dispatch; 0 WKB')
    call integer_variable('case_index', 'mode_sample', 'one-based index into case')
    call integer_variable('mode_index', 'mode_sample', 'positive harmonic m')
    call integer_variable('valid_wkb', 'mode_sample', &
                      '1 action/resolution checks pass; 0 outside numerical diagnostic')
    call nc%write_real_1d('s_tor', s_tor)
    call nc%write_real_1d('radius_dkes', radius_dkes)
    call nc%write_real('particle_speed', particle_speed)
    call nc%write_real('reference_origin_theta', reference_origin_theta)
    call nc%write_real('reference_origin_phi', reference_origin_phi)
    call nc%write_real('endpoint_cutoff_fraction', endpoint_cutoff_fraction)
    call nc%write_real_1d('nu_star', nu_star)
    call nc%write_real_1d('v_E_star', v_E_star)
    call nc%write_real_1d('lambda_off', lambda_off)
    call nc%write_real_1d('lambda_E77_anchor', e77_anchor)
    call nc%write_real_1d('Iv_transition_mode', transition_out)
    call nc%write_real_1d('Iv_eta_b_mode', boundary_out)
    call nc%write_real_1d('flux_mode', flux)
    call nc%write_real_1d('lambda_mode', lambda_mode)
    call nc%write_real_1d('F_E_mode', fe)
    call nc%write_real_1d('wkb_action', action)
    call nc%write_real_1d('g_mode', g_out)
    call nc%write_int_1d('surface_index', surface_index)
    call nc%write_int_1d('nu_index', nu_index)
    call nc%write_int_1d('vE_index', vE_index)
    call nc%write_int_1d('status', status)
    call nc%write_int_1d('zero_precession', zero_precession)
    call nc%write_int_1d('case_index', case_index)
    call nc%write_int_1d('mode_index', mode_index)
    call nc%write_int_1d('valid_wkb', valid)
    call nc%close()

contains

    subroutine real_variable(name, dimension, units, description)
        character(len=*), intent(in) :: name, dimension, units, description
        call nc%add_real_1d(name, dimension)
        call nc%add_attr(name, 'units', units)
        call nc%add_attr(name, 'long_name', description)
    end subroutine real_variable

    subroutine integer_variable(name, dimension, description)
        character(len=*), intent(in) :: name, dimension, description
        call nc%add_int_1d(name, dimension)
        call nc%add_attr(name, 'units', '1')
        call nc%add_attr(name, 'long_name', description)
    end subroutine integer_variable
end program sqrt_nu_rabe
