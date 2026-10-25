program test_sqrt_nu_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use netcdf_mod, only: netcdf_t
    implicit none
    integer, parameter :: ns = 2, ne = 1001, nm = 2, nn = 2, nv = 3
    integer :: n_surface = ns, n_eta = ne, n_mode = nm, n_nu = nn, n_vE = nv
    integer, parameter :: ncase = ns*nn*nv, nsample = ncase*nm
    real(dp) :: s_tor(ns), nu_star(nn), v_E_star(nv), normalization(ns)
    real(dp) :: eta(ne, ns), tau_b_seconds(ne, ns)
    real(dp) :: tau_d_per_nu(ne, ns), omega_theta_per_vE(ne, ns)
    real(dp) :: dIv_deta(ne, nm, ns), g_modes(nm, nn, ns)
    real(dp) :: flux(nsample), offset(ncase), lambda_mode(nsample), action(nsample)
    integer :: case_index(nsample), mode_index(nsample), zero_precession(ncase)
    integer :: surface_index(ncase), nu_index(ncase), vE_index(ncase), status(ncase)
    integer :: valid(nsample), unit, i, is, inu, ie, im, k, icase, exit_status
    real(dp) :: expected, expected_sum, q, decay, length
    logical :: magnetic_drift_included
    character(len=2048) :: executable
    character(len=256) :: convention, units
    type(netcdf_t) :: nc
    namelist /sqrt_nu_sizes/ n_surface, n_eta, n_mode, n_nu, n_vE
    namelist /sqrt_nu_config/ s_tor, nu_star, v_E_star, eta, tau_b_seconds, &
        tau_d_per_nu, omega_theta_per_vE, dIv_deta, g_modes, normalization, &
        magnetic_drift_included

    call get_command_argument(1, executable)
    s_tor = [0.2_dp, 0.7_dp]
    nu_star = [0.001_dp, 0.002_dp]
    v_E_star = [0.0_dp, -1.0_dp, 1.0_dp]
    normalization = [2.0_dp, 3.0_dp]
    do is = 1, n_surface
        do i = 1, n_eta
            eta(i, is) = 1.0_dp + real(i - 1, dp)/real(n_eta - 1, dp)
        end do
        tau_b_seconds(:, is) = real(is, dp)
        tau_d_per_nu(:, is) = 1.0_dp
        omega_theta_per_vE(:, is) = -1.0_dp
        do im = 1, n_mode
            dIv_deta(:, im, is) = real(im + is, dp)
            do inu = 1, n_nu
                g_modes(im, inu, is) = real(im, dp)/nu_star(inu)
            end do
        end do
    end do
    magnetic_drift_included = .false.
    open (newunit=unit, file='driver.in', status='replace')
    write (unit, nml=sqrt_nu_sizes)
    write (unit, nml=sqrt_nu_config)
    close (unit)
    call execute_command_line('"'//trim(executable)//'" driver.in driver.nc', &
                              exitstat=exit_status)
    if (exit_status /= 0) error stop 'synthetic driver run failed'
    call nc%open('driver.nc')
    call nc%read_real_1d('flux_mode', flux)
    call nc%read_real_1d('lambda_off', offset)
    call nc%read_real_1d('lambda_mode', lambda_mode)
    call nc%read_real_1d('wkb_action', action)
    call nc%read_int_1d('case_index', case_index)
    call nc%read_int_1d('mode_index', mode_index)
    call nc%read_int_1d('surface_index', surface_index)
    call nc%read_int_1d('nu_index', nu_index)
    call nc%read_int_1d('vE_index', vE_index)
    call nc%read_int_1d('status', status)
    call nc%read_int_1d('valid_wkb', valid)
    call nc%read_int_1d('zero_precession', zero_precession)
    call nc%read_global_attribute('Omega_hat_convention', convention)
    call nc%read_attr('F_E_mode', 'units', units)
    if (trim(units) /= 'T-1') error stop 'pitch units missing'
    call nc%read_attr('lambda_mode', 'units', units)
    if (trim(units) /= '1') error stop 'bootstrap units missing'
    call nc%close()
    if (index(convention, '=-v_E_star/r') == 0) error stop 'missing sign convention'
    if (any(status /= 0)) error stop 'unexpected numerical status'
    k = 0
    icase = 0
    do is = 1, n_surface
        do inu = 1, n_nu
            do ie = 1, n_vE
                icase = icase + 1
                if (surface_index(icase) /= is) error stop 'wrong surface axis'
                if (nu_index(icase) /= inu) error stop 'wrong collision axis'
                if (vE_index(icase) /= ie) error stop 'wrong electric axis'
                expected_sum = 0.0_dp
                do im = 1, n_mode
                    k = k + 1
                    if (case_index(k) /= icase) error stop 'wrong case mapping'
                    if (mode_index(k) /= im) error stop 'wrong harmonic mapping'
                    length = eta(n_eta, is) - eta(1, is)
                    if (ie == 1) then
                        expected = -dIv_deta(1, im, is)*length
                        if (zero_precession(icase) /= 1) error stop 'zero branch lost'
                        if (valid(k) /= 0) error stop 'zero incorrectly labelled WKB'
                    else
                        q = sqrt(real(im*is, dp)/(2.0_dp*nu_star(inu)))
                        decay = exp(-q*length)
                        expected = -dIv_deta(1, im, is)*(1.0_dp - decay* &
                                             (cos(q*length) - sin(q*length)))/(2.0_dp*q)
                        if (zero_precession(icase) /= 0) error stop 'wrong zero branch'
                        if (valid(k) /= 1) &
                            error stop 'resolved strong action unlabelled'
                        call assert_close(action(k), q*sqrt(2.0_dp)*length, 1.e-10_dp)
                    end if
                    call assert_close(flux(k), expected, 8.e-4_dp)
                    expected = normalization(is)*g_modes(im, inu, is)*expected
                    call assert_close(lambda_mode(k), expected, 8.e-4_dp)
                    expected_sum = expected_sum + expected
                end do
                call assert_close(offset(icase), expected_sum, 8.e-4_dp)
            end do
        end do
    end do

    call check_per_label_driver()
    magnetic_drift_included = .true.
    open (newunit=unit, file='magnetic.in', status='replace')
    write (unit, nml=sqrt_nu_sizes)
    write (unit, nml=sqrt_nu_config)
    close (unit)
    call execute_command_line('"'//trim(executable)// &
                  '" magnetic.in magnetic.nc > magnetic.log 2>&1', exitstat=exit_status)
    if (exit_status == 0) error stop 'magnetic precession was silently accepted'
contains
    ! B2 round trip: identical labels reproduce lambda_off, spread labels do not,
    ! and B1 (lambda_off) is unchanged by the optional per-label groups.
    subroutine check_per_label_driver()
        integer, parameter :: nl = 3
        integer :: n_line, st(ncase, nl), pass, ic, il
        real(dp) :: tau_b_line(ne, nl, ns), tau_d_per_nu_line(nl, ns)
        real(dp) :: per_label(ncase), line(ncase, nl), plain(ncase), c(nl)
        character(len=32) :: coefficient_construction
        character(len=256) :: attr
        namelist /sqrt_nu_per_label_size/ n_line
        namelist /sqrt_nu_per_label/ tau_b_line, tau_d_per_nu_line, &
            coefficient_construction

        n_line = nl
        coefficient_construction = 'per_label'
        plain = offset
        do pass = 1, 2
            c = merge(0.0_dp, [0.2_dp, -0.2_dp, 0.0_dp], pass == 1)
            do il = 1, nl
                tau_b_line(:, il, :) = tau_b_seconds*(1.0_dp + c(il))
                tau_d_per_nu_line(il, :) = 1.0_dp
            end do
            open (newunit=unit, file='b2.in', status='replace')
            write (unit, nml=sqrt_nu_sizes)
            write (unit, nml=sqrt_nu_config)
            write (unit, nml=sqrt_nu_per_label_size)
            write (unit, nml=sqrt_nu_per_label)
            close (unit)
            call execute_command_line('"'//trim(executable)//'" b2.in b2.nc', &
                                      exitstat=exit_status)
            if (exit_status /= 0) error stop 'per-label driver run failed'
            call nc%open('b2.nc')
            call nc%read_real_1d('lambda_off', offset)
            call nc%read_real_1d('lambda_off_per_label', per_label)
            call nc%read_real_2d('lambda_off_line', line)
            call nc%read_int_2d('status_line', st)
            call nc%read_global_attribute('coefficient_construction', attr)
            call nc%close()
            if (trim(attr) /= 'per_label') error stop 'construction attribute missing'
            if (any(st /= 0)) error stop 'per-label status'
            if (maxval(abs(offset - plain)) > 1.e-14_dp*maxval(abs(plain))) &
                error stop 'B1 lambda_off changed by per-label input'
            do ic = 1, ncase
                if (abs(per_label(ic) - sum(line(ic, :))/real(nl, dp)) > &
                    1.e-14_dp*abs(per_label(ic))) error stop 'per-label mean'
                if (pass == 1) then
                  if (maxval(abs(line(ic, :) - plain(ic))) > 1.e-12_dp*abs(plain(ic))) &
                        error stop 'identical labels must reproduce B1'
                else if (zero_precession(ic) == 0 .and. &
                         abs(per_label(ic) - plain(ic)) < 1.e-8_dp*abs(plain(ic))) then
                    error stop 'spread labels did not change B2'
                end if
            end do
        end do
    end subroutine check_per_label_driver

    subroutine assert_close(actual, reference, relative_tolerance)
        real(dp), intent(in) :: actual, reference, relative_tolerance
        if (.not. ieee_is_finite(actual)) error stop 'non-finite driver result'
        if (abs(actual - reference) > &
            relative_tolerance*max(abs(reference), 1.e-12_dp)) then
            print *, 'actual, reference: ', actual, reference
            error stop 'driver result does not match independent exponential integral'
        end if
    end subroutine assert_close
end program test_sqrt_nu_driver
