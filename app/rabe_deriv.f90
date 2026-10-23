program rabe_deriv
    !! Derivative mode of rabe (opt-in, separate executable): values of
    !! Lambda_A, Lambda_B, Lambda_S, nu_star_crit and eta_b per surface together
    !! with their exact Jacobian with respect to the Boozer spectrum, written to
    !! rabe_deriv.nc. Reads the usual `&rabe_config` from rabe.in and an
    !! optional `&rabe_derivative` group; see src/derivative/README.md. With
    !! `wout_derivative = .true.` (VMEC input only) all spectrum parameters
    !! come from libneo's differentiated VMEC -> Boozer transform, and the
    !! gradients of the outputs with respect to the wout harmonics, `phiedge`
    !! and `Rmajor_p` are written as well.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use constants, only: dp
    use netcdf_mod, only: netcdf_t
    use read_file, only: read_namelist, field_file, field_type, M_pol, N_tor, &
                         s_tor_config => s_tor, sign_sqrtg, max_n_fieldlines
    use boozer_field, only: boozer_field_t
    use spectrum_derivative, only: surface_spectrum_t, derivative_settings_t, &
                                   derivative_structure_t, &
                                   derivative_diagnostics_t, n_outputs, &
                                   n_params, pack_parameters, find_structure, &
                                   evaluate, lambda_off, off_iota, output_names, &
                                   parameter_layout
    use spectrum_io, only: spectralize_field, read_spectrum_file
    use wout_chain, only: wout_chain_t, boozer_wout_tangent_t, wout_chain_init, &
                          wout_chain_spectrum, wout_chain_vjp
    use git_version, only: git_hash

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    integer, parameter :: max_nu = 100
    integer, parameter :: n_diag = 12
    integer, parameter :: grid_refinement = 6

    integer :: mboz, nboz, n_gauss
    logical :: continuum_extrema, wout_derivative
    real(dp) :: nu_star(max_nu)
    character(len=256) :: spectrum_file, output_file
    character(len=8) :: lambda_half_to_full
    namelist /rabe_derivative/ mboz, nboz, n_gauss, nu_star, spectrum_file, &
        output_file, continuum_extrema, wout_derivative, lambda_half_to_full

    type(boozer_field_t) :: field
    type(surface_spectrum_t), allocatable :: spectra(:)
    type(derivative_settings_t) :: settings
    type(derivative_structure_t) :: structure
    type(derivative_diagnostics_t) :: diag
    real(dp), allocatable :: s_tor(:), params(:), values(:, :), jac(:, :, :)
    real(dp), allocatable :: residual(:), asymmetry(:), iota_input(:)
    real(dp), allocatable :: lam_off(:, :), lam_off_jac(:, :, :), nu_list(:)
    real(dp), allocatable :: diag_real(:, :)
    type(wout_chain_t) :: chain
    type(boozer_wout_tangent_t) :: wout_bar
    real(dp), allocatable :: grad_rmnc(:, :, :, :), grad_zmns(:, :, :, :)
    real(dp), allocatable :: grad_lmns(:, :, :, :), grad_iotaf(:, :, :)
    real(dp), allocatable :: grad_phiedge(:, :), grad_rmajor(:, :)
    real(dp) :: unit_cotangent(n_outputs, 1)
    integer :: i_out
    integer, allocatable :: q(:), p_rat(:), n_lines(:), n_split(:)
    integer :: n_surf, K, P, i, i_nu, n_nu
    real(dp) :: iota, B_I, B_G, nan_value
    integer(kind=8) :: t_start, t_end, t_rate

    call read_namelist(input_file)
    call read_derivative_namelist(input_file)
    nan_value = ieee_value(nan_value, ieee_quiet_nan)
    nu_list = pack(nu_star, .not. ieee_is_nan(nu_star))
    n_nu = size(nu_list)

    settings%M_pol = M_pol
    settings%N_tor = N_tor
    settings%max_n_fieldlines = max_n_fieldlines
    settings%n_gauss = n_gauss
    settings%sign_sqrtg = nint(sign_sqrtg)
    settings%continuum_extrema = continuum_extrema

    if (wout_derivative .and. (len_trim(spectrum_file) > 0 .or. &
                               trim(field_type) /= 'vmec_nc')) &
        error stop "rabe_deriv: wout_derivative needs field_type = 'vmec_nc' "// &
        "and no spectrum_file"

    if (trim(lambda_half_to_full) /= 'linear' .and. &
        trim(lambda_half_to_full) /= 'cubic') &
        error stop "rabe_deriv: lambda_half_to_full must be 'linear' or 'cubic'"

    if (len_trim(spectrum_file) > 0) then
        call read_spectrum_file(trim(spectrum_file), s_tor, spectra)
        n_surf = size(s_tor)
        allocate (residual(n_surf), asymmetry(n_surf))
        residual = 0.0_dp
        asymmetry = 0.0_dp
    else if (wout_derivative) then
        ! Every spectrum parameter comes from libneo's differentiated
        ! transform, so rabe's own field and projection are not needed.
        s_tor = s_tor_config
        n_surf = size(s_tor)
        allocate (spectra(n_surf), residual(n_surf), asymmetry(n_surf))
        residual = 0.0_dp
        asymmetry = 0.0_dp
        call set_lambda_map(lambda_half_to_full)
        call wout_chain_init(chain, field_file, mboz, nboz, grid_refinement)
        do i = 1, n_surf
            call wout_chain_spectrum(chain, s_tor(i), spectra(i))
        end do
    else
        s_tor = s_tor_config
        n_surf = size(s_tor)
        allocate (spectra(n_surf), residual(n_surf), asymmetry(n_surf))
        select case (trim(field_type))
        case ('vmec_nc')
            call set_lambda_map(lambda_half_to_full)
            call field%init_from_vmec(field_file, grid_refinement=grid_refinement)
        case ('booz_xform')
            call field%init_from_boozmn(field_file)
        case ('chartmap')
            call field%init_from_chartmap(field_file)
        case default
            error stop "rabe_deriv: unknown field_type"
        end select
        do i = 1, n_surf
            call field%fix_to_surface(s_tor(i))
            call field%get_iota(s_tor(i), iota)
            call field%get_covariant_components(B_I, B_G)
            spectra(i)%iota = iota
            spectra(i)%G = B_G
            spectra(i)%I = B_I
            spectra(i)%psi_edge = field%psi_tor_edge
            spectra(i)%R = field%R
            call spectralize_field(field, field%nfp, mboz, nboz, spectra(i), &
                                   residual(i), asymmetry(i))
        end do
    end if

    K = size(spectra(1)%xm)
    P = n_params(K)
    allocate (values(n_outputs, n_surf), jac(P, n_outputs, n_surf))
    allocate (iota_input(n_surf), q(n_surf), p_rat(n_surf), n_lines(n_surf), &
              n_split(n_surf), diag_real(n_diag, n_surf))
    allocate (lam_off(n_nu, n_surf), lam_off_jac(P, n_nu, n_surf))

    call system_clock(t_start, t_rate)
    do i = 1, n_surf
        if (size(spectra(i)%xm) /= K) error stop "rabe_deriv: mode sets differ"
        call find_structure(spectra(i), settings, structure)
        iota_input(i) = spectra(i)%iota
        spectra(i)%iota = structure%iota_line
        call pack_parameters(spectra(i), params)
        call evaluate(params, spectra(i), settings, structure, values(:, i), &
                      jac(:, :, i), diag)
        q(i) = structure%q
        p_rat(i) = structure%p
        n_lines(i) = structure%n_lines
        n_split(i) = structure%n_split_lines
        diag_real(:, i) = [diag%eta_b_gap, diag%nu_crit_gap, diag%max_newton_shift, &
                           diag%min_lambda_squared, diag%resonant_weight_A, &
                           diag%resonant_weight_B, diag%eta_b_shift, &
                           diag%nu_crit_shift, diag%eta_b_iota_kink, &
                           diag%nu_crit_iota_kink, &
                           merge(1.0_dp, 0.0_dp, diag%continuum_extrema), &
                           merge(1.0_dp, 0.0_dp, diag%nu_crit_crossing)]
        do i_nu = 1, n_nu
            call lambda_off(values(:, i), jac(:, :, i), nu_list(i_nu), &
                            lam_off(i_nu, i), lam_off_jac(:, i_nu, i))
        end do
        print "(A,F8.5,A,5ES14.6)", " s_tor ", s_tor(i), &
            "  L_A L_B L_S nu_crit eta_b ", values(:, i)
        print "(A,I4,A,I4,A,ES9.2,A,2ES9.2,A,ES9.2)", "   q ", q(i), "  lines ", &
            n_lines(i), "  eta_b gap ", diag%eta_b_gap, &
            "  continuum shift eta_b, nu_crit ", diag%eta_b_shift, &
            diag%nu_crit_shift, "  spectral residual ", residual(i)
    end do
    call system_clock(t_end)
    print "(A,F10.3,A,I6)", " derivative mode time [s] ", &
        real(t_end - t_start, dp)/real(t_rate, dp), "   n_params ", P

    if (wout_derivative) then
        call system_clock(t_start)
        do i = 1, n_surf
            do i_out = 1, n_outputs
                unit_cotangent = 0.0_dp
                unit_cotangent(i_out, 1) = 1.0_dp
                call wout_chain_vjp(chain, s_tor(i:i), jac(:, :, i:i), &
                                    unit_cotangent, wout_bar)
                if (.not. allocated(grad_rmnc)) then
                    allocate (grad_rmnc(size(wout_bar%rmnc, 1), &
                                        size(wout_bar%rmnc, 2), n_outputs, n_surf))
                    allocate (grad_zmns, grad_lmns, mold=grad_rmnc)
                    allocate (grad_iotaf(size(wout_bar%iotaf), n_outputs, n_surf))
                    allocate (grad_phiedge(n_outputs, n_surf), grad_rmajor(n_outputs, n_surf))
                end if
                grad_rmnc(:, :, i_out, i) = wout_bar%rmnc
                grad_zmns(:, :, i_out, i) = wout_bar%zmns
                grad_lmns(:, :, i_out, i) = wout_bar%lmns
                grad_iotaf(:, i_out, i) = wout_bar%iotaf
                grad_phiedge(i_out, i) = wout_bar%phiedge
                grad_rmajor(i_out, i) = wout_bar%rmajor
            end do
        end do
        call system_clock(t_end)
        print "(A,F10.3)", " wout gradient time [s] ", &
            real(t_end - t_start, dp)/real(t_rate, dp)
    end if

    call write_output()

contains

    subroutine read_derivative_namelist(filename)
        character(len=*), intent(in) :: filename
        integer :: unit, ios

        mboz = 16
        nboz = 16
        n_gauss = 64
        continuum_extrema = .true.
        wout_derivative = .false.
        lambda_half_to_full = 'linear'
        nu_star = ieee_value(1.0_dp, ieee_quiet_nan)
        spectrum_file = ""
        output_file = "rabe_deriv.nc"
        open (newunit=unit, file=filename, status="old", action="read")
        read (unit, nml=rabe_derivative, iostat=ios)
        close (unit)
        if (ios > 0) then
            print *, "rabe_deriv: error reading &rabe_derivative, iostat = ", ios
            error stop
        end if
    end subroutine read_derivative_namelist

    !> libneo's half -> full mesh map of the VMEC lambda harmonics
    !! ('linear' is libneo's default, 'cubic' is more accurate near the axis).
    subroutine set_lambda_map(option)
        use new_vmec_stuff_mod, only: lambda_map => lambda_half_to_full
        character(len=*), intent(in) :: option

        lambda_map = option
    end subroutine set_lambda_map

    subroutine write_output()
        type(netcdf_t) :: nc
        real(dp), allocatable :: bmnc(:, :), smnc(:, :)
        integer :: j

        allocate (bmnc(K, n_surf), smnc(K, n_surf))
        do j = 1, n_surf
            bmnc(:, j) = spectra(j)%bmnc
            smnc(:, j) = spectra(j)%nabla_s_mnc
        end do

        call nc%create(trim(output_file))
        call nc%add_global_attribute("title", &
                                     "rabe derivative mode: values and Jacobians")
        call nc%add_global_attribute("git_hash", git_hash)
        call nc%add_global_attribute("output_names", output_names)
        call nc%add_global_attribute("parameter_layout", parameter_layout)
        call nc%add_global_attribute("jacobian_convention", &
                                     "python: jacobian[surface, output, param] = "// &
                                     "d values[surface, output] / d param; "// &
                                     "B = sum bmnc cos(xm*theta - xn*phi), xn includes nfp")
        call nc%add_global_attribute("iota_convention", &
                                     "iota is the field-line rotational transform "// &
                                     "(rational approximation p/q of iota_input); its "// &
                                     "column is the frozen-rational surrogate")
        call nc%def_dim("surface", n_surf)
        call nc%def_dim("mode", K)
        call nc%def_dim("param", P)
        call nc%def_dim("output", n_outputs)
        call nc%def_dim("diagnostic", n_diag)
        if (n_nu > 0) call nc%def_dim("nu", n_nu)
        call nc%add_real("nfp")
        call nc%add_real_1d("s_tor", "surface")
        call nc%add_real_1d("xm", "mode")
        call nc%add_real_1d("xn", "mode")
        call nc%add_real_2d("bmnc", "mode", "surface")
        call nc%add_attr("bmnc", "unit", "T")
        call nc%add_real_2d("nabla_s_mnc", "mode", "surface")
        call nc%add_attr("nabla_s_mnc", "unit", "1/m")
        call nc%add_real_1d("iota", "surface")
        call nc%add_attr("iota", "long_name", "field-line iota used (p/q approximation)")
        call nc%add_real_1d("iota_input", "surface")
        call nc%add_real_1d("G", "surface")
        call nc%add_attr("G", "unit", "T*m")
        call nc%add_real_1d("I", "surface")
        call nc%add_attr("I", "unit", "T*m")
        call nc%add_real_1d("psi_edge", "surface")
        call nc%add_attr("psi_edge", "unit", "Wb")
        call nc%add_real_1d("R", "surface")
        call nc%add_attr("R", "unit", "m")
        call nc%add_real_2d("values", "output", "surface")
        call nc%add_attr("values", "long_name", output_names)
        call nc%add_real_3d("jacobian", "param", "output", "surface")
        call nc%add_int_1d("q", "surface")
        call nc%add_int_1d("p", "surface")
        call nc%add_int_1d("n_lines", "surface")
        call nc%add_int_1d("n_split_lines", "surface")
        call nc%add_real_2d("diagnostics", "diagnostic", "surface")
        call nc%add_attr("diagnostics", "long_name", &
                         "eta_b_gap nu_crit_gap max_newton_shift "// &
                         "min_lambda_squared resonant_weight_A resonant_weight_B "// &
                         "eta_b_shift nu_crit_shift eta_b_iota_kink "// &
                         "nu_crit_iota_kink continuum_extrema nu_crit_crossing")
        call nc%add_real_1d("spectral_residual", "surface")
        call nc%add_attr("spectral_residual", "long_name", &
                         "max |B_field - B_spectrum|/B_00 (0 for spectrum_file)")
        call nc%add_real_1d("spectral_asymmetry", "surface")
        if (n_nu > 0) then
            call nc%add_real_1d("nu_star", "nu")
            call nc%add_real_2d("lambda_off", "nu", "surface")
            call nc%add_attr("lambda_off", "long_name", &
                             "Lambda_A/sqrt(nu_star) + Lambda_B/nu_star")
            call nc%add_real_3d("lambda_off_jacobian", "param", "nu", "surface")
        end if
        if (wout_derivative) call define_wout_gradients(nc)

        call nc%write_real("nfp", spectra(1)%nfp)
        call nc%write_real_1d("s_tor", s_tor)
        call nc%write_real_1d("xm", spectra(1)%xm)
        call nc%write_real_1d("xn", spectra(1)%xn)
        call nc%write_real_2d("bmnc", bmnc)
        call nc%write_real_2d("nabla_s_mnc", smnc)
        call nc%write_real_1d("iota", [(spectra(j)%iota, j=1, n_surf)])
        call nc%write_real_1d("iota_input", iota_input)
        call nc%write_real_1d("G", [(spectra(j)%G, j=1, n_surf)])
        call nc%write_real_1d("I", [(spectra(j)%I, j=1, n_surf)])
        call nc%write_real_1d("psi_edge", [(spectra(j)%psi_edge, j=1, n_surf)])
        call nc%write_real_1d("R", [(spectra(j)%R, j=1, n_surf)])
        call nc%write_real_2d("values", values)
        call nc%write_real_3d("jacobian", jac)
        call nc%write_int_1d("q", q)
        call nc%write_int_1d("p", p_rat)
        call nc%write_int_1d("n_lines", n_lines)
        call nc%write_int_1d("n_split_lines", n_split)
        call nc%write_real_2d("diagnostics", diag_real)
        call nc%write_real_1d("spectral_residual", residual)
        call nc%write_real_1d("spectral_asymmetry", asymmetry)
        if (n_nu > 0) then
            call nc%write_real_1d("nu_star", nu_list)
            call nc%write_real_2d("lambda_off", lam_off)
            call nc%write_real_3d("lambda_off_jacobian", lam_off_jac)
        end if
        if (wout_derivative) then
            call nc%write_real_4d("grad_rmnc", grad_rmnc)
            call nc%write_real_4d("grad_zmns", grad_zmns)
            call nc%write_real_4d("grad_lmns", grad_lmns)
            call nc%write_real_3d("grad_iotaf", grad_iotaf)
            call nc%write_real_2d("grad_phiedge", grad_phiedge)
            call nc%write_real_2d("grad_rmajor", grad_rmajor)
        end if
        call nc%close()
    end subroutine write_output

    subroutine define_wout_gradients(nc)
        type(netcdf_t), intent(inout) :: nc

        call nc%add_global_attribute("wout_gradient_convention", &
                                     "python: grad_X[surface, output, radius, mn_mode] = "// &
                                     "d values[surface, output] / d wout X[radius, mn_mode] "// &
                                     "for X = rmnc, zmns, lmns (wout half-mesh lmns), "// &
                                     "grad_iotaf[surface, output, radius] w.r.t. wout iotaf, "// &
                                     "grad_phiedge[surface, output] w.r.t. phi_edge "// &
                                     "(whole phi profile scaled), grad_rmajor[surface, "// &
                                     "output] w.r.t. Rmajor_p; libneo boozer_vmec_ad "// &
                                     "chained with the full spectrum Jacobian")
        call nc%add_global_attribute("lambda_half_to_full", trim(lambda_half_to_full))
        call nc%def_dim("radius", size(grad_rmnc, 2))
        call nc%def_dim("mn_mode", size(grad_rmnc, 1))
        call nc%add_real_4d("grad_rmnc", "mn_mode", "radius", "output", "surface")
        call nc%add_real_4d("grad_zmns", "mn_mode", "radius", "output", "surface")
        call nc%add_real_4d("grad_lmns", "mn_mode", "radius", "output", "surface")
        call nc%add_real_3d("grad_iotaf", "radius", "output", "surface")
        call nc%add_real_2d("grad_phiedge", "output", "surface")
        call nc%add_real_2d("grad_rmajor", "output", "surface")
    end subroutine define_wout_gradients

end program rabe_deriv
