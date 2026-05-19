program rabe
    use constants, only: dp, pi
    use utils, only: linspace
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_labels, only: get_labels
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients, only: calc_nu_star_crit
    use coefficients, only: calc_finite_boundary_layer_correction
    use shaing_callen_mod, only: calc_trapped_fraction
    use shaing_callen_mod, only: get_non_omnigenous_remainder
    use netcdf_mod, only: netcdf_t
    use git_version, only: git_hash

    use read_file, only: read_namelist
    use read_file, only: field_file, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         sign_sqrtg, &
                         max_n_fieldlines, &
                         should_calc_shaing_callen, &
                         n_eta

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "rabe.nc"
    character(len=*), parameter :: dat_file = "rabe.dat"

    type(boozer_field_t) :: field

    integer :: n_stor
    integer :: this

    real(dp) :: R ![m]
    type(surface_average_t) :: average
    real(dp) :: dr_dAtheta ![rad/Tm/]
    real(dp) :: iota, approx_iota
    real(dp) :: B_theta_covariant, B_phi_covariant
    real(dp) :: nfp

    real(dp), dimension(:), allocatable :: xi_0

    integer :: n_fieldlines
    type(fieldline_t), dimension(:), allocatable :: fieldlines
    logical :: too_strong_violation

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp), dimension(:), allocatable :: Lambda_A, Lambda_B
    real(dp), dimension(:), allocatable :: nu_star_crit
    real(dp), dimension(:), allocatable :: Lambda_S
    integer, dimension(:), allocatable :: err_flag
    real(dp) :: helical_factor

    real(dp) :: trapped_fraction
    real(dp), dimension(:), allocatable :: lambda_SC, remainder

    type(netcdf_t) :: nc_output
    character(len=*), parameter :: dim_name = "surface"
    character(len=1024) :: description

    call read_namelist(input_file)

    n_stor = size(s_tor)
    allocate (Lambda_A(n_stor))
    allocate (Lambda_B(n_stor))
    allocate (nu_star_crit(n_stor))
    allocate (Lambda_S(n_stor))
    allocate (err_flag(n_stor))
    if (should_calc_shaing_callen) then
        allocate (lambda_SC(n_stor))
        allocate (remainder(n_stor))
    end if

    call field%boozer_field_init(field_file, grid_refinement=6)
    do this = 1, n_stor
        call field%fix_to_surface(s_tor(this))
        call field%get_iota_and_covariant_components(s_tor(this), &
                                                     iota, &
                                                     B_theta_covariant, &
                                                     B_phi_covariant)
        nfp = field%nfp

        call get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfp, &
                        xi_0, approx_iota)
        n_fieldlines = size(xi_0)
        allocate (fieldlines(n_fieldlines))
        iota = approx_iota

        call make_flock_of_fieldlines(fieldlines, &
                                      xi_0, &
                                      iota, &
                                      field, &
                                      M_pol, &
                                      N_tor, &
                                      nfp, &
                                      err_flag(this))

        call calc_deviation(fieldlines, deviation_A, deviation_B)

        covariant_factor = (B_phi_covariant + B_theta_covariant*iota)
        call calc_surface_averages(fieldlines, average)
        dr_dAtheta = sign_sqrtg/(average%nabla_s*field%psi_tor_edge)
        R = field%R
        Lambda_A(this) = deviation_A*dr_dAtheta* &
                         sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
        Lambda_B(this) = deviation_B*0.5*R*pi*dr_dAtheta
        nu_star_crit(this) = calc_nu_star_crit(fieldlines, &
                                               R, &
                                               B_theta_covariant, &
                                               B_phi_covariant)
        Lambda_S(this) = calc_finite_boundary_layer_correction(fieldlines, &
                                                               R, &
                                                               dr_dAtheta, &
                                                               B_theta_covariant, &
                                                               B_phi_covariant)

        print *, "s_tor: ", s_tor(this)
        if (should_calc_shaing_callen) then
            trapped_fraction = calc_trapped_fraction(field, fieldlines, n_eta)
            helical_factor = (B_phi_covariant*M_pol + &
                              B_theta_covariant*N_tor)/(M_pol*iota - N_tor)
            lambda_SC(this) = helical_factor*trapped_fraction
            lambda_SC(this) = lambda_SC(this)*dr_dAtheta
            remainder(this) = get_non_omnigenous_remainder(field, fieldlines, n_eta)
            remainder(this) = remainder(this)*covariant_factor*dr_dAtheta* &
                              nfp/(M_pol*iota - N_tor)
            print *, "omnigenous lambda_SC_bB: ", lambda_SC(this)
            print *, "non-omnigneous remainder: ", remainder(this)
        end if
        print *, "1/sqrt(nu_star) factor: ", Lambda_A(this)
        print *, "1/nu_star factor: ", Lambda_B(this)
        print *, "Lambda_S: ", Lambda_S(this)

        if (allocated(fieldlines)) deallocate (fieldlines)
        if (allocated(xi_0)) deallocate (xi_0)

    end do

    call nc_output%create(output_file)
    call nc_output%add_global_attribute("title", &
                                        "asymptotic bootstrap coefficient lambda_bB")
    call nc_output%add_global_attribute("definition", &
                          "lambda^{off}_bB = Lambda_A/sqrt(nu_star) + Lambda_B/nu_star")
    call nc_output%add_global_attribute("git_hash", git_hash)
    call nc_output%def_dim(dim_name, n_stor)
    call nc_output%add_real_1d("s_tor", dim_name)
    call nc_output%add_attr("s_tor", "long_name", &
                            "normalized toroidal flux label")
    call nc_output%add_attr("s_tor", "unit", &
                            "[1]")
    call nc_output%add_real_1d("Lambda_A", dim_name)
    call nc_output%add_attr("Lambda_A", "long_name", &
                            "1/sqrt(nu_star) factor")
    call nc_output%add_attr("Lambda_A", "unit", &
                            "[1]")
    call nc_output%add_real_1d("Lambda_B", dim_name)
    call nc_output%add_attr("Lambda_B", "long_name", &
                            "1/nu_star factor")
    call nc_output%add_attr("Lambda_B", "unit", &
                            "[1]")
    call nc_output%add_real_1d("nu_star_crit", dim_name)
    call nc_output%add_attr("nu_star_crit", "long_name", &
                          "lower collisionality limit for validity of asymptotic model")
    call nc_output%add_attr("nu_star_crit", "unit", &
                            "[1]")
    call nc_output%add_real_1d("Lambda_S", dim_name)
    call nc_output%add_attr("Lambda_S", "long_name", &
                      "sqrt(nu_star) factor accounting for finite boundary layer width")
    call nc_output%add_attr("Lambda_S", "unit", &
                            "[1]")
    call nc_output%add_int_1d("err_flag", dim_name)
    call nc_output%add_attr("err_flag", "long_name", &
                            "1 if violation of omnigeneity is too strong, 0 otherwise")
    call nc_output%add_real("R")
    call nc_output%add_attr("R", "long_name", &
                            "major radius")
    call nc_output%add_attr("R", "unit", "[m]")
    write (description, "(A,A,A)") "defines the reference length scale to convert ", &
        "to dimensionless quantities i.e. defines coeffients in respect to ", &
        "nu_star = pi*R/(2*mean_free_path) = pi*R*deflection_frequency/particle_speed"
    call nc_output%add_attr("R", "definition", description)
    if (should_calc_shaing_callen) then
        call nc_output%add_real_1d("lambda_SC_bB", dim_name)
        call nc_output%add_attr("lambda_SC_bB", "long_name", &
                                "omnigenous Shaing-Callen coefficient")
        call nc_output%add_attr("lambda_SC_bB", "unit", "[1]")
        call nc_output%add_real_1d("remainder", dim_name)
        call nc_output%add_attr("remainder", "long_name", &
                                "non-omnigenous remainder of Shaing-Callen coefficient")
        call nc_output%add_attr("remainder", "unit", "[1]")

        call nc_output%write_real_1d("lambda_SC_bB", lambda_SC)
        call nc_output%write_real_1d("remainder", remainder)
    end if
    call nc_output%write_real_1d("Lambda_A", Lambda_A)
    call nc_output%write_real_1d("Lambda_B", Lambda_B)
    call nc_output%write_real_1d("nu_star_crit", nu_star_crit)
    call nc_output%write_real_1d("Lambda_S", Lambda_S)
    call nc_output%write_int_1d("err_flag", err_flag)
    call nc_output%write_real_1d("s_tor", s_tor)
    call nc_output%write_real("R", field%R)
    call nc_output%close()

    if (should_calc_shaing_callen) then
        call write_dat_output(dat_file, git_hash, field%R, n_stor, &
                              s_tor, Lambda_A, Lambda_B, &
                              nu_star_crit, Lambda_S, err_flag, &
                              lambda_SC=lambda_SC, remainder=remainder)
    else
        call write_dat_output(dat_file, git_hash, field%R, n_stor, &
                              s_tor, Lambda_A, Lambda_B, &
                              nu_star_crit, Lambda_S, err_flag)
    end if

    if (allocated(Lambda_A)) deallocate (Lambda_A)
    if (allocated(Lambda_B)) deallocate (Lambda_B)
    if (allocated(nu_star_crit)) deallocate (nu_star_crit)
    if (allocated(Lambda_S)) deallocate (Lambda_S)
    if (allocated(lambda_SC)) deallocate (lambda_SC)
    if (allocated(remainder)) deallocate (remainder)

contains

    subroutine write_dat_output(filename, git_hash, R, n, s_tor, &
                                Lambda_A, Lambda_B, nu_star_crit, &
                                Lambda_S, err_flag, &
                                lambda_SC, remainder)
        use constants, only: dp
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: git_hash
        real(dp), intent(in) :: R
        integer, intent(in) :: n
        real(dp), intent(in) :: s_tor(n), Lambda_A(n), Lambda_B(n)
        real(dp), intent(in) :: nu_star_crit(n), Lambda_S(n)
        integer, intent(in) :: err_flag(n)
        real(dp), optional, intent(in) :: lambda_SC(n), remainder(n)

        integer :: u, i
        logical :: with_sc

        with_sc = present(lambda_SC)

        open (newunit=u, file=filename, status="replace", action="write")
        write (u, "(A)") "# asymptotic bootstrap coefficient lambda_bB"
        write (u, "(A)") &
            "# lambda^{off}_bB = Lambda_A/sqrt(nu_star) + Lambda_B/nu_star"
        write (u, "(A,A)") "# git_hash: ", trim(git_hash)
        write (u, "(A,ES23.15)") "# R [m]: ", R
        if (with_sc) then
            write (u, "(A)") &
                "# s_tor  Lambda_A  Lambda_B  nu_star_crit" &
                //"  Lambda_S  err_flag  lambda_SC_bB  remainder"
        else
            write (u, "(A)") &
                "# s_tor  Lambda_A  Lambda_B  nu_star_crit" &
                //"  Lambda_S  err_flag"
        end if
        do i = 1, n
            if (with_sc) then
                write (u, "(ES23.15,4(1X,ES23.15),1X,I2,2(1X,ES23.15))") &
                    s_tor(i), Lambda_A(i), Lambda_B(i), &
                    nu_star_crit(i), Lambda_S(i), err_flag(i), &
                    lambda_SC(i), remainder(i)
            else
                write (u, "(ES23.15,4(1X,ES23.15),1X,I2)") &
                    s_tor(i), Lambda_A(i), Lambda_B(i), &
                    nu_star_crit(i), Lambda_S(i), err_flag(i)
            end if
        end do
        close (u)
    end subroutine write_dat_output

end program rabe
