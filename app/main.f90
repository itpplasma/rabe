program rabe
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use shaing_callen_mod, only: calc_trapped_fraction
    use shaing_callen_mod, only: get_non_omnigenous_remainder
    use netcdf_mod, only: netcdf_t
    use git_version, only: git_hash

    use read_file, only: read_namelist
    use read_file, only: bc_filename, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         sign_sqrtg, &
                         phi_tol, &
                         n_fieldlines, &
                         should_calc_shaing_callen, &
                         n_eta

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "rabe.nc"

    type(neo_field_t) :: field

    integer :: n_stor
    integer :: this

    real(dp) :: R ![m]
    type(surface_average_t) :: average
    real(dp) :: dr_dAtheta ![rad/Tm/]
    real(dp) :: iota
    real(dp) :: nfp

    real(dp), dimension(:), allocatable :: theta_0
    real(dp), dimension(:), allocatable :: temp

    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp), dimension(:), allocatable :: C_A, C_B
    real(dp), dimension(:), allocatable :: nu_star_limit

    real(dp) :: trapped_fraction
    real(dp), dimension(:), allocatable :: lambda_SC, remainder

    type(netcdf_t) :: nc_output
    character(len=*), parameter :: dim_name = "surface"

    call read_namelist(input_file)

    n_stor = size(s_tor)
    allocate (C_A(n_stor))
    allocate (C_B(n_stor))
    allocate (nu_star_limit(n_stor))
    if (should_calc_shaing_callen) then
        allocate (lambda_SC(n_stor))
        allocate (remainder(n_stor))
    end if

    allocate (fieldlines(n_fieldlines))
    allocate (theta_0(n_fieldlines))
    allocate (temp(n_fieldlines + 1))
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    do this = 1, n_stor
        if (this == 1) then
            call field%neo_field_init(bc_filename, s_tor(1))
        else
            call field%neo_change_stor(s_tor(this))
        end if
        iota = field%iota
        nfp = field%nfp

        call make_flock_of_fieldlines(fieldlines, &
                                      theta_0, &
                                      iota, &
                                      field, &
                                      M_pol, &
                                      N_tor, &
                                      nfp, &
                                      phi_tol)

        call calc_deviation(fieldlines, deviation_A, deviation_B)

        covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
        call calc_surface_averages(fieldlines, average)
        dr_dAtheta = sign_sqrtg*sign(1.0_dp, field%psi_tor_edge)/average%sqrt_g11
        R = field%R
        C_A(this) = deviation_A*dr_dAtheta*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
        C_B(this) = deviation_B*0.5*R*pi*dr_dAtheta
        nu_star_limit(this) = R/fieldlines(1)%I_ref*(fieldlines(1)%eta_b - &
                                          1.0_dp/minval(fieldlines%B_max(1)))**2.0_dp/ &
                              fieldlines(1)%eta_b*0.25_dp*pi/3.0_dp
        nu_star_limit(this) = nu_star_limit(this)/covariant_factor

        print *, "s_tor: ", s_tor(this)
        if (should_calc_shaing_callen) then
            trapped_fraction = calc_trapped_fraction(field, fieldlines, n_eta)
            lambda_SC(this) = (field%B_phi_covariant*M_pol + &
                               field%B_theta_covariant*N_tor)/ &
                              (M_pol*iota - N_tor)*trapped_fraction
            lambda_SC(this) = lambda_SC(this)*dr_dAtheta
            remainder(this) = get_non_omnigenous_remainder(field, fieldlines, n_eta)
            remainder(this) = remainder(this)*covariant_factor*dr_dAtheta* &
                              nfp/(M_pol*iota - N_tor)
            print *, "omnigenous lambda_SC_bB: ", lambda_SC(this)
            print *, "non-omnigneous remainder: ", remainder(this)
        end if
        print *, "1/sqrt(nu_star) factor: ", C_A(this)
        print *, "1/nu_star factor: ", C_B(this)

    end do

    call nc_output%create(output_file)
    call nc_output%add_global_attribute("title", &
                                        "asymptotic bootstrap coefficient lambda_bB")
    call nc_output%add_global_attribute("definition", &
                                    "lambda^{off}_bB = C_A/sqrt(nu_star) + C_B/nu_star")
    call nc_output%add_global_attribute("git_hash", git_hash)
    call nc_output%def_dim(dim_name, n_stor)
    call nc_output%add_real_1d("s_tor", dim_name)
    call nc_output%add_real_attr("s_tor", "long_name", &
                                 "normalized toroidal flux label")
    call nc_output%add_real_attr("s_tor", "unit", &
                                 "[1]")
    call nc_output%add_real_1d("C_A", dim_name)
    call nc_output%add_real_attr("C_A", "long_name", &
                                 "1/sqrt(nu_star) factor")
    call nc_output%add_real_attr("C_A", "unit", &
                                 "[1]")
    call nc_output%add_real_1d("C_B", dim_name)
    call nc_output%add_real_attr("C_B", "long_name", &
                                 "1/nu_star factor")
    call nc_output%add_real_attr("C_B", "unit", &
                                 "[1]")
    call nc_output%add_real_1d("nu_star_limit", dim_name)
    call nc_output%add_real_attr("nu_star_limit", "long_name", &
                                "collisionality limit for validity of asymptotic model")
    call nc_output%add_real_attr("nu_star_limit", "unit", &
                                 "[1]")
    if (should_calc_shaing_callen) then
        call nc_output%add_real_1d("lambda_SC_bB", dim_name)
        call nc_output%add_real_attr("lambda_SC_bB", "long_name", &
                                     "omnigenous Shaing-Callen coefficient")
        call nc_output%add_real_attr("lambda_SC_bB", "unit", "[1]")
        call nc_output%add_real_1d("remainder", dim_name)
        call nc_output%add_real_attr("remainder", "long_name", &
                                "non-omnigenous remainder of Shaing-Callen coefficient")
        call nc_output%add_real_attr("remainder", "unit", "[1]")

        call nc_output%write_real_1d("lambda_SC_bB", lambda_SC)
        call nc_output%write_real_1d("remainder", remainder)
    end if
    call nc_output%write_real_1d("C_A", C_A)
    call nc_output%write_real_1d("C_B", C_B)
    call nc_output%write_real_1d("nu_star_limit", nu_star_limit)
    call nc_output%write_real_1d("s_tor", s_tor)
    call nc_output%close()

    if (allocated(C_A)) deallocate (C_A)
    if (allocated(C_B)) deallocate (C_B)
    if (allocated(nu_star_limit)) deallocate (nu_star_limit)
    if (allocated(lambda_SC)) deallocate (lambda_SC)
    if (allocated(remainder)) deallocate (remainder)

    if (allocated(fieldlines)) deallocate (fieldlines)
    if (allocated(theta_0)) deallocate (theta_0)
    if (allocated(temp)) deallocate (temp)

end program rabe
