program rabe
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation
    use shaing_callen_mod, only: calc_trapped_fraction
    use shaing_callen_mod, only: get_non_omnigenous_remainder
    use shaing_callen_integration, only: get_eta_integration_grid
    use netcdf_mod, only: netcdf_t
    use git_version, only: git_hash

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "rabe.nc"

    character(len=100) :: bc_filename
    real(dp) :: M_pol
    real(dp) :: N_tor
    real(dp) :: s_tor
    real(dp) :: ds_dr ![1/cm]
    real(dp) :: sign_sqrtg
    real(dp) :: phi_tol
    integer :: n_fieldlines
    logical :: should_calc_shaing_callen
    integer :: n_eta

    type(neo_field_t) :: field

    real(dp) :: R ![m]
    real(dp) :: psi_edge ![Tm^2/rad]
    real(dp) :: dr_dpsi ![rad/Tm/]
    real(dp) :: dr_dAtheta ![rad/Tm/]
    real(dp) :: iota
    real(dp) :: nfp

    real(dp), dimension(:), allocatable :: theta_0
    real(dp), dimension(:), allocatable :: temp

    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: C_A, C_B

    real(dp), dimension(:), allocatable :: eta_grid
    real(dp) :: trapped_fraction
    real(dp) :: lambda_SC, remainder

    type(netcdf_t) :: nc_output

    namelist /rabe_config/ &
        bc_filename, &
        M_pol, &
        N_tor, &
        s_tor, &
        ds_dr, &
        sign_sqrtg, &
        phi_tol, &
        n_fieldlines, &
        should_calc_shaing_callen, &
        n_eta

    should_calc_shaing_callen = .false.
    n_eta = 100

    call read_namelist(input_file)

    call field%neo_field_init(bc_filename, s_tor)
    iota = field%iota
    nfp = field%nfp

    allocate (fieldlines(n_fieldlines))
    allocate (theta_0(n_fieldlines))
    allocate (temp(n_fieldlines + 1))

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

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
    dr_dpsi = 1.0_dp/(ds_dr*100.0_dp)/field%psi_tor_edge
    dr_dAtheta = dr_dpsi*sign_sqrtg
    R = field%R
    C_A = deviation_A*dr_dAtheta*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    C_B = deviation_B*0.5*R*pi*dr_dAtheta

    if (should_calc_shaing_callen) then
        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        trapped_fraction = calc_trapped_fraction(field, fieldlines, eta_grid)
        lambda_SC = (field%B_phi_covariant*M_pol + field%B_theta_covariant*N_tor)/ &
                    (M_pol*iota - N_tor)*trapped_fraction
        lambda_SC = lambda_SC*dr_dAtheta
        remainder = get_non_omnigenous_remainder(field, fieldlines, eta_grid)
        remainder = remainder*covariant_factor*dr_dAtheta
        print *, "lambda_SC: ", lambda_SC
        print *, "non-omnigneous remainder: ", remainder
    end if

    print *, "1/sqrt(nu_star) factor: ", C_A
    print *, "1/nu_star factor: ", C_B

    call nc_output%create(output_file)
    call nc_output%add_global_attribute("title", &
                                        "asymptotic bootstrap coefficient lambda_bB")
    call nc_output%add_global_attribute("definition", &
                                    "lambda^{off}_bB = C_A/sqrt(nu_star) + C_B/nu_star")
    call nc_output%add_global_attribute("git_hash", git_hash)
    call nc_output%add_real("C_A")
    call nc_output%add_real_attr("C_A", "long_name", &
                                 "1/sqrt(nu_star) factor")
    call nc_output%add_real_attr("C_A", "unit", &
                                 "[1]")
    call nc_output%add_real("C_B")
    call nc_output%add_real_attr("C_B", "unit", &
                                 "[1]")
    if (should_calc_shaing_callen) then
        call nc_output%add_real("lambda_SC")
        call nc_output%add_real_attr("lambda_SC", "long_name", &
                                     "omnigenous Shaing-Callen coefficient")
        call nc_output%add_real_attr("lambda_SC", "unit", "[1]")
        call nc_output%write_real("lambda_SC", lambda_SC)
    end if
    call nc_output%write_real("C_A", C_A)
    call nc_output%write_real("C_B", C_B)
    call nc_output%close()

    deallocate (fieldlines)
    deallocate (theta_0)
    deallocate (temp)

contains

    subroutine read_namelist(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: ios, unit
        logical :: file_exists

        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "error: file not found: ", filename
            stop
        end if

        open (newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "error: cannot open file: ", filename
            stop
        end if

        read (unit, nml=rabe_config, iostat=ios)
        if (ios /= 0) then
            print *, "error: failed to read namelist rabe_config"
            print *, "iostat = ", ios
            stop
        end if

        close (unit)
    end subroutine read_namelist

end program rabe
