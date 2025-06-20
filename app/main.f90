program rabe
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"

    character(len=100) :: bc_filename
    real(dp) :: M_pol
    real(dp) :: N_tor
    real(dp) :: s_tor
    real(dp) :: ds_dr ![1/cm]
    real(dp) :: J_pol_over_N_tor
    real(dp) :: I_tor
    real(dp) :: flux_edge ![Tm^2]
    real(dp) :: sign_sqrtg
    real(dp) :: R
    real(dp) :: phi_tol
    integer :: n_fieldlines

    type(neo_field_t) :: field

    real(dp) :: psi_edge ![Tm^2/rad]
    real(dp) :: dr_dpsi ![rad/Tm/]
    real(dp) :: iota

    real(dp), dimension(:), allocatable :: theta_0
    real(dp), dimension(:), allocatable :: temp

    type(fieldline_t), dimension(:), allocatable :: fieldlines

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    namelist /rabe_config/ &
        bc_filename, &
        M_pol, &
        N_tor, &
        s_tor, &
        ds_dr, &
        J_pol_over_N_tor, &
        I_tor, &
        flux_edge, &
        sign_sqrtg, &
        R, &
        phi_tol, &
        n_fieldlines

    call read_namelist(input_file)
    psi_edge = flux_edge/(2.0_dp*pi)
    dr_dpsi = 1.0_dp/(ds_dr*100.0_dp)/psi_edge

    call field%neo_field_init(bc_filename, s_tor)
    iota = field%iota

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
                                  phi_tol)

    call calc_deviation(fieldlines, field, deviation_A, deviation_B)

    covariant_factor = -2.0_dp*1e-7*(J_pol_over_N_tor*abs(N_tor) + I_tor*iota)
    off_factor_A = deviation_A*dr_dpsi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dpsi

    print *, "1/sqrt(nu_star) factor: ", off_factor_A
    print *, "1/nu_star factor: ", off_factor_B

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
