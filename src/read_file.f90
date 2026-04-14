module read_file
    use constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

    character(len=100), public, protected :: field_file
    real(dp), public, protected :: M_pol
    real(dp), public, protected :: N_tor
    real(dp), dimension(:), allocatable, public, protected :: s_tor
    real(dp), public, protected :: ds_dr ![1/cm]
    real(dp), public, protected :: sign_sqrtg
    real(dp), public, protected :: phi_tol
    integer, public, protected :: max_n_fieldlines
    logical, public, protected :: should_calc_shaing_callen
    integer, public, protected :: n_eta
    real(dp), dimension(:), allocatable, public, protected :: Omega_hat

    namelist /rabe_config/ &
        field_file, &
        M_pol, &
        N_tor, &
        s_tor, &
        sign_sqrtg, &
        phi_tol, &
        max_n_fieldlines, &
        should_calc_shaing_callen, &
        n_eta, &
        Omega_hat

contains

    subroutine read_namelist(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: ios, unit
        logical :: file_exists

        real(dp), dimension(:), allocatable :: s_tor_temp
        integer, parameter :: n_stor_max = 1000
        integer :: n_stor
        real(dp), dimension(:), allocatable :: Omega_hat_temp
        integer, parameter :: n_Omega_hat_max = 1000
        integer :: n_Omega_hat
        real(dp) :: nan_value

        inquire (file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "Error in read_namelist: file not found: ", filename
            error stop
        end if

        open (newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "Error in read_namelist: cannot open file: ", filename
            error stop
        end if

        ! Default values
        should_calc_shaing_callen = .false.
        n_eta = 100
        if (allocated(s_tor)) deallocate (s_tor)
        allocate (s_tor(n_stor_max))
        nan_value = ieee_value(nan_value, ieee_quiet_nan)
        s_tor = nan_value
        if (allocated(Omega_hat)) deallocate (Omega_hat)
        allocate (Omega_hat(n_Omega_hat_max))
        Omega_hat = nan_value

        read (unit, nml=rabe_config, iostat=ios)
        if (ios /= 0) then
            print *, "Error in read_namelist:"
            print *, "iostat = ", ios
            error stop
        end if

        n_stor = count(.not. ieee_is_nan(s_tor))
        if (n_stor == 0) then
            print *, "Error in read_namelist: no s_tor values provided in namelist"
            error stop
        end if
        if (n_stor >= n_stor_max) then
           print *, "Error in read_namelist: too many s_tor values provided in namelist"
            error stop
        end if
        ! Need temporary storage to reduce the size of s_tor to amount of inputs
        if (allocated(s_tor_temp)) deallocate (s_tor_temp)
        allocate (s_tor_temp(n_stor))
        s_tor_temp = pack(s_tor,.not. ieee_is_nan(s_tor))
        deallocate (s_tor)
        allocate (s_tor(n_stor))
        s_tor = s_tor_temp
        deallocate (s_tor_temp)

        n_Omega_hat = count(.not. ieee_is_nan(Omega_hat))
        if (n_Omega_hat >= n_stor_max) then
            print *, "Error in read_namelist: too many Omega_hat values ", &
                "provided in namelist"
            error stop
        end if
        if (allocated(Omega_hat_temp)) deallocate (Omega_hat_temp)
        allocate (Omega_hat_temp(n_Omega_hat))
        Omega_hat_temp = pack(Omega_hat,.not. ieee_is_nan(Omega_hat))
        deallocate (Omega_hat)
        allocate (Omega_hat(n_Omega_hat))
        Omega_hat = Omega_hat_temp
        deallocate (Omega_hat_temp)

        call check_if_valid_namelist()

        close (unit)
    end subroutine read_namelist

    subroutine check_if_valid_namelist()
        use make_fieldline, only: is_not_integer
        use utils, only: not_same

        real(dp), parameter :: tol = 1e-15
        logical :: is_valid

        is_valid = .true.

        if (len(trim(field_file)) == 0) then
            print *, "field_file is empty!"
        end if
        if (ieee_is_nan(M_pol)) then
            print *, "M_pol is NaN!"
            is_valid = .false.
        end if
        if (ieee_is_nan(N_tor)) then
            print *, "N_tor is NaN!"
            is_valid = .false.
        end if
        if (any(ieee_is_nan(s_tor))) then
            print *, "s_tor is NaN!"
            is_valid = .false.
        end if
        if (ieee_is_nan(sign_sqrtg)) then
            print *, "sign_sqrtg is NaN!"
            is_valid = .false.
        end if
        if (ieee_is_nan(phi_tol)) then
            print *, "phi_tol is NaN!"
            is_valid = .false.
        end if

        if (is_not_integer(M_pol, tol)) then
            print *, "M_pol must be integer"
            is_valid = .false.
        end if
        if (is_not_integer(N_tor, tol)) then
            print *, "N_tor must be integer"
            is_valid = .false.
        end if
        if (not_same(sign_sqrtg, 1.0_dp, reltol_in=0.0_dp, abstol_in=tol) .and. &
            not_same(sign_sqrtg, -1.0_dp, reltol_in=0.0_dp, abstol_in=tol)) then
            print *, "sign_sqrtg must be +-1"
            is_valid = .false.
        end if
        if (any(s_tor <= 0.0_dp) .or. any(s_tor >= 1.0_dp)) then
            print *, "s_tor must be in ]0.0, 1.0["
            is_valid = .false.
        end if
        if (phi_tol <= 0.0_dp) then
            print *, "phi_tol must be positive"
            is_valid = .false.
        end if
        if (max_n_fieldlines <= 1) then
            print *, "max_n_fieldlines must be bigger than 1"
            is_valid = .false.
        end if
        if (should_calc_shaing_callen .and. n_eta <= 3) then
            print *, "n_eta must be bigger than 3"
        end if

        if (.not. is_valid) error stop

    end subroutine check_if_valid_namelist

end module read_file
