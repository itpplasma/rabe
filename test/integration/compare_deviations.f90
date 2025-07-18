program compare_deviations
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines
    use deviation, only: calc_deviation

    use plot_quantities, only: external_data_t
    use readers, only: read_column

    use fieldline_integrals, only: modes_t
    use distribution_function, only: get_g_modes_from_fieldlines
    use distribution_function, only: get_modes
    use distribution_function, only: get_offset_from_distribution

    implicit none

    character(len=*), parameter :: bc_filename = "input/nautilus.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 4.0_dp
    real(dp), parameter :: sign_sqrtg = -1.0_dp ! theta goes counter-clockwise
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.25_dp

    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 8e-7
    integer, parameter :: n_fieldlines = 151

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    type(external_data_t) :: g_neo2
    real(dp), parameter :: nu_star = 6e-5
    real(dp), parameter :: scaling = 100.0_dp/5.884 !100/bmod0 -> g_(NEO-2) in [cm]
    integer :: n_columns
    integer :: n_neo2

    !------------------Taken from NEO-2 output---------------------------------!
    real(dp), parameter :: ds_dr = 0.00752119_dp*100.0_dp ! [1/m] called "avnabpsi"
    !--------------------------------------------------------------------------!

    real(dp) :: dr_dAphi
    real(dp), parameter :: offset_neo2 = -0.00447589_dp

    real(dp) :: offset_g_rabe, offset_g_neo2
    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor, R
    real(dp) :: off_factor_A, off_factor_B
    real(dp) :: offset_rabe

    type(modes_t) :: g_rabe_modes, g_neo2_modes

    real(dp), parameter :: reltol = 1e-14
    logical :: test_failed

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_of_fieldlines(fieldlines, &
                                  theta_0, &
                                  iota, &
                                  field, &
                                  M_pol, &
                                  N_tor, &
                                  phi_tol)

    g_neo2%label = "NEO-2: $\hat{g}_0 \frac{\mathrm{bmod0}}{100}$"
    n_columns = 1
    call read_column("input/theta.dat", g_neo2%x, n_columns, 1)
    n_columns = size(g_neo2%x, dim=1)
    call read_column("input/gvpar0.dat", &
                     g_neo2%y, &
                     n_columns, &
                     n_columns/8)
    g_neo2%y = g_neo2%y/scaling

    call get_g_modes_from_fieldlines(fieldlines, field, nu_star, g_rabe_modes)

    n_neo2 = size(g_neo2%x)
    call get_modes(g_neo2%x(1:n_neo2 - 1), g_neo2%y(1:n_neo2 - 1), g_neo2_modes)

    offset_g_rabe = get_offset_from_distribution(g_rabe_modes, fieldlines)
    offset_g_neo2 = get_offset_from_distribution(g_neo2_modes, fieldlines)

    covariant_factor = (field%B_phi_covariant + field%B_theta_covariant*iota)
    call calc_deviation(fieldlines, deviation_A, deviation_B)
    R = field%R
    off_factor_A = deviation_A*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5_dp*R*pi
    offset_rabe = off_factor_A/sqrt(nu_star) + off_factor_B/nu_star

    dr_dAphi = 1.0_dp/(ds_dr*field%psi_tor_edge)*sign_sqrtg

    offset_g_neo2 = offset_g_neo2*dr_dAphi
    offset_g_rabe = offset_g_rabe*dr_dAphi
    offset_rabe = offset_rabe*dr_dAphi

    test_failed = .false.

    if (not_same(offset_rabe, offset_g_rabe, &
                 reltol_in=reltol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "compare_deviations failed: compare to offset due to distribution"
        test_failed = .true.
    end if

    print *, "sum of serperate contributions (Delta A & Delta eta): ", offset_rabe
    print *, "directly from asymptotic distributions function: ", offset_g_rabe
    print *, "from NEO-2 distributions function: ", offset_g_neo2
    print *, "result from NEO-2: ", offset_neo2

    if (test_failed) error stop 1

end program compare_deviations
