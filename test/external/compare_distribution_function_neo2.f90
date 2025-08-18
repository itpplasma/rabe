program compare_distribution_function_neo2
    use constants, only: dp, pi
    use utils, only: linspace
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use make_fieldline, only: make_flock_of_fieldlines

    use plot_quantities, only: plot_distribution_function
    use plot_quantities, only: external_data_t
    use readers, only: read_column

    use fieldline_integrals, only: modes_t
    use distribution_function, only: get_g_modes_from_fieldlines
    use distribution_function, only: get_modes
    use plot_quantities, only: compare_modes

    implicit none

    character(len=*), parameter :: bc_filename = "input/nautilus.bc"

    !------------------Taken from .bc file-------------------------------------!
    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 4.0_dp
    !--------------------------------------------------------------------------!

    real(dp), parameter :: stor = 0.25_dp

    type(neo_field_t) :: field

    real(dp), parameter :: phi_tol = 8e-7
    integer, parameter :: n_fieldlines = 151

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota, nfp
    type(fieldline_t), dimension(n_fieldlines) :: fieldlines

    type(external_data_t) :: g_neo2
    real(dp), parameter :: nu_star = 6e-5
    real(dp), parameter :: scaling = 100.0_dp/5.884 !100/bmod0 -> g_(NEO-2) in [cm]
    integer :: n_columns
    integer :: n_neo2

    type(modes_t), dimension(2) :: modes
    character(len=1024), dimension(2) :: labels

    call field%neo_field_init(bc_filename, stor)
    iota = field%iota
    nfp = field%nfp
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

    g_neo2%label = "NEO-2: $\hat{g}_0 \frac{\mathrm{bmod0}}{100}$"
    n_columns = 1
    call read_column("input/theta.dat", g_neo2%x, n_columns, 1)
    n_columns = size(g_neo2%x, dim=1)
    call read_column("input/gvpar0.dat", &
                     g_neo2%y, &
                     n_columns, &
                     n_columns/8)
    g_neo2%y = g_neo2%y/scaling
    call plot_distribution_function(fieldlines, field, nu_star, g_neo2)

    call get_g_modes_from_fieldlines(fieldlines, field, nu_star, modes(1))
    labels(1) = "rabe: $g_\mathrm{off}$"

    n_neo2 = size(g_neo2%x)
    call get_modes(g_neo2%x(1:n_neo2 - 1), g_neo2%y(1:n_neo2 - 1), modes(2))
    labels(2) = g_neo2%label

    call compare_modes(modes, labels)

end program compare_distribution_function_neo2
