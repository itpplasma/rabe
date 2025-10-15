module test_calc_eta_integrand_mod
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none
    private

    type :: F_t
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: value
    end type F_t

    class(field_t), pointer, private :: field_to_test => null()
    type(fieldline_t), private :: this_fieldline
    real(dp), private :: this_eta

    public :: average_eta_integrands_sum_equal_average_F
    public :: plot_eta_integrands

contains

    function average_eta_integrands_sum_equal_average_F(field, &
                                                        fieldlines, &
                                                        M_pol, &
                                                        N_tor, &
                                                        eta_integrands) result(is_equal)
        use utils, only: not_same
        use shaing_callen_mod, only: eta_integrand_t
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: M_pol, N_tor
        type(eta_integrand_t), dimension(:) :: eta_integrands
        logical :: is_equal

        real(dp), parameter :: reltol = 1e-4, abstol = 1e-10
        integer, parameter :: n_eta = 10
        integer :: this
        integer :: n_fieldlines
        real(dp), dimension(:), allocatable :: eta_grid
        type(F_t), dimension(:), allocatable :: Fs
        real(dp), dimension(:), allocatable :: average_F, average_integrands_sum

        is_equal = .true.

        n_fieldlines = size(fieldlines)
        allocate (Fs(n_fieldlines))
        allocate (eta_grid, source=eta_integrands(1)%eta_grid)
        allocate (field_to_test, source=field)
        do this = 1, n_fieldlines
            call calc_F(fieldlines(this), eta_grid, M_pol, N_tor, Fs(this))
        end do
        deallocate (field_to_test)
        deallocate (eta_grid)

        average_F = calc_average_F(Fs)
        average_integrands_sum = calc_average_eta_integrands_sum(eta_integrands, &
                                                                 M_pol, &
                                                                 N_tor, &
                                                                 fieldlines(1)%iota)

        if (not_same(average_F, average_integrands_sum, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "average_eta_integrands_sum not equal average_F:"
            print *, "average_integrands_sum: ", average_integrands_sum
            print *, "average_F: ", average_F
            print *, "maximum error: ", maxval(abs(average_F - &
                                                   average_integrands_sum) &
                                               /max(abs(average_F), 1.0_dp))
            is_equal = .false.
        end if

        deallocate (Fs, average_F, average_integrands_sum)

    end function average_eta_integrands_sum_equal_average_F

    function calc_average_eta_integrands_sum(eta_integrands, &
                                             M_pol, &
                                             N_tor, &
                                             iota) result(average_sum)
        use shaing_callen_mod, only: eta_integrand_t
        type(eta_integrand_t), dimension(:), intent(in) :: eta_integrands
        real(dp), intent(in) :: M_pol, N_tor, iota
        real(dp), dimension(:), allocatable :: average_sum

        integer :: this, n_fieldlines, n_eta

        n_fieldlines = size(eta_integrands)
        n_eta = size(eta_integrands(n_fieldlines)%eta_grid)
        allocate (average_sum(n_eta))
        average_sum = 0.0_dp
        do this = 1, n_fieldlines
            average_sum = average_sum &
                          + eta_integrands(this)%F1*M_pol/(M_pol*iota - N_tor) &
                          + eta_integrands(this)%F2 &
                          + eta_integrands(this)%F3
        end do
    end function calc_average_eta_integrands_sum

    function calc_average_F(Fs) result(average_F)
        type(F_t), dimension(:), intent(in) :: Fs
        real(dp), dimension(:), allocatable :: average_F

        integer :: this, n_fieldlines, n_eta

        n_fieldlines = size(Fs)
        n_eta = size(Fs(n_fieldlines)%eta_grid)
        allocate (average_F(n_eta))
        average_F = 0.0_dp
        do this = 1, n_fieldlines
            average_F = average_F + Fs(this)%value
        end do
    end function calc_average_F

    subroutine calc_F(fieldline, eta_grid, M_pol, N_tor, F)
        use fieldline_integrands, only: calc_lambda_squared
        use shaing_callen_mod, only: get_eta_integration_grid
        use shaing_callen_mod, only: get_phi_integration_grid
        use shaing_callen_mod, only: cumint
        use integrate, only: sum_trapez_1d

        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), intent(in) :: M_pol, N_tor
        type(F_t), intent(out) :: F

        real(dp), dimension(:), allocatable :: phi_grid
        real(dp) :: phi

        integer :: this, that, n_phi

        real(dp), dimension(:), allocatable :: antider_dBdtheta_over_lambda_cubed
        real(dp), dimension(:), allocatable :: phi_integrand_F

        this_fieldline = fieldline

        phi_grid = get_phi_integration_grid(fieldline)

        call allocate_F(eta_grid, F)
        n_phi = size(phi_grid)
        allocate (antider_dBdtheta_over_lambda_cubed(n_phi))
        allocate (phi_integrand_F(n_phi))

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            antider_dBdtheta_over_lambda_cubed = cumint(phi_grid, &
                                                    integral_dBdtheta_over_lambda_cubed)

            phi_integrand_F = wrapper_one_over_lambda(fieldline%phi_max(1)) &
                              *M_pol/(N_tor - M_pol*fieldline%iota) &
                              - 0.5_dp*this_eta*antider_dBdtheta_over_lambda_cubed
            do that = 1, n_phi
                phi_integrand_F(that) = wrapper_lambda_over_B_squared(phi_grid(that)) &
                                        *phi_integrand_F(that)
            end do
            F%value(this) = sum_trapez_1d(phi_grid, phi_integrand_F)
        end do
        F%value = -eta_grid*F%value

        deallocate (phi_grid)
        deallocate (antider_dBdtheta_over_lambda_cubed, phi_integrand_F)

    end subroutine calc_F

    subroutine allocate_F(eta_grid, F)
        real(dp), dimension(:), intent(in) :: eta_grid
        type(F_t), intent(inout) :: F

        integer :: n_eta

        if (allocated(F%eta_grid)) deallocate (F%eta_grid)
        if (allocated(F%value)) deallocate (F%value)

        n_eta = size(eta_grid)
        allocate (F%eta_grid, source=eta_grid)
        allocate (F%value(n_eta))
    end subroutine allocate_F

    function wrapper_one_over_lambda(phi)
        use fieldline_integrands, only: calc_lambda_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_one_over_lambda

        real(dp) :: theta, B, lambda

        theta = this_fieldline%get_theta(phi)
        call field_to_test%compute_B_mod(theta, phi, B)
        lambda = sqrt(calc_lambda_squared(B, this_eta))
        wrapper_one_over_lambda = 1.0_dp/lambda
    end function wrapper_one_over_lambda

    function integral_dBdtheta_over_lambda_cubed(phi_start, phi_end)
        use integrate, only: integrate_1d_substituted

        real(dp), intent(in) :: phi_start, phi_end
        real(dp) :: integral_dBdtheta_over_lambda_cubed

        call integrate_1d_substituted(wrapper_dBdtheta_over_lambda_cubed, &
                                      phi_start, &
                                      phi_end, &
                                      integral_dBdtheta_over_lambda_cubed)

    end function integral_dBdtheta_over_lambda_cubed

    function wrapper_dBdtheta_over_lambda_cubed(phi)
        use fieldline_integrands, only: calc_lambda_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_dBdtheta_over_lambda_cubed

        real(dp) :: theta, B, dB_dx(3), dB_dtheta, lambda

        theta = this_fieldline%get_theta(phi)
        call field_to_test%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        lambda = sqrt(calc_lambda_squared(B, this_eta))
        dB_dtheta = dB_dx(2)

        wrapper_dBdtheta_over_lambda_cubed = dB_dtheta/(lambda**3.0_dp)
    end function wrapper_dBdtheta_over_lambda_cubed

    function wrapper_lambda_over_B_squared(phi)
        use fieldline_integrands, only: lambda_over_B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_lambda_over_B_squared

        real(dp) :: theta

        theta = this_fieldline%get_theta(phi)
        wrapper_lambda_over_B_squared = lambda_over_B_squared(field_to_test, &
                                                              theta, &
                                                              phi, &
                                                              this_eta)
    end function wrapper_lambda_over_B_squared

    subroutine plot_eta_integrands(eta_integrands, M_pol, N_tor, iota)
        use myplot_module, only: myplot
        use shaing_callen_mod, only: eta_integrand_t

        type(eta_integrand_t), dimension(:), intent(in) :: eta_integrands
        real(dp), intent(in) :: M_pol, N_tor, iota

        type(myplot) :: plt
        real(dp), dimension(:), allocatable :: average_F1, average_F2, average_F3
        real(dp), dimension(:), allocatable :: eta_grid
        integer :: this, n_fieldlines, n_eta

        n_fieldlines = size(eta_integrands)
        n_eta = size(eta_integrands(1)%eta_grid)

        allocate (average_F1(n_eta))
        allocate (average_F2(n_eta))
        allocate (average_F3(n_eta))
        allocate (eta_grid(n_eta))

        eta_grid = eta_integrands(1)%eta_grid

        average_F1 = 0.0_dp
        average_F2 = 0.0_dp
        average_F3 = 0.0_dp

        do this = 1, n_fieldlines
            average_F1 = average_F1 + eta_integrands(this)%F1*M_pol/(M_pol*iota - N_tor)
            average_F2 = average_F2 + eta_integrands(this)%F2
            average_F3 = average_F3 + eta_integrands(this)%F3
        end do

        call plt%initialize(xlabel="$\eta$ [T$^{-1}$]", &
                            ylabel="$F(\eta)$ contributions", &
                            legend=.true., &
                            figsize=[10, 8])

        call plt%add_plot(eta_grid, &
                          average_F1, &
                          label="$F_1 \cdot M/(M\iota-N)$", &
                          linestyle="-")

        call plt%add_plot(eta_grid, &
                          average_F2, &
                          label="$F_2$", &
                          linestyle="--")

        call plt%add_plot(eta_grid, &
                          average_F3, &
                          label="$F_3$", &
                          linestyle=":")

        call plt%show()

        deallocate (average_F1, average_F2, average_F3, eta_grid)

    end subroutine plot_eta_integrands

end module test_calc_eta_integrand_mod
