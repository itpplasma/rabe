module test_calc_eta_integrand_mod
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    use shaing_callen_mod, only: calc_alternative_F
    implicit none
    private

    public :: average_eta_integrands_sum_equal_F
    public :: plot_eta_integrands

contains

    function average_eta_integrands_sum_equal_F(field, &
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
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: alternative_F, F
        real(dp), dimension(:), allocatable :: average_integrands_sum

        is_equal = .true.

        allocate (eta_grid, source=eta_integrands(1)%eta_grid)

        alternative_F = calc_alternative_F(field, &
                                           fieldlines, &
                                           eta_grid)

        average_integrands_sum = calc_average_eta_integrands_sum(eta_integrands, &
                                                                 M_pol, &
                                                                 N_tor, &
                                                                 fieldlines(1)%iota)

        F = average_integrands_sum/sum(fieldlines%integral_one_over_B_squared)
        if (not_same(alternative_F, F, &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "-------------------------------------------------------------"
            print *, "F not equal alternative_F:"
            print *, "F: ", F
            print *, "alternative_F: ", alternative_F
            print *, "maximum error: ", maxval(abs(alternative_F - F) &
                                               /max(abs(alternative_F), 1.0_dp))
            is_equal = .false.
        end if

        deallocate (alternative_F, average_integrands_sum)

    end function average_eta_integrands_sum_equal_F

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

    subroutine plot_eta_integrands(eta_integrands, fieldlines)
        use myplot_module, only: myplot
        use shaing_callen_mod, only: eta_integrand_t

        type(eta_integrand_t), dimension(:), intent(in) :: eta_integrands
        type(fieldline_t), dimension(:), intent(in) :: fieldlines

        type(myplot) :: plt
        real(dp), dimension(:), allocatable :: average_F1, average_F2, average_F3
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp) :: M_pol, N_tor, iota
        integer :: this, n_fieldlines, n_eta

        n_fieldlines = size(eta_integrands)
        n_eta = size(eta_integrands(1)%eta_grid)
        M_pol = fieldlines(1)%M_pol
        N_tor = fieldlines(1)%N_tor
        iota = fieldlines(1)%iota

        allocate (average_F1(n_eta))
        allocate (average_F2(n_eta))
        allocate (average_F3(n_eta))
        allocate (eta_grid(n_eta))

        eta_grid = eta_integrands(1)%eta_grid

        average_F1 = 0.0_dp
        average_F2 = 0.0_dp
        average_F3 = 0.0_dp

        do this = 1, n_fieldlines
            average_F1 = average_F1 + eta_integrands(this)%F1
            average_F2 = average_F2 + eta_integrands(this)%F2
            average_F3 = average_F3 + eta_integrands(this)%F3
        end do
        average_F1 = average_F1/sum(fieldlines%integral_one_over_B_squared)
        average_F2 = average_F2/sum(fieldlines%integral_one_over_B_squared)
        average_F3 = average_F3/sum(fieldlines%integral_one_over_B_squared)

        call plt%initialize(xlabel="$\eta$ [T$^{-1}$]", &
                            ylabel="$F(\eta)$ contributions", &
                            legend=.true., &
                            figsize=[15, 12])

        call plt%add_plot(eta_grid, &
                          average_F1*M_pol/(M_pol*iota - N_tor), &
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

        call plt%add_plot(eta_grid, &
                          average_F1*M_pol/(M_pol*iota - N_tor) + &
                          average_F2 + &
                          average_F3, &
                          label="sum", &
                          linestyle="o")
        call plt%add_plot(eta_grid, &
                          eta_grid*M_pol/(M_pol*iota - N_tor), &
                          label="analytic sum", &
                          linestyle="-")

        call plt%show()

        deallocate (average_F1, average_F2, average_F3, eta_grid)

    end subroutine plot_eta_integrands

end module test_calc_eta_integrand_mod
