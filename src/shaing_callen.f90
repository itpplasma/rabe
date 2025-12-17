module shaing_callen_mod
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_integration, only: integrate_over_eta_grid
    implicit none

contains

    function calc_trapped_fraction(field, &
                                   fieldlines, &
                                   n_eta) result(trapped_fraction)
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        integer, intent(in) :: n_eta
        real(dp) :: trapped_fraction

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp), dimension(:), allocatable :: integrand

        eta_grid = get_eta_integration_grid(fieldlines(1)%eta_b, n_eta)
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        allocate (integrand(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           fieldlines, &
                                                                           eta_grid)
        integrand = eta_grid*avg_B_squared_over_avg_lambda
        trapped_fraction = 1.0_dp - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)
        deallocate (integrand)
        deallocate (avg_B_squared_over_avg_lambda)

    end function calc_trapped_fraction

    function calc_avg_B_squared_over_avg_lambda(field, fieldlines, eta_grid) &
        result(avg_B_squared_over_avg_lambda)
        use shaing_callen_wrappers, only: this_field
        class(field_t), intent(in) :: field
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda

        integer :: this
        integer :: n_eta, n_fieldlines
        real(dp) :: avg_well_length
        real(dp), dimension(:), allocatable :: avg_lambda_over_B_squared

        n_eta = size(eta_grid)
        n_fieldlines = size(fieldlines)
        allocate (avg_lambda_over_B_squared(n_eta))
        avg_well_length = 0.0_dp
        avg_lambda_over_B_squared = 0.0_dp
        allocate (this_field, source=field)
        do this = 1, n_fieldlines
            avg_lambda_over_B_squared = avg_lambda_over_B_squared + &
                              calc_avg_lambda_over_B_squared(fieldlines(this), eta_grid)
            avg_well_length = avg_well_length + &
                              fieldlines(this)%phi_max(2) - &
                              fieldlines(this)%phi_max(1)
        end do
        deallocate (this_field)
        allocate (avg_B_squared_over_avg_lambda(n_eta))

        if (any(avg_lambda_over_B_squared <= 0.0_dp)) then
            print *, "Error in calc_avg_B_squared_over_avg_lambda: ", &
                "average lambda/B^2 is not positiv"
            print *, "avg_lambda_over_B_squared: ", avg_lambda_over_B_squared
            error stop
        end if

        avg_B_squared_over_avg_lambda = avg_well_length/ &
                                        avg_lambda_over_B_squared

        deallocate (avg_lambda_over_B_squared)

    end function calc_avg_B_squared_over_avg_lambda

    function calc_avg_lambda_over_B_squared(fieldline, eta_grid) &
        result(avg_lambda_over_B_squared)
        use integrate, only: integrate_1d_substituted
        use shaing_callen_wrappers, only: wrapper_lambda_over_B_squared
        use shaing_callen_wrappers, only: this_eta, null_eta
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: avg_lambda_over_B_squared

        integer :: this

        this_fieldline = fieldline

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            call integrate_1d_substituted(wrapper_lambda_over_B_squared, &
                                          fieldline%phi_max(1), &
                                          fieldline%phi_max(2), &
                                          avg_lambda_over_B_squared(this))
            this_eta = null_eta
        end do

        this_fieldline = null_fieldline

    end function calc_avg_lambda_over_B_squared

end module shaing_callen_mod
