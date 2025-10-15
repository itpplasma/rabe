module test_shaing_callen_mod
    use constants, only: dp, pi

    implicit none

contains

    function calc_quasi_symmetric_trapped_fraction(field, &
                                                   eta_b) result(trapped_fraction)
        use field_base, only: field_t
        use shaing_callen_mod, only: get_eta_integration_grid
        use shaing_callen_mod, only: integrate_over_eta_grid
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: eta_b
        real(dp) :: trapped_fraction

        integer, parameter :: n_eta = 1000
        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(n_eta) :: integrand
        real(dp) :: eta, average

        integer :: this

        eta_grid = get_eta_integration_grid(eta_b, n_eta)
        do this = 1, n_eta
            eta = eta_grid(this)
            average = get_theta_average_lambda_over_B_squared(field, eta)
            integrand(this) = eta/average
        end do
        trapped_fraction = 1.0_dp - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)
    end function calc_quasi_symmetric_trapped_fraction

    function get_theta_average_lambda_over_B_squared(field, eta) result(average)
        use field_base, only: field_t
        use fieldline_integrands, only: lambda_over_B_squared
        use utils, only: linspace

        class(field_t), intent(in) :: field
        real(dp), intent(in) :: eta
        real(dp) :: average

        integer, parameter :: n_theta = 100
        real(dp), dimension(n_theta) :: thetas
        real(dp), dimension(n_theta - 1) :: integrand
        real(dp), parameter :: phi = 0.0_dp

        integer :: this

        call linspace(0.0_dp, 2.0_dp*pi, n_theta, thetas)
        do this = 1, n_theta - 1
            integrand(this) = lambda_over_B_squared(field, &
                                                    thetas(this), &
                                                    phi, &
                                                    eta)
        end do
        average = sum(integrand)/real(n_theta - 1, kind=dp)
    end function get_theta_average_lambda_over_B_squared

end module test_shaing_callen_mod
