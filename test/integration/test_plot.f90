program test_plot
    use plot_mod, only: plot
    use constants, only: dp, pi
    use utils, only: linspace
    implicit none

    integer, parameter :: n_points = 100
    real(dp), dimension(n_points) :: x, f

    call linspace(0.0_dp, 2.0_dp*pi, n_points, x)
    f = sin(x)

    call plot(x, f)

end program test_plot
