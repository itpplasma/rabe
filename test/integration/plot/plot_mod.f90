module plot_mod
    use pyplot_module, only: pyplot
    use constants, only: dp

    implicit none

contains

    subroutine plot(x, f)

        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x)), intent(in) :: f

        type(pyplot) :: plt

        call plt%initialize(grid=.true., legend=.true.)
        call plt%add_plot(x, f, label="f", linestyle='b-o', markersize=5, linewidth=2)
        call plt%showfig()
    end subroutine plot

end module plot_mod
