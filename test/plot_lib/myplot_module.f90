module myplot_module
    use pyplot_module, only: pyplot
    use constants, only: dp

    implicit none

    type :: myplot
        type(pyplot) :: plt
    contains
        procedure :: initialize
        procedure :: add_plot
        procedure :: add_contour
        procedure :: show
    end type myplot

contains

    subroutine initialize(self, xlabel, ylabel, legend, figsize, title, fontsize)
        class(myplot), intent(out) :: self

        character(len=*), intent(in), optional :: xlabel, ylabel
        logical, intent(in), optional :: legend
        integer, dimension(2), intent(in), optional :: figsize
        character(len=*), intent(in), optional :: title
        integer, intent(in), optional :: fontsize

        integer :: fontsize_internal
        integer, dimension(2) :: figsize_internal
        logical, parameter :: tight_layout = .true.

        if (present(fontsize)) then
            fontsize_internal = fontsize
        else
            fontsize_internal = 30
        end if

        if (present(figsize)) then
            figsize_internal = figsize
        else
            figsize_internal = [10, 8]
        end if

        call self%plt%initialize(grid=.true., &
                                 xlabel=xlabel, &
                                 ylabel=ylabel, &
                                 legend=legend, &
                                 figsize=figsize_internal, &
                                 title=title, &
                                 raw_strings=.true., &
                                 font_size=fontsize_internal, &
                                 legend_fontsize=fontsize_internal, &
                                 axes_labelsize=fontsize_internal, &
                                 xtick_labelsize=fontsize_internal, &
                                 ytick_labelsize=fontsize_internal, &
                                 tight_layout=tight_layout)
    end subroutine initialize

    subroutine add_plot(self, x, f, label, linestyle, linewidth, markersize, &
                        xscale, yscale, xlim, ylim)
        class(myplot), intent(inout) :: self
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(in) :: f
        character(len=*), intent(in) :: label, linestyle

        integer, intent(in), optional :: linewidth
        integer, intent(in), optional :: markersize
        character(len=*), intent(in), optional :: xscale, yscale
        real(dp), dimension(2), optional :: xlim, ylim

        integer :: my_linewidth

        if (present(linewidth)) then
            my_linewidth = linewidth
        else
            my_linewidth = 1
        end if
        call self%plt%add_plot(x, f, label, linestyle, &
                               linewidth=my_linewidth, &
                               markersize=markersize, &
                               xscale=xscale, &
                               yscale=yscale, &
                               xlim=xlim, &
                               ylim=ylim)
    end subroutine add_plot

    subroutine add_contour(self, x, y, f, levels, colorbar, filled)
        use utils, only: linspace
        class(myplot), intent(inout) :: self
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), dimension(:, :), intent(in) :: f
        integer, optional :: levels
        logical, optional :: colorbar
        logical, optional :: filled

        real(dp), dimension(:), allocatable :: levels_pyplot

        if (present(levels)) then
            allocate (levels_pyplot(levels))
            call linspace(minval(f), maxval(f), levels, levels_pyplot)
            call self%plt%add_contour(x, y, f, linestyle="-", &
                                      levels=levels_pyplot, &
                                      colorbar=colorbar, &
                                      filled=filled)
        else
            call self%plt%add_contour(x, y, f, linestyle="-", &
                                      colorbar=colorbar, &
                                      filled=filled)
        end if
    end subroutine add_contour

    subroutine show(self)
        class(myplot), intent(inout) :: self

        call self%plt%showfig()
    end subroutine show

end module myplot_module
