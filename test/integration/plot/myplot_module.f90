module myplot_module
    use pyplot_module, only: pyplot
    use constants, only: dp

    implicit none

    type :: myplot
        type(pyplot) :: plt
    contains
        procedure :: initialize
        procedure :: add_plot
        procedure :: show
    end type myplot

contains

    subroutine initialize(self)
        class(myplot), intent(out) :: self

        call self%plt%initialize()
    end subroutine initialize

    subroutine add_plot(self, x, f, label, linestyle)
        class(myplot), intent(inout) :: self
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x)), intent(in) :: f
        character(len=*), intent(in) :: label, linestyle

        call self%plt%add_plot(x, f, label, linestyle)
    end subroutine add_plot

    subroutine show(self)
        class(myplot), intent(inout) :: self

        call self%plt%showfig()
    end subroutine show

end module myplot_module
