module utils
    use constants, only: dp

    implicit none

contains
    subroutine linspace(a, b, n, x)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp), dimension(:), intent(out) :: x

        real(dp) :: dx
        integer :: i

        dx = (b - a)/(n - 1)
        do i = 1, n
            x(i) = a + (i - 1)*dx
        end do
    end subroutine linspace
end module utils
