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

    function is_same(array_1, array_2, reltol_in, abstol_in)
        real(dp), dimension(:), intent(in) :: array_1, array_2
        real(dp), intent(in), optional :: reltol_in, abstol_in
        logical :: is_same

        real(dp) :: reltol, abstol

        if (present(reltol_in)) then
            reltol = reltol_in
        else
            reltol = 0.0_dp
        end if
        if (present(abstol_in)) then
            abstol = abstol_in
        else
            abstol = reltol
        end if

        is_same = any(abs(array_1 - array_2) > reltol*array_1 + abstol)
    end function is_same

end module utils
