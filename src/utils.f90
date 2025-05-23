module utils
    use constants, only: dp

    implicit none

    interface is_same
        module procedure is_same_scalar
        module procedure is_same_array
    end interface

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

    function is_same_scalar(scalar_1, scalar_2, reltol_in, abstol_in)
        use ieee_arithmetic, only: ieee_is_nan
        real(dp), intent(in) :: scalar_1, scalar_2
        real(dp), intent(in), optional :: reltol_in, abstol_in
        logical :: is_same_scalar

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

        if (ieee_is_nan(scalar_1) .or. ieee_is_nan(scalar_2)) then
            is_same_scalar = .true.
        else
            is_same_scalar = abs(scalar_1 - scalar_2) > (reltol*abs(scalar_1) + abstol)
        end if
    end function is_same_scalar

    function is_same_array(array_1, array_2, reltol_in, abstol_in)
        real(dp), dimension(:), intent(in) :: array_1, array_2
        real(dp), intent(in), optional :: reltol_in, abstol_in
        logical :: is_same_array

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

        if (contains_nan(array_1) .or. contains_nan(array_2)) then
            is_same_array = .true.
        else
            is_same_array = any(abs(array_1 - array_2) > (reltol*abs(array_1) + abstol))
        end if
    end function is_same_array

    function contains_nan(array)
        use ieee_arithmetic, only: ieee_is_nan
        implicit none
        real(dp), intent(in) :: array(:)
        logical :: contains_nan
        integer :: i

        contains_nan = .false.
        do i = 1, size(array)
            if (ieee_is_nan(array(i))) then
                contains_nan = .true.
                return
            end if
        end do
    end function contains_nan

end module utils
