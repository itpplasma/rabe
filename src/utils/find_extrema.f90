module find_extrema
    use constants, only: dp
    implicit none

    abstract interface
        subroutine func1d(x, value)
            use constants, only: dp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(:), intent(out) :: value
        end subroutine func1d
    end interface

contains

    subroutine find_local_minima(func, interval, location, achieved_error)
        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:), allocatable, intent(out) :: location
        real(dp), dimension(:), allocatable, intent(out), optional :: achieved_error

        call find_local_maxima(negative_func, interval, location, achieved_error)

    contains
        subroutine negative_func(x, value)
            use constants, only: dp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(:), intent(out) :: value
            call func(x, value)
            value = -value
        end subroutine negative_func

    end subroutine find_local_minima

    subroutine find_local_maxima(func, interval, location, achieved_error)
        use utils, only: linspace, not_same

        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), dimension(:), allocatable, intent(out) :: location
        real(dp), dimension(:), allocatable, intent(out), optional :: achieved_error

        integer, parameter :: n_coarse = 1001
        integer, parameter :: n_refine = 11
        real(dp) :: x(n_coarse), value(n_coarse)
        real(dp), dimension(:, :), allocatable :: brackets
        real(dp), dimension(:), allocatable :: errors
        integer :: i, n_maxima

        call linspace(interval(1), interval(2), n_coarse, x)
        call func(x, value)

        n_maxima = 0
        do i = 2, n_coarse - 1
            if (value(i - 1) <= value(i) .and. value(i + 1) <= value(i)) &
                n_maxima = n_maxima + 1
        end do

        allocate (brackets(2, n_maxima))
        allocate (location(n_maxima))
        allocate (errors(n_maxima))

        n_maxima = 0
        do i = 2, n_coarse - 1
            if (value(i - 1) <= value(i) .and. value(i + 1) <= value(i)) then
                n_maxima = n_maxima + 1
                brackets(1, n_maxima) = x(i - 1)
                brackets(2, n_maxima) = x(i + 1)
            end if
        end do

        do i = 1, n_maxima
            call refine_maximum(brackets(:, i), location(i), errors(i))
        end do

        if (present(achieved_error)) achieved_error = errors

    contains

        subroutine refine_maximum(bracket, loc, err)
            real(dp), intent(in) :: bracket(2)
            real(dp), intent(out) :: loc, err

            real(dp) :: x_sub(n_refine), val_sub(n_refine)
            real(dp) :: machine_spacing
            real(dp) :: subinterval(2)
            integer :: i_max

            subinterval = bracket
            err = subinterval(2) - subinterval(1)
            do
                call linspace(subinterval(1), subinterval(2), n_refine, x_sub)
                call func(x_sub, val_sub)
                i_max = maxloc(val_sub, dim=1)
                machine_spacing = spacing(max(abs(subinterval(1)), abs(subinterval(2))))
                if (i_max == 1) then
                    err = x_sub(2) - x_sub(1)
                    if (cannot_resolve_at_edge(val_sub(1:2))) exit
                    subinterval = [x_sub(1), x_sub(2)]
                else if (i_max == n_refine) then
                    err = x_sub(n_refine) - x_sub(n_refine - 1)
                    if (cannot_resolve_at_edge(val_sub(n_refine - 1:n_refine))) exit
                    subinterval = [x_sub(n_refine - 1), x_sub(n_refine)]
                else
                    err = x_sub(i_max + 1) - x_sub(i_max - 1)
                    if (cannot_resolve(val_sub(i_max - 1:i_max + 1))) exit
                    subinterval = [x_sub(i_max - 1), x_sub(i_max + 1)]
                end if
                if (err <= machine_spacing) exit
            end do
            loc = x_sub(i_max)
        end subroutine refine_maximum

    end subroutine find_local_maxima

    function cannot_resolve(values)
        use utils, only: not_same
        real(dp), dimension(3), intent(in) :: values
        logical :: cannot_resolve

        logical :: can_resolve
        real(dp) :: tol

        tol = epsilon(values(1))*2.0_dp
        can_resolve = not_same(values(1), values(2), abstol_in=tol) .and. &
                      not_same(values(3), values(2), abstol_in=tol)
        cannot_resolve = .not. can_resolve

    end function cannot_resolve

    function cannot_resolve_at_edge(values)
        use utils, only: not_same
        real(dp), dimension(2), intent(in) :: values
        logical :: cannot_resolve_at_edge

        logical :: can_resolve
        real(dp) :: tol

        tol = epsilon(values(1))*2.0_dp
        can_resolve = not_same(values(1), values(2), abstol_in=tol)
        cannot_resolve_at_edge = .not. can_resolve

    end function cannot_resolve_at_edge

    function find_global_extrema(func, interval, abstol) result(extrema)
        use utils, only: linspace

        procedure(func1d) :: func
        real(dp), intent(in) :: interval(2)
        real(dp), intent(in), optional :: abstol
        real(dp) :: extrema(2)

        integer :: n_steps
        real(dp), dimension(:), allocatable :: x, value

        if (present(abstol)) then
            n_steps = int(abs(interval(2) - interval(1))/abstol) + 2
        else
            n_steps = 1000
        end if

        allocate (x(n_steps), value(n_steps))
        call linspace(interval(1), interval(2), n_steps, x)
        call func(x, value)
        extrema(1) = minval(value)
        extrema(2) = maxval(value)
    end function find_global_extrema

end module find_extrema
