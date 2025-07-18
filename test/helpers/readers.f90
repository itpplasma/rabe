module readers
    use constants, only: dp
    implicit none

contains
    subroutine read_column(filename, data, n_columns, select_column)
        character(len=*), intent(in) :: filename
        real(dp), dimension(:), allocatable, intent(out) :: data
        integer, intent(in) :: n_columns
        integer, intent(in), optional :: select_column

        integer :: i, n, unit, column, ios
        real(dp), dimension(:), allocatable :: row
        character(len=1024) :: line

        column = 1
        if (present(select_column)) column = select_column

        open (newunit=unit, file=filename, status='old', action='read')
        n = 0
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            n = n + 1
        end do
        rewind (unit)

        read (unit, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, 'Error reading file.'
            stop
        end if

        if (column > n_columns) then
          print *, 'Requested column ', column, ' exceeds number of columns ', n_columns
            error stop
        end if
        rewind (unit)

        allocate (data(n))
        allocate (row(n_columns))
        do i = 1, n
            read (unit, *, iostat=ios) row
            if (ios /= 0) then
                print *, 'Error reading line ', i
                error stop
            end if
            data(i) = row(column)
        end do

        close (unit)
    end subroutine read_column

end module readers
