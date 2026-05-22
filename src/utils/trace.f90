module trace_mod
    use constants, only: dp
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    character(len=:), allocatable, save :: trace_dir
    logical, save :: trace_enabled = .false.
    logical, save :: trace_initialized = .false.
    integer, save :: surface_id = 0

    public :: trace_init
    public :: trace_active
    public :: trace_set_surface
    public :: trace_scalar
    public :: trace_int
    public :: trace_real_1d
    public :: trace_int_1d
    public :: trace_real_2d
    public :: trace_surface_open

contains

    subroutine trace_init()
        character(len=512) :: buf
        integer :: stat, n

        if (trace_initialized) return
        call get_environment_variable("RABE_TRACE_DIR", buf, length=n, status=stat)
        if (stat == 0 .and. n > 0) then
            trace_dir = trim(buf(1:n))
            trace_enabled = .true.
            call execute_command_line("mkdir -p '"//trace_dir//"'", wait=.true.)
        end if
        trace_initialized = .true.
    end subroutine trace_init

    function trace_active() result(active)
        logical :: active
        if (.not. trace_initialized) call trace_init()
        active = trace_enabled
    end function trace_active

    subroutine trace_set_surface(idx)
        integer, intent(in) :: idx
        surface_id = idx
    end subroutine trace_set_surface

    function tag_path(tag) result(path)
        character(len=*), intent(in) :: tag
        character(len=:), allocatable :: path
        character(len=32) :: buf
        write (buf, "(I6.6)") surface_id
        path = trim(trace_dir)//"/s"//trim(adjustl(buf))//"_"//tag//".txt"
    end function tag_path

    subroutine open_tag(tag, unit)
        character(len=*), intent(in) :: tag
        integer, intent(out) :: unit
        integer :: ios
        open (newunit=unit, file=tag_path(tag), status="replace", action="write", &
              iostat=ios)
        if (ios /= 0) then
            print *, "trace_mod: failed to open ", tag_path(tag)
            error stop
        end if
    end subroutine open_tag

    subroutine trace_surface_open(idx)
        integer, intent(in) :: idx
        call trace_set_surface(idx)
    end subroutine trace_surface_open

    subroutine trace_scalar(name, value)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: value
        integer :: u
        if (.not. trace_active()) return
        call open_tag(name, u)
        write (u, "(ES23.15)") value
        close (u)
    end subroutine trace_scalar

    subroutine trace_int(name, value)
        character(len=*), intent(in) :: name
        integer, intent(in) :: value
        integer :: u
        if (.not. trace_active()) return
        call open_tag(name, u)
        write (u, "(I12)") value
        close (u)
    end subroutine trace_int

    subroutine trace_real_1d(name, arr)
        character(len=*), intent(in) :: name
        real(dp), dimension(:), intent(in) :: arr
        integer :: u, i
        if (.not. trace_active()) return
        call open_tag(name, u)
        do i = 1, size(arr)
            write (u, "(ES23.15)") arr(i)
        end do
        close (u)
    end subroutine trace_real_1d

    subroutine trace_int_1d(name, arr)
        character(len=*), intent(in) :: name
        integer, dimension(:), intent(in) :: arr
        integer :: u, i
        if (.not. trace_active()) return
        call open_tag(name, u)
        do i = 1, size(arr)
            write (u, "(I12)") arr(i)
        end do
        close (u)
    end subroutine trace_int_1d

    subroutine trace_real_2d(name, arr)
        character(len=*), intent(in) :: name
        real(dp), dimension(:, :), intent(in) :: arr
        integer :: u, i, j
        if (.not. trace_active()) return
        call open_tag(name, u)
        do j = 1, size(arr, 2)
            do i = 1, size(arr, 1)
                write (u, "(ES23.15)") arr(i, j)
            end do
        end do
        close (u)
    end subroutine trace_real_2d

end module trace_mod
