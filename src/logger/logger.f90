module logger
    use constants, only: dp
    implicit none

    type :: log_level_t
        integer, private :: value
        character(len=5), private :: name
    end type

    type :: log_levels_t
        type(log_level_t) :: DEBUG = log_level_t(0, 'DEBUG')
        type(log_level_t) :: INFO = log_level_t(1, 'INFO ')
        type(log_level_t) :: WARN = log_level_t(2, 'WARN ')
        type(log_level_t) :: ERROR = log_level_t(3, 'ERROR')
    end type

    type(log_levels_t), parameter :: log_lvl = log_levels_t()

    interface log_val
        module procedure log_val_real
        module procedure log_val_real_array
        module procedure log_val_int
        module procedure log_val_logical
    end interface log_val

    private
    public :: log_init, log_msg, log_val, log_finalize, error_stop, debug_probe
    public :: log_lvl, log_level_t, log_levels_t

    type(log_level_t) :: log_level = log_level_t(1, 'INFO ')

    integer, parameter :: stdout_unit = 6
    logical, parameter :: unsafe_mode = .false.
    integer :: log_unit = stdout_unit

    integer, parameter :: BUFFER_SIZE = 100
    character(len=256) :: buffer(BUFFER_SIZE)
    integer :: buffer_count = 0

    integer, parameter :: MAX_PROBE_FILES = 32
    character(len=256) :: probe_names(MAX_PROBE_FILES)
    integer :: probe_units(MAX_PROBE_FILES)
    integer :: n_probes = 0

contains

    subroutine log_init(log_file, level_name)
        character(len=*), intent(in), optional :: log_file
        character(len=*), intent(in), optional :: level_name

        if (present(log_file)) then
            open (newunit=log_unit, file=log_file, status="replace", action="write")
        else
            log_unit = stdout_unit
        end if
        if (present(level_name)) then
            select case (trim(level_name))
            case (log_lvl%DEBUG%name); log_level = log_lvl%DEBUG
            case (log_lvl%INFO%name); log_level = log_lvl%INFO
            case (log_lvl%WARN%name); log_level = log_lvl%WARN
            case (log_lvl%ERROR%name); log_level = log_lvl%ERROR
            case default
                print *, "Unknown log level: ", trim(level_name)
                error stop
            end select
        else
            log_level = log_lvl%INFO
        end if
    end subroutine log_init

    subroutine log_msg(level, msg)
        type(log_level_t), intent(in) :: level
        character(len=*), intent(in) :: msg
        character(len=8) :: date
        character(len=10) :: time
        character(len=8) :: unsafe_str
        character(len=256) :: entry

        if (level%value < log_level%value) return

        call date_and_time(date, time)

        unsafe_str = merge("UNSAFE  ", "        ", unsafe_mode)

        write (entry, '(A,1X,A,1X,A,1X,A)') &
            date//'T'//time(1:6), &
            level%name, &
            unsafe_str, &
            trim(msg)

        buffer_count = buffer_count + 1
        buffer(buffer_count) = entry

        if (level%value >= log_lvl%ERROR%value .or. buffer_count >= BUFFER_SIZE) then
            call log_flush()
        end if
    end subroutine log_msg

    subroutine log_val_real(level, label, val)
        type(log_level_t), intent(in) :: level
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: val
        character(len=256) :: buf
        write (buf, '(A,G0)') label, val
        call log_msg(level, trim(buf))
    end subroutine log_val_real

    subroutine log_val_real_array(level, label, val)
        type(log_level_t), intent(in) :: level
        character(len=*), intent(in) :: label
        real(dp), dimension(:), intent(in) :: val
        character(len=256) :: buf
        integer :: i
        do i = 1, size(val)
            write (buf, '(A,"(",I0,"): ",G0)') label, i, val(i)
            call log_msg(level, trim(buf))
        end do
    end subroutine log_val_real_array

    subroutine log_val_int(level, label, val)
        type(log_level_t), intent(in) :: level
        character(len=*), intent(in) :: label
        integer, intent(in) :: val
        character(len=256) :: buf
        write (buf, '(A,I0)') label, val
        call log_msg(level, trim(buf))
    end subroutine log_val_int

    subroutine log_val_logical(level, label, val)
        type(log_level_t), intent(in) :: level
        character(len=*), intent(in) :: label
        logical, intent(in) :: val
        character(len=256) :: buf
        write (buf, '(A,L1)') label, val
        call log_msg(level, trim(buf))
    end subroutine log_val_logical

    subroutine log_flush()
        integer :: i
        do i = 1, buffer_count
            write (log_unit, '(A)') trim(buffer(i))
        end do
        buffer_count = 0
    end subroutine log_flush

    subroutine log_finalize()
        integer :: i
        call log_flush()
        if (log_unit /= stdout_unit) close (log_unit)
        do i = 1, n_probes
            close (probe_units(i))
        end do
    end subroutine log_finalize

    subroutine error_stop(msg)
        character(len=*), intent(in) :: msg
        call log_finalize()
        error stop msg
    end subroutine error_stop

    subroutine debug_probe(value, filename)
        real, intent(in) :: value
        character(len=*), intent(in) :: filename
        integer :: unit, i
        if (log_lvl%DEBUG%value < log_level%value) return
        do i = 1, n_probes
            if (trim(probe_names(i)) == trim(filename)) then
                write (probe_units(i), '(ES20.10)') value
                return
            end if
        end do
        n_probes = n_probes + 1
        open (newunit=unit, file=filename)
        probe_names(n_probes) = filename
        probe_units(n_probes) = unit
        write (unit, '(ES20.10)') value
    end subroutine debug_probe

end module logger
