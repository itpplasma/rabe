module logger
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

    type(log_levels_t), parameter :: LOG = log_levels_t()
    type(log_level_t) :: log_level = LOG%INFO

    integer, parameter :: stdout_unit = 6
    integer :: log_unit = stdout_unit
    logical :: unsafe_mode = .false.

    integer, parameter :: BUFFER_SIZE = 100
    character(len=256) :: buffer(BUFFER_SIZE)
    integer :: buf_count = 0

    integer, parameter :: MAX_PROBE_FILES = 32
    character(len=256) :: probe_names(MAX_PROBE_FILES)
    integer :: probe_units(MAX_PROBE_FILES)
    integer :: n_probes = 0

contains

    subroutine log_init(unit, level)
        integer, intent(in), optional :: unit
        character(len=*), intent(in), optional :: level

        if (present(unit)) log_unit = unit
        if (present(level)) then
            select case (trim(level))
            case (LOG%DEBUG%name); log_level = LOG%DEBUG
            case (LOG%INFO%name); log_level = LOG%INFO
            case (LOG%WARN%name); log_level = LOG%WARN
            case (LOG%ERROR%name); log_level = LOG%ERROR
            case default
                print *, "Unknown log level: ", trim(level)
                error stop
            end select
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

        buf_count = buf_count + 1
        buffer(buf_count) = entry

        if (level%value >= LOG%ERROR%value .or. buf_count >= BUFFER_SIZE) then
            call log_flush()
        end if
    end subroutine log_msg

    subroutine log_flush()
        integer :: i
        do i = 1, buf_count
            write (log_unit, '(A)') trim(buffer(i))
            if (log_unit /= 6) write (6, '(A)') trim(buffer(i))
        end do
        buf_count = 0
    end subroutine log_flush

    subroutine log_finalize()
        integer :: i
        call log_flush()
        if (log_unit /= 6) close (log_unit)
        do i = 1, n_probes
            close (probe_units(i))
        end do
    end subroutine log_finalize

    subroutine error_stop(msg)
        character(len=*), intent(in) :: msg
        call log_flush()
        error stop msg
    end subroutine error_stop

    subroutine debug_probe(value, filename)
        real, intent(in) :: value
        character(len=*), intent(in) :: filename
        integer :: unit, i
        if (LOG%DEBUG%value < log_level%value) return
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
