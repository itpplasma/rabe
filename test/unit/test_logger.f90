program test_logger
    use logger, only: log_init, log_msg, log_finalize, LOG
    implicit none

    character(len=*), parameter :: log_file = "test_logger_smoke.log"
    integer :: unit, ios, n_lines
    character(len=512) :: line
    logical :: found_debug, found_info, found_warn, found_error
    logical :: test_failed

    test_failed = .false.

    ! --- smoke: all four levels appear when configured at DEBUG ---
    call log_init(log_file=log_file, level_name="DEBUG")
    call log_msg(LOG%DEBUG, "debug message")
    call log_msg(LOG%INFO, "info message")
    call log_msg(LOG%WARN, "warn message")
    call log_msg(LOG%ERROR, "error message")
    call log_finalize()

    open (newunit=unit, file=log_file, status="old", action="read")
    found_debug = .false.; found_info = .false.
    found_warn = .false.; found_error = .false.
    n_lines = 0
    do
        read (unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n_lines = n_lines + 1
        if (index(line, "DEBUG") > 0 .and. index(line, "debug message") > 0) &
            found_debug = .true.
        if (index(line, "INFO") > 0 .and. index(line, "info message") > 0) &
            found_info = .true.
        if (index(line, "WARN") > 0 .and. index(line, "warn message") > 0) &
            found_warn = .true.
        if (index(line, "ERROR") > 0 .and. index(line, "error message") > 0) &
            found_error = .true.
    end do
    close (unit)

    if (n_lines /= 4) then
        print *, "smoke test failed: expected 4 lines, got", n_lines
        test_failed = .true.
    end if
    if (.not. found_debug) then
        print *, "smoke test failed: DEBUG line missing or malformed"
        test_failed = .true.
    end if
    if (.not. found_info) then
        print *, "smoke test failed: INFO line missing or malformed"
        test_failed = .true.
    end if
    if (.not. found_warn) then
        print *, "smoke test failed: WARN line missing or malformed"
        test_failed = .true.
    end if
    if (.not. found_error) then
        print *, "smoke test failed: ERROR line missing or malformed"
        test_failed = .true.
    end if

    ! --- level filtering: DEBUG suppressed when configured at INFO ---
    call log_init(log_file=log_file, level_name="INFO")
    call log_msg(LOG%DEBUG, "suppressed debug")
    call log_msg(LOG%INFO, "visible info")
    call log_finalize()

    open (newunit=unit, file=log_file, status="old", action="read")
    n_lines = 0
    do
        read (unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n_lines = n_lines + 1
        if (index(line, "suppressed debug") > 0) then
            print *, "level filtering failed: DEBUG message appeared at INFO level"
            test_failed = .true.
        end if
    end do
    close (unit)

    if (n_lines /= 1) then
        print *, "level filtering failed: expected 1 line, got", n_lines
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_logger
