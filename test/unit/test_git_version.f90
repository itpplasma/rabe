program test_git_version
    use netcdf_mod, only: netcdf_t
    use git_version, only: git_hash

    implicit none

    type(netcdf_t) :: nc_out, nc_in
    character(len=*), parameter :: test_file = "test_git_version.nc"
    character(len=100) :: read_git_hash

    logical :: file_exists
    integer :: unit, iostat

    call nc_out%create(test_file)
    call nc_out%add_global_attribute("git_hash", git_hash)
    call nc_out%close()

    inquire (file=test_file, exist=file_exists)
    if (.not. file_exists) then
        print *, "-------------------------------------------------------------"
        print *, "test_git_version failed: file does not exist after creation"
        error stop
    end if

    call nc_in%open(test_file)
    call nc_in%read_global_attribute("git_hash", read_git_hash)
    call nc_in%close()

    if (len_trim(git_hash) /= 40 .and. trim(git_hash) /= "unknown") then
        print *, "-------------------------------------------------------------"
        print *, "test_git_version failed: git_hash should be 40 chars or 'unknown'"
        print *, "found: '", trim(git_hash), "' with length ", len_trim(git_hash)
        error stop
    end if

    if (trim(read_git_hash) /= trim(git_hash)) then
        print *, "-------------------------------------------------------------"
        print *, "test_git_version failed: git_hash mismatch"
        print *, "found: ", trim(read_git_hash)
        print *, "expected: ", trim(git_hash)
        error stop
    end if

    open (newunit=unit, file=test_file, status="old", iostat=iostat)
    if (iostat == 0) then
        close (unit, status="delete")
    else
        print *, "-------------------------------------------------------------"
        print *, "test_git_version failed: could not clean up test file"
        error stop
    end if

    print *, "test_git_version passed"

end program test_git_version
