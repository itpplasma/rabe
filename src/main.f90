program rabe
use, intrinsic :: iso_fortran_env, only: dp => real64
use read_file, only: read_field_file, read_boozer_file

call printer(read_boozer_file, "test.bc")

contains

subroutine printer(reader, field_file)
    procedure(read_field_file)  :: reader
    character(len=*) :: field_file

    character(len=20) :: message
    real(dp) :: B 

    message = "Starting rabe!"
    print *, message
    call reader(field_file, B)
    print *, B
end subroutine printer

end program rabe
