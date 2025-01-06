module read_file
use, intrinsic :: iso_fortran_env, only: dp => real64


    abstract interface
        subroutine read_field_file(file_name, B)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            character(len=256), intent(in) :: file_name
            real(dp), intent(out) :: B 
        end subroutine read_field_file
    end interface


    contains


        subroutine read_boozer_file(file_name, B)
            character(len=256), intent(in) :: file_name
            real(dp), intent(out) :: B 
            
            integer :: file_unit

            open(newunit=file_unit, file=file_name)
            read(file_unit, *) B
            close(file_unit)
        end subroutine read_boozer_file

end module read_file
