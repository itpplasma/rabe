program rabe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use read_file, only: read_field_file, read_boozer_file
    use read_file, only: modes
    use neo_magfie, only: neo_magfie_a

    implicit none

    real(dp), dimension(3) :: x = (/0.5d0, 0.0d0, 0.d0/)
    real(dp) :: bmod, sqrtg, dB_dx(3)

    call neo_magfie_a(x, bmod, sqrtg, dB_dx)
    print *, x
    print *, "bmod"
    print *, bmod
    print *, "sqrtg"
    print *, sqrtg
    print *, "bder"
    print *, dB_dx

contains

    subroutine printer(reader, field_file)
        procedure(read_field_file) :: reader
        character(len=*) :: field_file

        character(len=20) :: message
        type(modes) :: B

        message = "Starting rabe!"
        print *, message
        call reader(field_file, B)
        print *, B%s_tor(1:3)
        print *, B%coef(1:3, 1:3)
        print *, B%m(1:3)
        print *, B%n(1:3)
    end subroutine printer

end program rabe
