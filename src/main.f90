program redl

    call printer()

    contains

        subroutine printer()
            character(len=20) :: message

            message = "Starting redl!"
            print *, message
        end subroutine printer

end program redl
