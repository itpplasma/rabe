
module field_instance
    use constants, only: dp
    use boozer_field, only: boozer_field_t

    implicit none

    private
    type(boozer_field_t) :: field
    logical :: field_initialized = .false.

    public :: magfie

contains

    subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), dimension(3), intent(in) :: x
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), dimension(3), intent(out) :: bder, hcovar, hctrvr, hcurl

        if (.not. field_initialized) then
            print *, "Error: magfie called before field is initialized."
            error stop
        end if
        call field%evaluate(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end subroutine magfie

    subroutine initialize_field_instance(field_in)
        type(boozer_field_t), intent(in) :: field_in
        if (.not. field_initialized) then
            field = field_in
            field_initialized = .true.
        else
            error stop "Error: initialize_field_instance called more than once."
        end if
    end subroutine initialize_field_instance

end module field_instance
