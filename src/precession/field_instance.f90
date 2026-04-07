
module field_instance
    use constants, only: dp
    use boozer_field, only: boozer_field_t

    implicit none

    private
    type(boozer_field_t) :: field
    logical :: field_initialized = .false.
    real(dp), parameter :: cm2m = 1e-2, gauss2tesla = 1e-4

    public :: magfie, initialize_field_instance

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
        bmod = bmod/gauss2tesla
        sqrtg = sqrtg/cm2m**3.0_dp
        hcovar = hcovar/cm2m
        hctrvr = hctrvr*cm2m
        hcurl = hcurl*cm2m**2.0_dp
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
