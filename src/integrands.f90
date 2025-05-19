module integrands
    use constants, only: dp
    use field_base, only: field_t

    implicit none

contains

    function B_mod_squared(field, theta, phi)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi
        real(dp) :: B_mod_squared

        real(dp) :: B_mod

        call field%compute_B_mod(theta, phi, B_mod)
        B_mod_squared = B_mod**2
    end function B_mod_squared

    function lambda_squared(field, theta, phi, eta)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi, eta
        real(dp) :: lambda_squared

        real(dp) :: B_mod

        call field%compute_B_mod(theta, phi, B_mod)
        lambda_squared = 1 - B_mod*eta

        if (lambda_squared .lt. 0.0_dp) then
            print *, "Square of pitch parameter is negative!"
            print *, "B_mod: ", B_mod
            print *, "eta: ", eta
            error stop
        end if

    end function lambda_squared

    function radial_drift_velocity(field, theta, phi, eta)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi, eta
        real(dp) :: radial_drift_velocity

        real(dp) :: B_mod, lambda_sqared, dB_dtheta_0

        error stop "radial_drift_veloctiy not yet implemented!"

    end function radial_drift_velocity

end module integrands
