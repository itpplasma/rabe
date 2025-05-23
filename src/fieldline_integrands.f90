module fieldline_integrands
    use constants, only: dp
    use field_base, only: field_t

    implicit none

contains

    function radial_drift_velocity(field, theta, phi, eta)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi, eta
        real(dp) :: radial_drift_velocity

        real(dp) :: B, dB_dx(3), lambda_squared, dB_dtheta_0

        call field%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        lambda_squared = pitchparameter_squared(B, eta)

        dB_dtheta_0 = dB_dx(2)

        radial_drift_velocity = 0.5_dp*sqrt(lambda_squared)/B**3*(3 + lambda_squared) &
                                *dB_dtheta_0
    end function radial_drift_velocity

    function lambda_over_B_squared(field, theta, phi, eta)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi, eta
        real(dp) :: lambda_over_B_squared

        real(dp) :: B

        call field%compute_B_mod(theta, phi, B)
        lambda_over_B_squared = sqrt(pitchparameter_squared(B, eta))/B**2
    end function lambda_over_B_squared

    function B_squared(field, theta, phi)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi
        real(dp) :: B_squared

        real(dp) :: B

        call field%compute_B_mod(theta, phi, B)
        B_squared = B**2
    end function B_squared

    function pitchparameter_squared(B, eta)
        real(dp), intent(in) :: B, eta
        real(dp) :: pitchparameter_squared

        pitchparameter_squared = 1 - B*eta

        if (pitchparameter_squared .lt. 0.0_dp) then
            print *, "Square of pitch parameter (1 - B*eta) is negative!"
            print *, "B: ", B
            print *, "eta: ", eta
            error stop
        end if
    end function pitchparameter_squared

end module fieldline_integrands
