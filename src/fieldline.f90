module fieldline_mod
    use constants, only: dp, pi

    implicit none

    type :: fieldline_t
        real(dp) :: theta_0
        real(dp) :: phi_0
        real(dp) :: iota
        real(dp) :: phi_max(2)
        real(dp) :: B_max(2)
        real(dp) :: iota_p

        real(dp) :: eta_b
        real(dp) :: delta_eta
        real(dp) :: integral_lambda_b_over_B_squared
        real(dp) :: delta_aspect_ratio
        real(dp) :: integral_one_over_B_squared
        real(dp) :: radial_drift
    contains
        generic :: get_theta => get_theta_scalar, get_theta_array
        procedure, private :: get_theta_scalar
        procedure, private :: get_theta_array
    end type fieldline_t

contains

    function get_theta_scalar(self, phi) result(theta)
        class(fieldline_t), intent(in) :: self
        real(dp) :: phi

        real(dp) :: theta

        theta = (phi - self%phi_0)*self%iota + self%theta_0
        theta = modulo(theta, 2.0_dp*pi)
    end function get_theta_scalar

    function get_theta_array(self, phi) result(theta)
        class(fieldline_t), intent(in) :: self
        real(dp), dimension(:) :: phi

        real(dp), dimension(size(phi)) :: theta

        theta = (phi - self%phi_0)*self%iota + self%theta_0
        theta = modulo(theta, 2.0_dp*pi)
    end function get_theta_array

end module fieldline_mod
