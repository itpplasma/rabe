module constants
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197169399_dp
    real(dp), parameter :: machine_eps = epsilon(1.0_dp)
end module constants
