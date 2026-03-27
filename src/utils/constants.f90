module constants
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197169399_dp
    real(dp), parameter :: eps = epsilon(pi)
end module constants
