module neo_input
    ! Input from data files (Boozer)
    use nrtype, only: dp
    integer, dimension(:), allocatable :: ixm, ixn
    integer, dimension(:), allocatable :: pixm, pixn
    integer, dimension(:), allocatable :: i_m, i_n

    real(kind=dp), dimension(:), allocatable :: es
    real(kind=dp), dimension(:), allocatable :: pprime
    real(kind=dp), dimension(:), allocatable :: sqrtg00

    real(kind=dp), dimension(:, :), allocatable :: rmnc, zmns, vmns
    real(kind=dp), dimension(:, :), allocatable :: bmnc
    real(kind=dp), dimension(:), allocatable :: b00

    real(kind=dp) :: flux, psi_pr

    integer :: m0b, n0b
    integer :: ns, mnmax, nfp
    integer :: m_max, n_max
end module neo_input
