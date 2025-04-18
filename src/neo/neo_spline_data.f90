module neo_spline_data
    ! Splines along s
    use nrtype, only: dp, I4B
    real(kind=dp), dimension(:, :), allocatable :: a_bmnc, b_bmnc, c_bmnc, d_bmnc

    real(kind=dp), dimension(:), allocatable :: a_iota, b_iota
    real(kind=dp), dimension(:), allocatable :: c_iota, d_iota
    real(kind=dp), dimension(:), allocatable :: a_curr_tor, b_curr_tor
    real(kind=dp), dimension(:), allocatable :: c_curr_tor, d_curr_tor
    real(kind=dp), dimension(:), allocatable :: a_curr_pol, b_curr_pol
    real(kind=dp), dimension(:), allocatable :: c_curr_pol, d_curr_pol

    real(kind=dp), dimension(:), allocatable :: r_m, r_mhalf
    integer(I4B), dimension(:), allocatable :: sp_index

    logical, save :: lsw_linear_boozer

end module neo_spline_data
