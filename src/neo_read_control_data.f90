module sizey_bo
! Definition for rk4d_bo also used in main routine neo
    integer :: npart
    integer :: multra
    integer :: ndim
    integer, parameter :: npq = 7
end module sizey_bo

module sizey_cur
! Definition for rk4d_bo also used in main routine neo
    integer :: npart_cur
    integer :: ndim_cur
    integer, parameter :: npq_cur = 11
    integer :: alpha_cur
end module sizey_cur

module sizey_pla
    use nrtype, only: dp
    ! Definition for rk4d_pla
    integer :: npart_pla
    integer :: ndim_pla
    integer, parameter :: npq_pla = 3
    real(kind=dp) :: lamup_pla
    real(kind=dp) :: lambda_alpha
    real(kind=dp) :: nufac_pla
end module sizey_pla

module neo_van
    use nrtype, only: dp
    real(kind=dp) :: v_phi0, v_theta0
    real(kind=dp) :: bmin_tol
    integer :: v_nper, v_steps
    integer :: v_num_mm
    integer :: no_minima
    integer, dimension(:), allocatable :: li_minima
    integer :: no_gamma
    integer :: tau_num
    integer :: tau_max_iter
    real(kind=dp) :: lambda_fac
    real(kind=dp) :: temp_e
    real(kind=dp) :: gamma_eps
    real(kind=dp) :: phi_eps
end module neo_van
