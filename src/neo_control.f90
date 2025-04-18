!> Control parameters from input file
module neo_control
    use nrtype, only: dp
    character(40) :: in_file
    integer :: theta_n
    integer :: phi_n
    integer :: s_ind_in
    integer :: write_progress
    integer :: write_output_files
    integer :: calc_fourier
    integer :: spline_test
    integer :: max_m_mode, max_n_mode
    integer :: lab_swi, inp_swi, ref_swi, eout_swi
    integer :: chk_swi
    integer :: fluxs_interp
    integer :: s_num
    real(kind=dp) :: s_start, s_end
    integer :: g11_swi
    integer :: eval_mode
    integer :: no_fluxs, no_fluxs_s
    integer, dimension(:), allocatable :: fluxs_arr

    !> Controls setting of rt0 and bmref in neo_sub.
    logical :: set_rt0_from_rmnc_for_zero_mode = .true.
end module neo_control
