program plot_I_ref_compare
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use neo_field, only: neo_field_t
    use fieldline_mod, only: fieldline_t
    use fieldline_labels, only: get_labels
    use make_fieldline, only: make_flock_of_fieldlines

    implicit none

    character(len=*), parameter :: vmec_file = &
                      "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    character(len=*), parameter :: bc_file = &
                                   "input/landreman_paul_qh_flipped_vmns.bc"

    type(boozer_field_t) :: bfield
    type(neo_field_t) :: nfield

    real(dp), parameter :: stor = 0.5_dp
    real(dp), parameter :: M_pol = -1.0_dp, N_tor = 4.0_dp
    integer, parameter :: max_n_fieldlines = 100

    type(fieldline_t), dimension(:), allocatable :: fl_b, fl_n
    real(dp), dimension(:), allocatable :: xi_0
    real(dp) :: iota, approx_iota, B_theta, B_phi
    real(dp), dimension(2) :: x_range, I_ref_line
    type(myplot) :: plt

    call bfield%boozer_field_init(vmec_file, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=6)
    call bfield%fix_to_surface(stor)
    call bfield%get_iota_and_covariant_components(stor, iota, B_theta, B_phi)
    call get_labels(max_n_fieldlines, iota, M_pol, N_tor, bfield%nfp, &
                    xi_0, approx_iota)
    allocate (fl_b(size(xi_0)))
    call make_flock_of_fieldlines(fl_b, xi_0, approx_iota, bfield, &
                                  M_pol, N_tor, bfield%nfp)

    call nfield%neo_field_init(bc_file, stor)
    iota = nfield%iota
    call get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfield%nfp, &
                    xi_0, approx_iota)
    allocate (fl_n(size(xi_0)))
    call make_flock_of_fieldlines(fl_n, xi_0, approx_iota, nfield, &
                                  M_pol, N_tor, nfield%nfp)

    call plt%initialize( &
        xlabel="$\xi_0$ [$\pi$]", &
        ylabel="$I$", &
        legend=.true., &
        title="Boundary layer width I")
    call plt%add_plot(fl_b%xi_0/pi, &
                      fl_b%integral_lambda_b_over_B_squared, &
                      label="boozer_field", linestyle="bo")
    call plt%add_plot(fl_n%xi_0/pi, &
                      fl_n%integral_lambda_b_over_B_squared, &
                      label="neo_field", linestyle="rx")

    x_range = [0.0_dp, 2.0_dp]
    I_ref_line = fl_b(1)%I_ref
    call plt%add_plot(x_range, I_ref_line, &
                      label="$I_\mathrm{ref}$ boozer", linestyle="b--")
    I_ref_line = fl_n(1)%I_ref
    call plt%add_plot(x_range, I_ref_line, &
                      label="$I_\mathrm{ref}$ neo", linestyle="r--")
    call plt%show()

    print *, (nfield%R - bfield%R)/bfield%R

end program plot_I_ref_compare
