program plot_nabla_s_compare
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use utils, only: linspace
    use boozer_field, only: boozer_field_t
    use neo_field, only: neo_field_t

    implicit none

    character(len=*), parameter :: vmec_file = &
                      "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    character(len=*), parameter :: bc_file = &
                                   "input/landreman_paul_qh_flipped_vmns.bc"

    type(boozer_field_t) :: bfield
    type(neo_field_t) :: nfield

    real(dp), parameter :: stor = 0.5_dp
    integer, parameter :: n_theta = 100, n_phi = 100
    real(dp), dimension(n_theta) :: theta
    real(dp), dimension(n_phi) :: phi
    real(dp), dimension(n_theta, n_phi) :: rel_diff

    real(dp) :: nabla_s_b, nabla_s_n
    type(myplot) :: plt
    integer :: i, j

    call bfield%boozer_field_init(vmec_file, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=6)
    call bfield%fix_to_surface(stor)

    call nfield%neo_field_init(bc_file, stor)

    call linspace(0.0_dp, 2.0_dp*pi, n_theta, theta)
    call linspace(0.0_dp, 2.0_dp*pi, n_phi, phi)

    do j = 1, n_phi
        do i = 1, n_theta
            call bfield%compute_nabla_s(theta(i), phi(j), nabla_s_b)
            call nfield%compute_nabla_s(theta(i), phi(j), nabla_s_n)
            rel_diff(i, j) = (nabla_s_b - nabla_s_n)/nabla_s_n
        end do
    end do

    call plt%initialize( &
        xlabel="$\vartheta$ [$\pi$]", &
        ylabel="$\varphi$ [$\pi$]", &
        title="Relative difference of $\nabla s$: Landreman-Paul QH")
    call plt%add_contour(theta/pi, phi/pi, rel_diff, &
                         levels=20, colorbar=.true., filled=.true.)
    call plt%show()

end program plot_nabla_s_compare
