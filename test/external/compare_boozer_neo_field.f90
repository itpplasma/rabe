program compare_boozer_neo_field
    use constants, only: dp, pi
    use utils, only: linspace
    use boozer_field, only: boozer_field_t
    use neo_field, only: neo_field_t
    use myplot_module, only: myplot

    implicit none

    integer, parameter :: n_theta = 100, n_phi = 100
    real(dp), parameter :: stor = 0.397959184_dp

    character(len=*), parameter :: nc_filename = &
        "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
    character(len=*), parameter :: bc_filename = &
        "input/landreman_paul_qa.bc"

    type(boozer_field_t) :: bfield
    type(neo_field_t) :: nfield

    real(dp) :: theta(n_theta), phi(n_phi)
    real(dp) :: B_boozer, B_neo, dB_dx_boozer(3), dB_dx_neo(3)
    real(dp) :: rel_diff_B(n_theta, n_phi)
    real(dp) :: rel_diff_dBdtheta(n_theta, n_phi)
    real(dp) :: dBdtheta_neo(n_theta, n_phi)
    integer :: i, j
    type(myplot) :: plt

    call bfield%boozer_field_init(nc_filename, &
        radial_spline_order=5, &
        angular_spline_order=5, &
        grid_refinement=3)
    call bfield%fix_to_surface(stor)

    call nfield%neo_field_init(bc_filename, stor)

    call linspace(0.0_dp, 2.0_dp*pi, n_theta, theta, &
        include_endpoint=.false.)
    call linspace(0.0_dp, 2.0_dp*pi/bfield%nfp, n_phi, phi, &
        include_endpoint=.false.)

    do j = 1, n_phi
        do i = 1, n_theta
            call bfield%compute_B_and_dB_dx(theta(i), phi(j), &
                B_boozer, dB_dx_boozer)
            call nfield%compute_B_and_dB_dx(theta(i), phi(j), &
                B_neo, dB_dx_neo)
            rel_diff_B(i, j) = (B_boozer - B_neo)/B_neo
            dBdtheta_neo(i, j) = dB_dx_neo(2)
            rel_diff_dBdtheta(i, j) = &
                dB_dx_boozer(2) - dB_dx_neo(2)
        end do
    end do

    call plt%initialize( &
        xlabel="$\varphi$", &
        ylabel="$\vartheta$", &
        title="$(B_{boozer} - B_{neo}) / B_{neo}$")
    call plt%add_contour(phi, theta, rel_diff_B, &
        levels=20, colorbar=.true., filled=.true.)
    call plt%show()

    rel_diff_dBdtheta = rel_diff_dBdtheta &
        /maxval(abs(dBdtheta_neo))

    call plt%initialize( &
        xlabel="$\varphi$", &
        ylabel="$\vartheta$", &
        title= &
        "$(\partial_\vartheta B_{boozer} - \partial_\vartheta B_{neo})" &
        //" / \max|\partial_\vartheta B_{neo}|$")
    call plt%add_contour(phi, theta, rel_diff_dBdtheta, &
        levels=20, colorbar=.true., filled=.true.)
    call plt%show()

end program compare_boozer_neo_field
