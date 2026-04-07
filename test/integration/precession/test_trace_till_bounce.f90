program test_trace_till_bounce
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use utils, only: not_same
    use field_instance, only: initialize_field_instance
    use bounce, only: trace_orbit_till_bounce
    use params, only: params_init

    implicit none

    real(dp), parameter :: reltol = 1.4e-2, abstol = 1e-10
    real(dp), parameter :: abstol_for_zero = 1e-4
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield

    real(dp), dimension(5) :: z_start, z_end
    logical :: test_failed

    call bfield%boozer_field_init(nc_filename, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=6)

    test_failed = .false.

    call initialize_field_instance(bfield)
    call params_init(nfperiods=bfield%nfp, rmajor=bfield%R)

    z_start = [0.5_dp, 0.0_dp, 0.314_dp, 1.0_dp, -0.1_dp]
    call trace_orbit_till_bounce(z_start, z_end)
    print *, "z_end at bounce:", z_end
    z_start = z_end
    call trace_orbit_till_bounce(z_start, z_end)
    print *, "z_end at bounce:", z_end

    if (test_failed) error stop

end program test_trace_till_bounce
