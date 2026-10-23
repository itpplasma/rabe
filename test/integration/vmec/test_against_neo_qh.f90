!> Test that boozer_field_t produces the same results as the NEO code
!> for a quasi-helical equilibrium.
!> The reference values were obtained from a field evalution of a .bc file
!> converted from the original .nc file
program test_against_neo_qh
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.2e-2, abstol = 1e-10
    real(dp), parameter :: abstol_for_zero = 1e-3
    ! Symmetry-zero dB_dx components (exact 0 in the NEO reference) come out
    ! as spline roundoff, measured <= 1.6e-10 with an optimized libneo build.
    real(dp), parameter :: abstol_roundoff = 1e-9
    ! NEO values below are for the ns101 fixture (wout SHA256 a8fb3ac3...).
    ! They were regenerated with the route that reproduces the original ns50
    ! values to 7.2e-9: booz_xform 0.1.0 (mboz = nboz = 32) on all half-grid
    ! surfaces, written with the libneo boozer.py stellarator-symmetric .bc
    ! writer with vmns = +numns_b*nfp/(2*pi), and evaluated by neo_field_t of
    ! rabe_internal 3b37cf4. The 100-surface .bc has SHA256 a8e7d0d3....
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QH_reactorScale_ns101_reference.nc"

    type(boozer_field_t) :: bfield

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)

    real(dp) :: bmod_b, sqrtg_b, dB_dx_b(3), nabla_s
    real(dp) :: iota_b, B_theta_b, B_phi_b

    ! Neo reference values (from generate_neo_ref_qh)
    real(dp), parameter :: cm2m = 1e-2_dp
    real(dp), parameter :: psi_tor_edge_ref = 6.66284343900593610e+00_dp
    real(dp), parameter :: nfp_ref = 4.0_dp

    real(dp) :: iota_ref(n_cases)
    real(dp) :: B_theta_ref(n_cases), B_phi_ref(n_cases)
    real(dp) :: bmod_ref(n_cases), sqrtg_ref(n_cases)
    real(dp) :: dB_dx_ref(n_cases, 3)
    real(dp) :: nabla_s_ref(n_cases)

    integer :: case
    logical :: test_failed

    call bfield%init_from_vmec(nc_filename, &
                               radial_spline_order=5, &
                               angular_spline_order=5, &
                               grid_refinement=3)

    test_failed = .false.
    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta = [0.0_dp, 1.0_dp, 3.14_dp, 0.5_dp, 2.0_dp]
    phi = [0.0_dp, 0.5_dp, 1.57_dp, 2.0_dp, 3.0_dp]

    iota_ref = [-1.24190205199088588e+00_dp, &
                -1.24282935621631685e+00_dp, &
                -1.24348360775440425e+00_dp, &
                -1.24410842130342236e+00_dp, &
                -1.24469547847326711e+00_dp]

    B_theta_ref = [2.65738812255297193e-16_dp, &
                   -3.20754212570409648e-16_dp, &
                   -6.84778672128444935e-16_dp, &
                   -4.30291922957897257e-17_dp, &
                   -2.18349235602578885e-15_dp]

    B_phi_ref = [1.02238469207864199e+02_dp, &
                 1.02238456374797025e+02_dp, &
                 1.02238451199997868e+02_dp, &
                 1.02238448000010081e+02_dp, &
                 1.02238445599999793e+02_dp]

    bmod_ref = [5.58457421386239261e+00_dp, &
                6.44751295416938319e+00_dp, &
                6.64012823602645064e+00_dp, &
                6.42993108804714275e+00_dp, &
                5.80699374768080379e+00_dp]

    sqrtg_ref = [-2.18420709451071694e+07_dp, &
                 -1.63866223697137162e+07_dp, &
                 -1.54497315248996336e+07_dp, &
                 -1.64763579998579454e+07_dp, &
                 -2.02009191106078215e+07_dp]
    sqrtg_ref = sqrtg_ref*(cm2m**3.0_dp)

    dB_dx_ref(1, :) = [-1.43832229148989721e+00_dp, &
                       0.00000000000000000e+00_dp, &
                       0.00000000000000000e+00_dp]
    dB_dx_ref(2, :) = [1.01144187047813827e+00_dp, &
                       8.28688681473852135e-02_dp, &
                       3.32286861749094331e-01_dp]
    dB_dx_ref(3, :) = [8.23617368255439630e-01_dp, &
                       3.83248854998836062e-03_dp, &
                       1.51487651151810329e-02_dp]
    dB_dx_ref(4, :) = [4.34479001471452964e-01_dp, &
                       6.94134340983508569e-01_dp, &
                       2.78178410488925065e+00_dp]
    dB_dx_ref(5, :) = [-2.77241453503382224e-02_dp, &
                       9.01822519643726728e-01_dp, &
                       3.60561428573554688e+00_dp]

    nabla_s_ref = [7.54027404433375570e-01_dp, &
                   6.72419478593518916e-01_dp, &
                   2.01580646707398969e+00_dp, &
                   8.04748159252873796e-01_dp, &
                   7.80362731036828627e-01_dp]

    ! Compare global quantities
    if (not_same(bfield%psi_tor_edge, psi_tor_edge_ref, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "psi_tor_edge mismatch:"
        print *, "  boozer: ", bfield%psi_tor_edge
        print *, "  neo:    ", psi_tor_edge_ref
        print *, "Relative error: ", abs(bfield%psi_tor_edge - psi_tor_edge_ref) &
            /abs(psi_tor_edge_ref)
        print *, "Absolute error: ", abs(bfield%psi_tor_edge - psi_tor_edge_ref)
        test_failed = .true.
    end if

    if (not_same(bfield%nfp, nfp_ref, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "nfp mismatch:"
        print *, "  boozer: ", bfield%nfp
        print *, "  neo:    ", nfp_ref
        print *, "Relative error: ", abs(bfield%nfp - nfp_ref)/abs(nfp_ref)
        print *, "Absolute error: ", abs(bfield%nfp - nfp_ref)
        test_failed = .true.
    end if

    do case = 1, n_cases
        call bfield%fix_to_surface(stor(case))

        call bfield%get_iota(stor(case), iota_b)
        call bfield%get_covariant_components(B_theta_b, B_phi_b)

        if (not_same(iota_b, iota_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "iota mismatch at case ", case
            print *, "  boozer: ", iota_b
            print *, "  neo:    ", iota_ref(case)
            print *, "Relative error: ", abs(iota_b - iota_ref(case)) &
                /abs(iota_ref(case))
            print *, "Absolute error: ", abs(iota_b - iota_ref(case))
            test_failed = .true.
        end if

        if (not_same(B_theta_b, B_theta_ref(case), &
                     reltol_in=0.0_dp, abstol_in=abstol_for_zero)) then
            print *, "B_theta_covariant mismatch at case ", case
            print *, "  boozer: ", B_theta_b
            print *, "  neo:    ", B_theta_ref(case)
            print *, "Relative error: ", abs(B_theta_b - B_theta_ref(case)) &
                /abs(B_theta_ref(case))
            print *, "Absolute error: ", abs(B_theta_b - B_theta_ref(case))
            test_failed = .true.
        end if

        if (not_same(B_phi_b, B_phi_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_phi_covariant mismatch at case ", case
            print *, "  boozer: ", B_phi_b
            print *, "  neo:    ", B_phi_ref(case)
            print *, "Relative error: ", abs(B_phi_b - B_phi_ref(case)) &
                /abs(B_phi_ref(case))
            print *, "Absolute error: ", abs(B_phi_b - B_phi_ref(case))
            test_failed = .true.
        end if

        call bfield%compute_B_sqrtg_dB_dx(theta(case), phi(case), &
                                          bmod_b, sqrtg_b, dB_dx_b)

        if (not_same(bmod_b, bmod_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_mod mismatch at case ", case
            print *, "  boozer: ", bmod_b
            print *, "  neo:    ", bmod_ref(case)
            print *, "Relative error: ", abs(bmod_b - bmod_ref(case)) &
                /abs(bmod_ref(case))
            print *, "Absolute error: ", abs(bmod_b - bmod_ref(case))
            test_failed = .true.
        end if

        if (not_same(sqrtg_b, sqrtg_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "sqrtg mismatch at case ", case
            print *, "  boozer: ", sqrtg_b
            print *, "  neo:    ", sqrtg_ref(case)
            print *, "Relative error: ", abs(sqrtg_b - sqrtg_ref(case)) &
                /abs(sqrtg_ref(case))
            print *, "Absolute error: ", abs(sqrtg_b - sqrtg_ref(case))
            test_failed = .true.
        end if

        if (not_same(dB_dx_b, dB_dx_ref(case, :), &
                     reltol_in=reltol, abstol_in=abstol_roundoff)) then
            print *, "dB_dx mismatch at case ", case
            print *, "  boozer: ", dB_dx_b
            print *, "  neo:    ", dB_dx_ref(case, :)
            print *, "Relative error: ", abs(dB_dx_b - dB_dx_ref(case, :)) &
                /abs(dB_dx_ref(case, :))
            print *, "Absolute error: ", abs(dB_dx_b - dB_dx_ref(case, :))
            test_failed = .true.
        end if

        call bfield%compute_nabla_s(theta(case), phi(case), nabla_s)

        if (not_same(nabla_s, nabla_s_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "nabla_s mismatch at case ", case
            print *, "  boozer: ", nabla_s
            print *, "  neo:    ", nabla_s_ref(case)
            print *, "Relative error: ", abs(nabla_s - nabla_s_ref(case)) &
                /abs(nabla_s_ref(case))
            print *, "Absolute error: ", abs(nabla_s - nabla_s_ref(case))
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_neo_qh
