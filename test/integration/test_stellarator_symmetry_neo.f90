!> Test that neo_field_t respects stellarator symmetry and periodicity:
!>   B(theta,phi)             = B(-theta,-phi)   (stellarator symmetry)
!>   B(theta + 2pi, phi)      = B(theta,phi)      (poloidal periodicity)
!>   B(theta, phi + 2pi/nfp)  = B(theta,phi)      (toroidal periodicity)
program test_stellarator_symmetry_neo
    use constants, only: dp, pi
    use neo_field, only: neo_field_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol_symmetry = 1.e-15_dp, abstol = 0.0_dp
    real(dp), parameter :: reltol_periodicity = 10.0_dp*epsilon(1.0_dp)
    character(len=*), parameter :: bc_filename = "input/single_mode_m_2_n_minus4.bc"

    type(neo_field_t) :: nfield

    integer, parameter :: n_stor = 3
    integer, parameter :: n_pts = 5
    real(dp), parameter :: stors(n_stor) = [0.2_dp, 0.5_dp, 0.8_dp]
    real(dp), parameter :: theta(n_pts) = [0.3_dp, 1.0_dp, 2.0_dp, -0.5_dp, 3.0_dp]
    real(dp), parameter :: phi(n_pts) = [0.1_dp, 0.7_dp, 1.5_dp, 0.4_dp, 2.8_dp]

    real(dp) :: bmod_ref
    real(dp) :: bmod_sym
    real(dp) :: max_sym_err, phi_period
    integer :: is, ip
    logical :: test_failed

    call nfield%neo_field_init(bc_filename, stors(1))

    test_failed = .false.
    max_sym_err = 0.0_dp
    phi_period = 2.0_dp*pi/nfield%nfp

    do is = 1, n_stor
        print *, "Testing stellarator symmetry and periodicity at s_tor =", stors(is)
        call nfield%neo_change_stor(stors(is))

        do ip = 1, n_pts
            call nfield%compute_B_mod( &
                theta(ip), phi(ip), bmod_ref)

            ! --- stellarator symmetry ---
            call nfield%compute_B_mod( &
                -theta(ip), -phi(ip), bmod_sym)

            max_sym_err = max(max_sym_err, abs(bmod_ref - bmod_sym)/abs(bmod_ref))

            if (not_same(bmod_ref, bmod_sym, &
                         reltol_in=reltol_symmetry, abstol_in=abstol)) then
                print *, "Stellarator symmetry violated: B(theta,phi) /= B(-theta,-phi)"
                print *, "  s_tor =", stors(is), &
                    "  theta =", theta(ip), "  phi =", phi(ip)
                print *, "  B(theta,phi)   =", bmod_ref
                print *, "  B(-theta,-phi) =", bmod_sym
                print *, "  Relative error =", abs(bmod_ref - bmod_sym)/abs(bmod_ref)
                test_failed = .true.
            end if

            ! --- poloidal periodicity ---
            call nfield%compute_B_mod( &
                theta(ip) + 2.0_dp*pi, phi(ip), bmod_sym)

            if (not_same(bmod_ref, bmod_sym, &
                         reltol_in=reltol_periodicity, abstol_in=abstol)) then
                print *, "Poloidal periodicity violated:"
                print *, "  s_tor =", stors(is), &
                    "  theta =", theta(ip), "  phi =", phi(ip)
                print *, "  B(theta,phi)       =", bmod_ref
                print *, "  B(theta+2pi,phi)   =", bmod_sym
                print *, "  Relative error =", abs(bmod_ref - bmod_sym)/abs(bmod_ref)
                test_failed = .true.
            end if

            ! --- toroidal periodicity ---
            call nfield%compute_B_mod( &
                theta(ip), phi(ip) + phi_period, bmod_sym)

  if (not_same(bmod_ref, bmod_sym, reltol_in=reltol_periodicity, abstol_in=abstol)) then
                print *, "Toroidal periodicity violated:"
                print *, "  s_tor =", stors(is), &
                    "  theta =", theta(ip), "  phi =", phi(ip), &
                    "  nfp =", nfield%nfp
                print *, "  B(theta,phi)           =", bmod_ref
                print *, "  B(theta,phi+2pi/nfp)   =", bmod_sym
                print *, "  Relative error =", abs(bmod_ref - bmod_sym)/abs(bmod_ref)
                test_failed = .true.
            end if
        end do
    end do

    print *, "Max stellarator symmetry relative error:", max_sym_err

    if (test_failed) error stop

end program test_stellarator_symmetry_neo
