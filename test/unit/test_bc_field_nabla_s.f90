!> Validate bc_field_t against an analytically consistent single-surface
!> .bc file: a circular-cross-section torus with concentric surfaces
!> r(s) = a*sqrt(s) and a purely toroidal field B = B0*sqrt(R0/R).
!>
!> Choosing flux = pi*a^2*B0 makes the Boozer Jacobian identity
!> sqrt(g) = psi_pr*(B_phi + iota*B_theta)/B^2 hold exactly for the
!> geometric Jacobian r*R*dr/ds, so |nabla s| = 2*sqrt(s)/a exactly.
program test_bc_field_nabla_s
    use constants, only: dp, pi
    use bc_field, only: bc_field_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: R0 = 10.0_dp, a = 1.0_dp, B0 = 5.0_dp
    real(dp), parameter :: s0 = 0.49_dp
    real(dp), parameter :: r = a*sqrt(s0)
    real(dp), parameter :: eps = r/R0
    real(dp), parameter :: mu0_over_two_pi = 2.0e-7_dp
    integer, parameter :: m_max = 24

    type(bc_field_t) :: field, field_neo2, field_vflip
    real(dp) :: bmnc(0:m_max)
    real(dp) :: theta, phi, B_mod, B_expected, nabla_s, nabla_s_expected
    real(dp) :: nabla_s_neo2, nabla_s_vflip
    integer :: i_theta, i_phi
    logical :: test_failed

    test_failed = .false.

    call compute_bmnc_harmonics(bmnc)
    call write_test_bc("nabla_s_test.bc", bmnc, 0.0_dp)
    call write_test_bc("nabla_s_test_vflip.bc", bmnc, -0.05_dp)

    call field%bc_field_init("nabla_s_test.bc")
    call field_neo2%bc_field_init("nabla_s_test.bc", neo2_nabla_s=.true.)
    call field_vflip%bc_field_init("nabla_s_test_vflip.bc")

    if (not_same(field%psi_tor_edge, 0.5_dp*a**2*B0, reltol_in=1e-8_dp)) then
        print *, "psi_tor_edge mismatch: ", field%psi_tor_edge
        test_failed = .true.
    end if
    if (not_same(field%R, R0, reltol_in=1e-8_dp)) then
        print *, "R mismatch: ", field%R
        test_failed = .true.
    end if

    call field%fix_to_surface(s0)
    nabla_s_expected = 2.0_dp*sqrt(s0)/a

    do i_theta = 0, 7
        do i_phi = 0, 3
            theta = real(i_theta, dp)*2.0_dp*pi/8.0_dp
            phi = real(i_phi, dp)*2.0_dp*pi/4.0_dp

            call field%compute_B_mod(theta, phi, B_mod)
            B_expected = B0/sqrt(1.0_dp + eps*cos(theta))
            if (not_same(B_mod, B_expected, reltol_in=1e-12_dp)) then
                print *, "B_mod mismatch at theta =", theta, ":", &
                    B_mod, B_expected
                test_failed = .true.
            end if

            call field%compute_nabla_s(theta, phi, nabla_s)
            if (not_same(nabla_s, nabla_s_expected, reltol_in=1e-10_dp)) then
                print *, "nabla_s mismatch at theta =", theta, ":", &
                    nabla_s, nabla_s_expected
                test_failed = .true.
            end if
        end do
    end do
    print *, "nabla_s on circular torus: ", nabla_s, &
        " expected: ", nabla_s_expected

    ! The v-sign switch: reading vmns with neo2_nabla_s=.true. must equal
    ! reading -vmns with the default convention.
    call field_neo2%fix_to_surface(s0)
    call field_vflip%fix_to_surface(s0)
    theta = 0.7_dp
    phi = 0.3_dp
    call field%compute_nabla_s(theta, phi, nabla_s)
    call field_neo2%compute_nabla_s(theta, phi, nabla_s_neo2)
    if (not_same(nabla_s, nabla_s_neo2, reltol_in=1e-14_dp)) then
        print *, "neo2_nabla_s must not change the result for vmns = 0:", &
            nabla_s, nabla_s_neo2
        test_failed = .true.
    end if

    call field_vflip%compute_nabla_s(theta, phi, nabla_s_vflip)
    call read_neo2_of_vflipped_file(bmnc, nabla_s_neo2, theta, phi)
    if (not_same(nabla_s_vflip, nabla_s_neo2, reltol_in=1e-14_dp)) then
        print *, "v-sign switch is not equivalent to flipping vmns:", &
            nabla_s_vflip, nabla_s_neo2
        test_failed = .true.
    end if

    if (test_failed) error stop "test_bc_field_nabla_s failed"
    print *, "test_bc_field_nabla_s passed"

contains

    !> nabla_s of the +vmns file read with the legacy NEO-2 convention.
    subroutine read_neo2_of_vflipped_file(bmnc, nabla_s, theta, phi)
        real(dp), intent(in) :: bmnc(0:), theta, phi
        real(dp), intent(out) :: nabla_s

        type(bc_field_t) :: field_tmp

        call write_test_bc("nabla_s_test_vplus.bc", bmnc, 0.05_dp)
        call field_tmp%bc_field_init("nabla_s_test_vplus.bc", &
                                     neo2_nabla_s=.true.)
        call field_tmp%fix_to_surface(s0)
        call field_tmp%compute_nabla_s(theta, phi, nabla_s)
    end subroutine read_neo2_of_vflipped_file

    !> Fourier cosine harmonics of B0/sqrt(1 + eps*cos(theta)) by trapezoid
    !> rule (spectrally accurate for periodic integrands).
    subroutine compute_bmnc_harmonics(bmnc)
        real(dp), intent(out) :: bmnc(0:m_max)

        integer, parameter :: n_grid = 4096
        integer :: i, m
        real(dp) :: theta_i, f_i

        bmnc = 0.0_dp
        do i = 0, n_grid - 1
            theta_i = real(i, dp)*2.0_dp*pi/real(n_grid, dp)
            f_i = B0/sqrt(1.0_dp + eps*cos(theta_i))
            do m = 0, m_max
                bmnc(m) = bmnc(m) + f_i*cos(real(m, dp)*theta_i)
            end do
        end do
        bmnc = bmnc*2.0_dp/real(n_grid, dp)
        bmnc(0) = 0.5_dp*bmnc(0)
    end subroutine compute_bmnc_harmonics

    subroutine write_test_bc(filename, bmnc, vmns_1)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: bmnc(0:m_max)
        real(dp), intent(in) :: vmns_1

        integer :: iunit, m
        real(dp) :: flux, Jpol_over_nper, rmnc_m, zmns_m, vmns_m

        flux = pi*a**2*B0
        Jpol_over_nper = -B0*R0/mu0_over_two_pi

        open (newunit=iunit, file=filename, status='replace', action='write')
        write (iunit, '(A)') 'CC analytic circular-torus test file'
        write (iunit, '(A)') ' m0b   n0b  nsurf  nper    flux [Tm^2]'// &
            '        a [m]          R [m]'
        write (iunit, '(4I6, 3ES23.15)') m_max, 0, 1, 1, flux, a, R0
        write (iunit, '(A)') '        s               iota'// &
            '           Jpol/nper          Itor            pprime'// &
            '         sqrt g(0,0)'
        write (iunit, '(A)') '                                          [A]'// &
            '           [A]             [Pa]         (dV/ds)/nper'
        write (iunit, '(6ES23.15)') s0, 0.0_dp, Jpol_over_nper, 0.0_dp, &
            0.0_dp, 0.0_dp
        write (iunit, '(A)') '    m    n      rmnc [m]'// &
            '         zmns [m]         vmns [ ]         bmnc [T]'
        do m = 0, m_max
            rmnc_m = 0.0_dp
            zmns_m = 0.0_dp
            vmns_m = 0.0_dp
            if (m == 0) rmnc_m = R0
            if (m == 1) then
                rmnc_m = r
                zmns_m = r
                vmns_m = vmns_1
            end if
            write (iunit, '(2I6, 4ES23.15)') m, 0, rmnc_m, zmns_m, vmns_m, &
                bmnc(m)
        end do
        close (iunit)
    end subroutine write_test_bc

end program test_bc_field_nabla_s
