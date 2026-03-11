module boozer_field

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: field_t
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use vector_potentail_mod, only: torflux

    implicit none
    private

    public :: boozer_field_t

    type, extends(field_t) :: boozer_field_t
        logical :: initialized = .false.
        real(dp) :: s
        real(dp) :: iota
        real(dp) :: nfp
        real(dp) :: psi_tor_edge
        real(dp) :: B_theta_covariant, B_phi_covariant
    contains
        procedure :: boozer_field_init
        procedure :: set_stor
        procedure :: evaluate
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_sqrt_g11
    end type boozer_field_t

contains

    subroutine boozer_field_init(self, vmec_file, s, &
                                 radial_spline_order, &
                                 angular_spline_order, grid_refinement)
        class(boozer_field_t), intent(inout) :: self
        character(len=*), intent(in) :: vmec_file
        real(dp), intent(in) :: s
        integer, intent(in), optional :: radial_spline_order
        integer, intent(in), optional :: angular_spline_order
        integer, intent(in), optional :: grid_refinement

        call get_boozer_coordinates(vmec_file, &
                                    radial_spline_order, &
                                    angular_spline_order, &
                                    grid_refinement)
        self%initialized = .true.
        call self%set_stor(s)
    end subroutine boozer_field_init

    subroutine set_stor(self, s)
        use new_vmec_stuff_mod, only: nper

        class(boozer_field_t), intent(inout) :: self
        real(dp), intent(in) :: s

        real(dp) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr
        real(dp) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: B_vartheta_B, dB_vartheta_B, d2B_vartheta_B
        real(dp) :: B_varphi_B, dB_varphi_B, d2B_varphi_B
        real(dp) :: Bmod_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        self%s = s

        call splint_boozer_coord(s, 0.0_dp, 0.0_dp, 0, &
                                 A_theta, A_phi, dA_theta_dr, &
                                 dA_phi_dr, d2A_phi_dr2, &
                                 d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, &
                                 d2B_vartheta_B, &
                                 B_varphi_B, dB_varphi_B, &
                                 d2B_varphi_B, &
                                 Bmod_B, dBmod_B, d2Bmod_B, &
                                 B_r, dB_r, d2B_r)

        self%iota = -dA_phi_dr/dA_theta_dr
        self%B_theta_covariant = B_vartheta_B
        self%B_phi_covariant = B_varphi_B
        self%psi_tor_edge = torflux
        self%nfp = real(nper, dp)
    end subroutine set_stor

    subroutine evaluate(self, x, bmod, sqrtg, bder, &
                        hcovar, hctrvr, hcurl)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp) :: r, vartheta_B, varphi_B, &
                    A_phi, A_theta, dA_phi_dr, dA_theta_dr, &
                    d2A_phi_dr2, d3A_phi_dr3, &
                    B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                    B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                    Bmod_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        real(dp) :: aiota, Bctrvr_theta, Bctrvr_phi, sqrtgbmod

        integer, parameter :: mode_secders = 0

        if (.not. self%initialized) then
            error stop "boozer_field_evaluate: field not initialized"
        end if

        r = x(1)
        vartheta_B = x(2)
        varphi_B = x(3)

        call splint_boozer_coord(r, vartheta_B, varphi_B, &
                                 mode_secders, &
                                 A_theta, A_phi, dA_theta_dr, &
                                 dA_phi_dr, d2A_phi_dr2, &
                                 d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, &
                                 d2B_vartheta_B, &
                                 B_varphi_B, dB_varphi_B, &
                                 d2B_varphi_B, &
                                 Bmod_B, dBmod_B, d2Bmod_B, &
                                 B_r, dB_r, d2B_r)

        aiota = -dA_phi_dr/dA_theta_dr

        bmod = Bmod_B
        bder = dBmod_B/Bmod_B

        sqrtg = (aiota*B_vartheta_B + B_varphi_B)/bmod**2*torflux

        Bctrvr_phi = dA_theta_dr/sqrtg
        Bctrvr_theta = aiota*Bctrvr_phi
        hctrvr(1) = 0.d0
        hctrvr(2) = Bctrvr_theta/bmod
        hctrvr(3) = Bctrvr_phi/bmod

        hcovar(1) = B_r/bmod
        hcovar(2) = B_vartheta_B/bmod
        hcovar(3) = B_varphi_B/bmod

        sqrtgbmod = sqrtg*bmod
        hcurl(1) = (B_vartheta_B*bder(3) &
                    - B_varphi_B*bder(2))/sqrtgbmod
        hcurl(2) = (B_varphi_B*bder(1) - B_r*bder(3) &
                    + dB_r(3) - dB_varphi_B)/sqrtgbmod
        hcurl(3) = (B_r*bder(2) - B_vartheta_B*bder(1) &
                    + dB_vartheta_B - dB_r(2))/sqrtgbmod

    end subroutine evaluate

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, &
                                     dB_dx)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        real(dp) :: x(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [self%s, theta, phi]
        call self%evaluate(x, B_mod, sqrtg, dB_dx, hcovar, &
                           hctrvr, hcurl)
    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: x(3), dummy_sqrtg, hcovar(3), hctrvr(3), hcurl(3)

        x = [self%s, theta, phi]
        call self%evaluate(x, B_mod, dummy_sqrtg, dB_dx, hcovar, &
                           hctrvr, hcurl)
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: x(3), dummy_sqrtg, dummy_dB_dx(3)
        real(dp) :: hcovar(3), hctrvr(3), hcurl(3)

        x = [self%s, theta, phi]
        call self%evaluate(x, B_mod, dummy_sqrtg, dummy_dB_dx, &
                           hcovar, hctrvr, hcurl)
    end subroutine compute_B_mod

    subroutine compute_sqrt_g11(self, theta, phi, sqrt_g11)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: sqrt_g11

        sqrt_g11 = 0.0_dp
        error stop "compute_sqrt_g11 not implemented for &
            &boozer_field_t"
    end subroutine compute_sqrt_g11

end module boozer_field
