module boozer_field

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use vector_potentail_mod, only: torflux

    implicit none
    private

    public :: boozer_field_t

    type :: boozer_field_t
        logical :: initialized = .false.
    contains
        procedure :: init => boozer_field_init
        procedure :: evaluate => boozer_field_evaluate
    end type boozer_field_t

contains

    subroutine boozer_field_init(self, vmec_file, radial_spline_order, &
                                 angular_spline_order, grid_refinement)
        class(boozer_field_t), intent(inout) :: self
        character(len=*), intent(in) :: vmec_file
        integer, intent(in), optional :: radial_spline_order
        integer, intent(in), optional :: angular_spline_order
        integer, intent(in), optional :: grid_refinement

        call get_boozer_coordinates(vmec_file, &
                                    radial_spline_order, &
                                    angular_spline_order, &
                                    grid_refinement)
        self%initialized = .true.
    end subroutine boozer_field_init

    subroutine boozer_field_evaluate(self, x, bmod, sqrtg, bder, hcovar, &
                                     hctrvr, hcurl)
        ! Computes magnetic field quantities in Boozer coordinates.
        !
        ! Input:  x(1)=s (normalized toroidal flux),
        !         x(2)=vartheta_B (Boozer poloidal angle),
        !         x(3)=varphi_B (Boozer toroidal angle).
        !
        ! Output: bmod    - magnetic field module
        !         sqrtg   - Jacobian (sqrt of metric determinant)
        !         bder    - derivatives of log(B)
        !         hcovar  - covariant components of unit vector along B
        !         hctrvr  - contravariant components of unit vector along B
        !         hcurl   - contravariant components of curl of unit vector

        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp) :: r, vartheta_B, varphi_B, &
                    A_phi, A_theta, dA_phi_dr, dA_theta_dr, d2A_phi_dr2, &
                    d3A_phi_dr3, &
                    B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                    B_varphi_B, dB_varphi_B, d2B_varphi_B, Bmod_B, B_r
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

        call splint_boozer_coord(r, vartheta_B, varphi_B, mode_secders, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, &
                                 d2B_vartheta_B, &
                                 B_varphi_B, dB_varphi_B, d2B_varphi_B, &
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

    end subroutine boozer_field_evaluate

end module boozer_field
