module magfie_can_boozer_sub

    use libneo_kinds, only: dp
    use boozer_sub, only: splint_boozer_coord
    use vector_potentail_mod, only: torflux

    implicit none

contains

    subroutine magfie_boozer(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        ! Computes magnetic field module in units of the magnetic code  - bmod,
        ! square root of determinant of the metric tensor               - sqrtg,
        ! derivatives of the logarythm of the magnetic field module
        ! over coordinates                                              - bder,
        ! covariant componets of the unit vector of the magnetic
        ! field direction                                               - hcovar,
        ! contravariant components of this vector                       - hctrvr,
        ! contravariant component of the curl of this vector            - hcurl
        ! Order of coordinates is the following: x(1)=s (normalized toroidal flux),
        ! x(2)=vartheta_B (Boozer's poloidal angle), x(3)=varphi_B (Boozer's toroidal angle).
        !
        !  Input parameters:
        !            formal:  x(3)             - array of Boozer coordinates
        !  Output parameters:
        !            formal:  bmod
        !                     sqrtg
        !                     bder(3)          - derivatives of $\log(B)$
!                     hcovar(3)        - covariant components of unit vector $\bh$
!                                        along $\bB$
!                     hctrvr(3)        - contra-variant components of unit vector $\bh$
!                                        along $\bB$
        !                     hcurl(3)         - contra-variant components of curl of $\bh$
        !
        !  Called routines: splint_boozer_coord

        implicit none

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

        r = x(1)
        vartheta_B = x(2)
        varphi_B = x(3)

        call splint_boozer_coord(r, vartheta_B, varphi_B, mode_secders, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, &
                                 d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
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
        hcurl(1) = (B_vartheta_B*bder(3) - B_varphi_B*bder(2))/sqrtgbmod
        hcurl(2) = (B_varphi_B*bder(1) - B_r*bder(3) + dB_r(3) - dB_varphi_B)/sqrtgbmod
        hcurl(3) = (B_r*bder(2) - B_vartheta_B*bder(1) + dB_vartheta_B - &
                    dB_r(2))/sqrtgbmod

    end subroutine magfie_boozer

end module magfie_can_boozer_sub
