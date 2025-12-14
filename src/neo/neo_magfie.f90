module neo_magfie

    use nrtype
    use neo_input, &
        only: es, ixm, ixn, mnmax, psi_pr, pixm, pixn, nfp
    use neo_control, &
        only: fluxs_interp, phi_n, theta_n
    use neo_sub_mod, &
        only: neo_read_control, neo_init
    use neo_spline_data, &
        only: r_mhalf, &
              a_bmnc, b_bmnc, c_bmnc, d_bmnc, &
              a_rmnc, b_rmnc, c_rmnc, d_rmnc, &
              a_zmns, b_zmns, c_zmns, d_zmns, &
              a_vmns, b_vmns, c_vmns, d_vmns, &
              a_iota, b_iota, c_iota, d_iota, &
              a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor, &
              a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol
    use neo_work, &
        only: cosmth, cosnph, sinmth, sinnph, theta_int, phi_int, &
              theta_start, theta_end, phi_start, phi_end
    use neo_spline_mod, only: spl2d, poi2d, eva2d

    implicit none

    integer :: magfie_spline = 1
    integer :: magfie_newspline = 1
    integer :: magfie_sarray_len

    real(dp), dimension(:), allocatable :: curr_tor_array
    real(dp), dimension(:), allocatable :: curr_tor_s_array
    real(dp), dimension(:), allocatable :: curr_pol_array
    real(dp), dimension(:), allocatable :: curr_pol_s_array
    real(dp), dimension(:), allocatable :: iota_array
    real(dp), dimension(:), allocatable :: iota_s_array

    real(dp), dimension(:, :, :, :, :), allocatable, target :: bmod_spl
    real(dp), dimension(:, :, :, :, :), allocatable, target :: bb_s_spl
    real(dp), dimension(:, :, :, :, :), allocatable, target :: bb_tb_spl
    real(dp), dimension(:, :, :, :, :), allocatable, target :: bb_pb_spl
    real(dp), dimension(:, :, :, :, :), allocatable, target :: gval_spl

contains

    !> \brief Calculate magnetic field quantities.
    !>
    !> input:
    !> ------
    !> x: vector of floats, 3 elements?, coordinates?
    !>
    !> output:
    !> -------
    !> bmod: float, magnetic field modulus at given location?
    !> sqrtg: float, square root of metric determinant at given location?
    !> dB_dx: vector of floats, same size as x.
    subroutine neo_magfie_a(x, bmod, sqrtg, dB_dx, iota, &
                            B_theta_covariant, B_phi_covariant, sqrt_g11)

        use neo_exchange, only: b_min, b_max, theta_bmin, theta_bmax, &
            & phi_bmin, phi_bmax
        use neo_work, only: phi_arr, theta_arr
        use neo_spline_over_s, only: spline_1d

        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(out) :: bmod
        real(dp), intent(out) :: sqrtg
        real(dp), dimension(size(x)), intent(out) :: dB_dx
        real(dp), intent(out) :: iota

        real(dp), intent(out), optional :: B_theta_covariant, B_phi_covariant
        real(dp), intent(out), optional :: sqrt_g11

        integer(i4b) :: swd = 1
        integer :: i, m, n
        integer :: npsi
        real(dp) :: m0 = 0.0_dp
        real(dp) :: yp, ypp, yppp

        real(dp) :: sinv, cosv
        real(dp) :: curr_tor
        real(dp) :: curr_pol
        real(dp) :: bb_s, bb_tb, bb_pb

        integer :: k_es = 1
        integer :: imn
        integer :: it, ip, im, in
        integer :: mt = 1
        integer :: mp = 1
        integer :: theta_ind, phi_ind
        integer :: ierr
        integer, dimension(2) :: b_minpos, b_maxpos

        real(dp) :: s
        real(dp) :: bi, bi_s, ri, zi, vi

        real(dp) :: theta_d, phi_d

        real(dp), dimension(:), allocatable :: s_bmnc, s_bmnc_s
        real(dp), dimension(:), allocatable :: s_rmnc, s_zmns, s_vmns

        real(dp), dimension(:, :), allocatable :: bmod_a
        real(dp), dimension(:, :), allocatable :: bb_s_a
        real(dp), dimension(:, :), allocatable :: bb_tb_a
        real(dp), dimension(:, :), allocatable :: bb_pb_a

        real(dp), dimension(:, :), allocatable :: r
        real(dp), dimension(:, :), allocatable :: r_tb, z_tb, v_tb, p_tb
        real(dp), dimension(:, :), allocatable :: r_pb, z_pb, v_pb, p_pb
        real(dp), dimension(:, :), allocatable :: gtbtb, gpbpb, gtbpb
        real(dp), dimension(:, :), allocatable :: sqrg11_met

        real(dp), dimension(:, :, :, :), pointer :: p_spl
        real(dp), dimension(1) :: magfie_sarray

        magfie_sarray(1) = x(1)
        !*******************************************************************
        ! Initialisation if necessary
        !*******************************************************************
        if (.not. allocated(es)) then
            call neo_read_control()
            fluxs_interp = 1
            call neo_init(npsi)
            print *, 'theta_start,theta_end,phi_start,phi_end'
            print *, theta_start, theta_end, phi_start, phi_end
        end if

        !*******************************************************************
        ! Spline of surfaces in magfie_sarray
        !*******************************************************************
        if (magfie_spline .eq. 1 .and. magfie_newspline .eq. 1) then
            magfie_sarray_len = size(magfie_sarray)

            if (allocated(bmod_spl)) deallocate (bmod_spl)
            if (allocated(bb_s_spl)) deallocate (bb_s_spl)
            if (allocated(bb_tb_spl)) deallocate (bb_tb_spl)
            if (allocated(bb_pb_spl)) deallocate (bb_pb_spl)
            if (allocated(gval_spl)) deallocate (gval_spl)
            allocate (bmod_spl(4, 4, theta_n, phi_n, magfie_sarray_len))
            allocate (bb_s_spl(4, 4, theta_n, phi_n, magfie_sarray_len))
            allocate (bb_tb_spl(4, 4, theta_n, phi_n, magfie_sarray_len))
            allocate (bb_pb_spl(4, 4, theta_n, phi_n, magfie_sarray_len))
            allocate (gval_spl(4, 4, theta_n, phi_n, magfie_sarray_len))

            if (allocated(curr_tor_array)) deallocate (curr_tor_array)
            if (allocated(curr_pol_array)) deallocate (curr_pol_array)
            if (allocated(iota_array)) deallocate (iota_array)
            allocate (curr_tor_array(magfie_sarray_len))
            allocate (curr_pol_array(magfie_sarray_len))
            allocate (iota_array(magfie_sarray_len))

            s = magfie_sarray(k_es)

            allocate (s_bmnc(mnmax))
            allocate (s_bmnc_s(mnmax))
            allocate (s_rmnc(mnmax))
            allocate (s_zmns(mnmax))
            allocate (s_vmns(mnmax))

            do imn = 1, mnmax
                ! Switch swd turns on (1) / off (0) the computation of the
                ! radial derivatives within spline_1d
                swd = 1
                call spline_1d(es, &
                & a_bmnc(:, imn), b_bmnc(:, imn), &
                & c_bmnc(:, imn), d_bmnc(:, imn), &
                & swd, r_mhalf(imn),            &
                & s, s_bmnc(imn), s_bmnc_s(imn), ypp, yppp)
                call spline_1d(es, &
                & a_rmnc(:, imn), b_rmnc(:, imn), &
                & c_rmnc(:, imn), d_rmnc(:, imn), &
                & swd, r_mhalf(imn),            &
                & s, s_rmnc(imn), yp, ypp, yppp)
                call spline_1d(es, &
                & a_zmns(:, imn), b_zmns(:, imn), &
                & c_zmns(:, imn), d_zmns(:, imn), &
                & swd, r_mhalf(imn),            &
                & s, s_zmns(imn), yp, ypp, yppp)
                call spline_1d(es, &
                & a_vmns(:, imn), b_vmns(:, imn), &
                & c_vmns(:, imn), d_vmns(:, imn), &
                & swd, r_mhalf(imn),            &
                & s, s_vmns(imn), yp, ypp, yppp)
            end do

            !*************************************************************
            ! Fourier summation for the full theta-phi array
            !*************************************************************
            allocate (bmod_a(theta_n, phi_n))
            allocate (bb_s_a(theta_n, phi_n))
            allocate (bb_tb_a(theta_n, phi_n))
            allocate (bb_pb_a(theta_n, phi_n))
            bmod_a = 0.0_dp
            bb_s_a = 0.0_dp
            bb_tb_a = 0.0_dp
            bb_pb_a = 0.0_dp

            allocate (r(theta_n, phi_n))
            allocate (r_tb(theta_n, phi_n))
            allocate (z_tb(theta_n, phi_n))
            allocate (v_tb(theta_n, phi_n))
            allocate (p_tb(theta_n, phi_n))
            allocate (r_pb(theta_n, phi_n))
            allocate (z_pb(theta_n, phi_n))
            allocate (v_pb(theta_n, phi_n))
            allocate (p_pb(theta_n, phi_n))
            r = 0.0_dp
            r_tb = 0.0_dp
            z_tb = 0.0_dp
            v_tb = 0.0_dp
            p_tb = 0.0_dp
            r_pb = 0.0_dp
            z_pb = 0.0_dp
            v_pb = 0.0_dp
            p_pb = 0.0_dp

            do imn = 1, mnmax
                bi = s_bmnc(imn)
                bi_s = s_bmnc_s(imn)
                ri = s_rmnc(imn)
                zi = s_zmns(imn)
                vi = s_vmns(imn)

                m = ixm(imn)
                n = ixn(imn)
                im = pixm(imn)
                in = pixn(imn)
                do ip = 1, phi_n
                    do it = 1, theta_n

                        cosv = cosmth(it, im)*cosnph(ip, in) &
                               + sinmth(it, im)*sinnph(ip, in)
                        sinv = sinmth(it, im)*cosnph(ip, in) &
                               - cosmth(it, im)*sinnph(ip, in)

                        bmod_a(it, ip) = bmod_a(it, ip) + bi*cosv
                        bb_s_a(it, ip) = bb_s_a(it, ip) + bi_s*cosv
                        bb_tb_a(it, ip) = bb_tb_a(it, ip) - m*bi*sinv
                        bb_pb_a(it, ip) = bb_pb_a(it, ip) + n*bi*sinv

                        r(it, ip) = r(it, ip) + ri*cosv
                        r_tb(it, ip) = r_tb(it, ip) - m*ri*sinv
                        r_pb(it, ip) = r_pb(it, ip) + n*ri*sinv
                        z_tb(it, ip) = z_tb(it, ip) + m*zi*cosv
                        z_pb(it, ip) = z_pb(it, ip) - n*zi*cosv
                        v_tb(it, ip) = v_tb(it, ip) + m*vi*cosv
                        v_pb(it, ip) = v_pb(it, ip) - n*vi*cosv

                    end do
                end do
            end do
            deallocate (s_bmnc)
            deallocate (s_bmnc_s)
            deallocate (s_rmnc)
            deallocate (s_zmns)
            deallocate (s_vmns)

            ! **********************************************************************
            ! Ensure periodicity boundaries to be the same
            ! **********************************************************************
            bmod_a(theta_n, :) = bmod_a(1, :)
            bmod_a(:, phi_n) = bmod_a(:, 1)

            bb_tb_a(theta_n, :) = bb_tb_a(1, :)
            bb_tb_a(:, phi_n) = bb_tb_a(:, 1)

            bb_s_a(theta_n, :) = bb_s_a(1, :)
            bb_s_a(:, phi_n) = bb_s_a(:, 1)

            bb_pb_a(theta_n, :) = bb_pb_a(1, :)
            bb_pb_a(:, phi_n) = bb_pb_a(:, 1)

            r_tb(theta_n, :) = r_tb(1, :)
            r_tb(:, phi_n) = r_tb(:, 1)
            r_pb(theta_n, :) = r_pb(1, :)
            r_pb(:, phi_n) = r_pb(:, 1)
            z_tb(theta_n, :) = z_tb(1, :)
            z_tb(:, phi_n) = z_tb(:, 1)
            z_pb(theta_n, :) = z_pb(1, :)
            z_pb(:, phi_n) = z_pb(:, 1)
            v_tb(theta_n, :) = v_tb(1, :)
            v_tb(:, phi_n) = v_tb(:, 1)
            v_pb(theta_n, :) = v_pb(1, :)
            v_pb(:, phi_n) = v_pb(:, 1)

            ! phi_cyl = phi_boozer - 2pi/nfp * v
            p_tb = -v_tb*twopi/nfp
            p_pb = 1.0_dp - v_pb*twopi/nfp

            gtbtb = r_tb*r_tb + z_tb*z_tb + r*r*p_tb*p_tb
            gpbpb = r_pb*r_pb + z_pb*z_pb + r*r*p_pb*p_pb
            gtbpb = r_tb*r_pb + z_tb*z_pb + r*r*p_tb*p_pb
            sqrg11_met = sqrt(abs(gtbtb*gpbpb - gtbpb*gtbpb))

            p_spl => bmod_spl(:, :, :, :, k_es)
            call spl2d(theta_n, phi_n, theta_int, phi_int, mt, mp, &
                       bmod_a, p_spl)
            p_spl => bb_s_spl(:, :, :, :, k_es)
            call spl2d(theta_n, phi_n, theta_int, phi_int, mt, mp, &
                       bb_s_a, p_spl)
            p_spl => bb_tb_spl(:, :, :, :, k_es)
            call spl2d(theta_n, phi_n, theta_int, phi_int, mt, mp, &
                       bb_tb_a, p_spl)
            p_spl => bb_pb_spl(:, :, :, :, k_es)
            call spl2d(theta_n, phi_n, theta_int, phi_int, mt, mp, &
                       bb_pb_a, p_spl)
            p_spl => gval_spl(:, :, :, :, k_es)
            call spl2d(theta_n, phi_n, theta_int, phi_int, mt, mp, &
                       sqrg11_met, p_spl)

            deallocate (bb_s_a)
            deallocate (bb_tb_a)
            deallocate (bb_pb_a)
            deallocate (sqrg11_met)

            swd = 0
            call spline_1d(es, a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor, &
                           swd, m0, s, curr_tor_array(k_es), &
                           yp, ypp, yppp)
            call spline_1d(es, a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol, &
                           swd, m0, s, curr_pol_array(k_es), &
                           yp, ypp, yppp)
            call spline_1d(es, a_iota, b_iota, c_iota, d_iota, &
                           swd, m0, s, iota_array(k_es), yp, ypp, yppp)
            magfie_newspline = 0

            ! Minimum and Maximum in the new mode
            if (magfie_sarray_len .eq. 1) then
                ! **********************************************************************
                ! Calculate absolute minimum and maximum of b and its location (theta, phi)
                ! **********************************************************************
                b_minpos = minloc(bmod_a)
                b_min = bmod_a(b_minpos(1), b_minpos(2))
                theta_bmin = theta_arr(b_minpos(1))
                phi_bmin = phi_arr(b_minpos(2))

                b_maxpos = maxloc(bmod_a)
                b_max = bmod_a(b_maxpos(1), b_maxpos(2))
                theta_bmax = theta_arr(b_maxpos(1))
                phi_bmax = phi_arr(b_maxpos(2))
            end if
            deallocate (bmod_a)

        end if

        curr_tor = curr_tor_array(k_es)
        curr_pol = curr_pol_array(k_es)
        iota = iota_array(k_es)

        call poi2d(theta_int, phi_int, mt, mp, &
                   theta_start, theta_end, phi_start, phi_end, &
                   x(3), x(2), theta_ind, phi_ind, theta_d, phi_d, ierr)
        p_spl => bmod_spl(:, :, :, :, k_es)
        call eva2d(theta_n, phi_n, theta_ind, phi_ind, theta_d, phi_d, &
                   p_spl, bmod)
        p_spl => bb_s_spl(:, :, :, :, k_es)
        call eva2d(theta_n, phi_n, theta_ind, phi_ind, theta_d, phi_d, &
                   p_spl, bb_s)
        p_spl => bb_tb_spl(:, :, :, :, k_es)
        call eva2d(theta_n, phi_n, theta_ind, phi_ind, theta_d, phi_d, &
                   p_spl, bb_tb)
        p_spl => bb_pb_spl(:, :, :, :, k_es)
        call eva2d(theta_n, phi_n, theta_ind, phi_ind, theta_d, phi_d, &
                   p_spl, bb_pb)
        if (present(sqrt_g11)) then
            p_spl => gval_spl(:, :, :, :, k_es)
            call eva2d(theta_n, phi_n, theta_ind, phi_ind, theta_d, phi_d, &
                       p_spl, sqrt_g11)
            sqrt_g11 = sqrt_g11*bmod**2.0_dp/abs(curr_pol + iota*curr_tor)
        end if
        sqrtg = psi_pr*(curr_pol + iota*curr_tor)/bmod**2.0_dp*1d6

        dB_dx(1) = bb_s
        dB_dx(3) = bb_tb
        dB_dx(2) = bb_pb

        if (present(B_theta_covariant)) B_theta_covariant = -curr_tor
        if (present(B_phi_covariant)) B_phi_covariant = -curr_pol

    end subroutine neo_magfie_a

end module neo_magfie
