module boozer_sub
    use spl_three_to_five_sub
    use interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                           construct_batch_splines_1d, construct_batch_splines_3d, &
                           evaluate_batch_splines_1d_der2, &
                           evaluate_batch_splines_1d_der3, &
                           evaluate_batch_splines_3d_der, &
                           evaluate_batch_splines_3d_der2, &
                           destroy_batch_splines_1d, destroy_batch_splines_3d
    use math_constants, only: TWOPI
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    ! Public API
    public :: get_boozer_coordinates
    public :: get_boozer_coordinates_from_chartmap
    public :: get_boozer_coordinates_from_boozmn
    public :: splint_boozer_coord
    public :: reset_boozer_batch_splines

    integer, parameter :: MAX_FIELD3D_QUANTITIES = 3

    ! Batch spline data for 3D field quantities (Bmod, sqrt_g_ss, optionally B_r)
    type(BatchSplineData3D), save :: field3d_batch_spline
    logical, save :: field3d_batch_spline_ready = .false.
    integer, save :: field3d_num_quantities = 0
    real(dp), allocatable, save :: bmod_grid(:, :, :)
    real(dp), allocatable, save :: br_grid(:, :, :)
    real(dp), allocatable, save :: sqrt_g_ss_grid(:, :, :)

    ! Batch spline for A_phi (vector potential)
    type(BatchSplineData1D), save :: aphi_batch_spline
    logical, save :: aphi_batch_spline_ready = .false.

    ! Batch spline for B_theta, B_phi covariant components
    type(BatchSplineData1D), save :: bcovar_tp_batch_spline
    logical, save :: bcovar_tp_batch_spline_ready = .false.

contains

    !> Initialize Boozer coordinates using VMEC field (backward compatibility)
    subroutine get_boozer_coordinates(vmec_file, &
                                      radial_spline_order, &
                                      angular_spline_order, &
                                      grid_refinment)
        use new_vmec_stuff_mod, only: netcdffile, ns_s, ns_tp, multharm
        use spline_vmec_sub, only: spline_vmec_data

        character(len=*), intent(in) :: vmec_file
        integer, intent(in), optional :: radial_spline_order, angular_spline_order, grid_refinment

        netcdffile = vmec_file
        if (present(radial_spline_order)) then
            ns_s = radial_spline_order
        else
            ns_s = 5
        end if
        if (present(angular_spline_order)) then
            ns_tp = angular_spline_order
        else
            ns_tp = 5
        end if
        if (present(grid_refinment)) then
            multharm = grid_refinment
        else
            multharm = 3
        end if
        call spline_vmec_data()
        call reset_boozer_batch_splines()
        call get_boozer_coordinates_impl()

    end subroutine get_boozer_coordinates

    subroutine get_boozer_coordinates_impl

        use vector_potentail_mod, only: ns, hs
        use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, ns_s, ns_tp
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, use_B_r

        implicit none

        ns_s_B = ns_s
        ns_tp_B = ns_tp
        ns_B = ns
        n_theta_B = n_theta
        n_phi_B = n_phi

        hs_B = hs*real(ns - 1, dp)/real(ns_B - 1, dp)
        h_theta_B = h_theta*real(n_theta - 1, dp)/real(n_theta_B - 1, dp)
        h_phi_B = h_phi*real(n_phi - 1, dp)/real(n_phi_B - 1, dp)

        call compute_boozer_data

        call build_boozer_aphi_batch_spline
        call build_boozer_bcovar_tp_batch_spline
        call build_boozer_field3d_batch_spline

    end subroutine get_boozer_coordinates_impl

    subroutine splint_boozer_coord(r, vartheta_B, varphi_B, mode_secders, &
                                   A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                   d2A_phi_dr2, d3A_phi_dr3, &
                                   B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                                   B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                                   Bmod_B, dBmod_B, d2Bmod_B, &
                                   sqrt_g_ss_B, &
                                   B_r, dB_r, d2B_r)

        use boozer_coordinates_mod, only: use_B_r
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper
        use chamb_mod, only: rnegflag
        use diag_mod, only: dodiag, icounter

        implicit none

        integer, intent(in) :: mode_secders

        real(dp), intent(in) :: r, vartheta_B, varphi_B
        real(dp), intent(out) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr
        real(dp), intent(out) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp), intent(out) :: B_vartheta_B, dB_vartheta_B, d2B_vartheta_B
        real(dp), intent(out) :: B_varphi_B, dB_varphi_B, d2B_varphi_B
        real(dp), intent(out) :: Bmod_B, sqrt_g_ss_B, B_r
        real(dp), intent(out) :: dBmod_B(3), dB_r(3)
        real(dp), intent(out) :: d2Bmod_B(6), d2B_r(6)

        integer :: i_br
        real(dp) :: r_eval, rho_tor, drhods, drhods2, d2rhods2m
        real(dp) :: qua, dqua_dr, dqua_dt, dqua_dp
        real(dp) :: d2qua_dr2, d2qua_drdt, d2qua_drdp, d2qua_dt2, &
                    d2qua_dtdp, d2qua_dp2
        real(dp) :: x_eval(3)
        real(dp) :: y_eval(MAX_FIELD3D_QUANTITIES)
        real(dp) :: dy_eval(3, MAX_FIELD3D_QUANTITIES)
        real(dp) :: d2y_eval(6, MAX_FIELD3D_QUANTITIES)
        real(dp) :: theta_wrapped, phi_wrapped
        real(dp) :: y1d(2), dy1d(2), d2y1d(2)

        if (dodiag) then
!$omp atomic
            icounter = icounter + 1
        end if
        r_eval = r
        if (r_eval .le. 0.0_dp) then
            rnegflag = .true.
            r_eval = abs(r_eval)
        end if

        A_theta = torflux*r_eval
        dA_theta_dr = torflux

        ! Interpolate A_phi over s (batch spline 1D)
        if (.not. aphi_batch_spline_ready) then
            error stop "splint_boozer_coord: Aphi batch spline not initialized"
        end if

        if (mode_secders > 0) then
            ! Need third derivative - use der3 which computes all in one pass
            block
                real(dp) :: d3y1d(1)
                call evaluate_batch_splines_1d_der3(aphi_batch_spline, r_eval, &
                                                    y1d(1:1), dy1d(1:1), &
                                                    d2y1d(1:1), d3y1d)
                d3A_phi_dr3 = d3y1d(1)
            end block
        else
            call evaluate_batch_splines_1d_der2(aphi_batch_spline, r_eval, y1d(1:1), &
                                                dy1d(1:1), d2y1d(1:1))
            d3A_phi_dr3 = 0.0_dp
        end if
        A_phi = y1d(1)
        dA_phi_dr = dy1d(1)
        d2A_phi_dr2 = d2y1d(1)

        ! Interpolation of mod-B (and B_r if use_B_r)
        rho_tor = sqrt(r_eval)
        theta_wrapped = modulo(vartheta_B, TWOPI)
        phi_wrapped = modulo(varphi_B, TWOPI/real(nper, dp))

        if (.not. field3d_batch_spline_ready) then
            error stop "splint_boozer_coord: Bmod/Br batch spline not initialized"
        end if

        x_eval(1) = rho_tor
        x_eval(2) = theta_wrapped
        x_eval(3) = phi_wrapped

        i_br = field3d_num_quantities

        ! Chain rule coefficients for rho -> s conversion
        drhods = 0.5_dp/rho_tor
        drhods2 = drhods**2
        d2rhods2m = drhods2/rho_tor  ! -d2rho/ds2 (negative of second derivative)

        if (mode_secders == 2) then
            call evaluate_batch_splines_3d_der2(field3d_batch_spline, x_eval, &
                                                y_eval(1:field3d_num_quantities), &
                                                dy_eval(:, 1:field3d_num_quantities), &
                                                d2y_eval(:, 1:field3d_num_quantities))

            ! Extract Bmod (quantity 1)
            qua = y_eval(1)
            dqua_dr = dy_eval(1, 1)
            dqua_dt = dy_eval(2, 1)
            dqua_dp = dy_eval(3, 1)

            d2qua_dr2 = d2y_eval(1, 1)
            d2qua_drdt = d2y_eval(2, 1)
            d2qua_drdp = d2y_eval(3, 1)
            d2qua_dt2 = d2y_eval(4, 1)
            d2qua_dtdp = d2y_eval(5, 1)
            d2qua_dp2 = d2y_eval(6, 1)

            d2qua_dr2 = d2qua_dr2*drhods2 - dqua_dr*d2rhods2m
            dqua_dr = dqua_dr*drhods
            d2qua_drdt = d2qua_drdt*drhods
            d2qua_drdp = d2qua_drdp*drhods

            Bmod_B = qua

            dBmod_B(1) = dqua_dr
            dBmod_B(2) = dqua_dt
            dBmod_B(3) = dqua_dp

            d2Bmod_B(1) = d2qua_dr2
            d2Bmod_B(2) = d2qua_drdt
            d2Bmod_B(3) = d2qua_drdp
            d2Bmod_B(4) = d2qua_dt2
            d2Bmod_B(5) = d2qua_dtdp
            d2Bmod_B(6) = d2qua_dp2

            sqrt_g_ss_B = y_eval(2)

            ! Extract B_r (if present)
            if (use_B_r) then
                qua = y_eval(i_br)
                dqua_dr = dy_eval(1, i_br)
                dqua_dt = dy_eval(2, i_br)
                dqua_dp = dy_eval(3, i_br)

                d2qua_dr2 = d2y_eval(1, i_br)
                d2qua_drdt = d2y_eval(2, i_br)
                d2qua_drdp = d2y_eval(3, i_br)
                d2qua_dt2 = d2y_eval(4, i_br)
                d2qua_dtdp = d2y_eval(5, i_br)
                d2qua_dp2 = d2y_eval(6, i_br)

                d2qua_dr2 = d2qua_dr2*drhods2 - dqua_dr*d2rhods2m
                dqua_dr = dqua_dr*drhods
                d2qua_drdt = d2qua_drdt*drhods
                d2qua_drdp = d2qua_drdp*drhods

                B_r = qua*drhods

                dB_r(1) = dqua_dr*drhods - qua*d2rhods2m
                dB_r(2) = dqua_dt*drhods
                dB_r(3) = dqua_dp*drhods

                d2B_r(1) = d2qua_dr2*drhods - 2.0_dp*dqua_dr*d2rhods2m + &
                           qua*drhods*(3.0_dp/4.0_dp)/r_eval**2
                d2B_r(2) = d2qua_drdt*drhods - dqua_dt*d2rhods2m
                d2B_r(3) = d2qua_drdp*drhods - dqua_dp*d2rhods2m
                d2B_r(4) = d2qua_dt2*drhods
                d2B_r(5) = d2qua_dtdp*drhods
                d2B_r(6) = d2qua_dp2*drhods
            else
                B_r = 0.0_dp
                dB_r = 0.0_dp
                d2B_r = 0.0_dp
            end if
        else
            call evaluate_batch_splines_3d_der(field3d_batch_spline, x_eval, &
                                               y_eval(1:field3d_num_quantities), &
                                               dy_eval(:, 1:field3d_num_quantities))

            Bmod_B = y_eval(1)
            dBmod_B(1) = dy_eval(1, 1)*drhods
            dBmod_B(2) = dy_eval(2, 1)
            dBmod_B(3) = dy_eval(3, 1)

            sqrt_g_ss_B = y_eval(2)

            d2Bmod_B = 0.0_dp

            if (mode_secders == 1) then
                call evaluate_batch_splines_3d_der2(field3d_batch_spline, x_eval, &
                                                    y_eval(1:field3d_num_quantities), &
                                                    dy_eval(:, &
                                                            1:field3d_num_quantities), &
                                                    d2y_eval(:, &
                                                             1:field3d_num_quantities))
                d2Bmod_B(1) = d2y_eval(1, 1)*drhods2 - dy_eval(1, 1)*d2rhods2m
            end if

            if (use_B_r) then
                qua = y_eval(i_br)
                dqua_dr = dy_eval(1, i_br)
                dqua_dt = dy_eval(2, i_br)
                dqua_dp = dy_eval(3, i_br)

                dqua_dr = dqua_dr*drhods
                B_r = qua*drhods

                dB_r(1) = dqua_dr*drhods - qua*d2rhods2m
                dB_r(2) = dqua_dt*drhods
                dB_r(3) = dqua_dp*drhods

                d2B_r = 0.0_dp
                if (mode_secders == 1) then
                    d2qua_dr2 = d2y_eval(1, i_br)*drhods2 - dy_eval(1, i_br)*d2rhods2m
                    d2B_r(1) = d2qua_dr2*drhods - 2.0_dp*dqua_dr*d2rhods2m + &
                               qua*drhods*(3.0_dp/4.0_dp)/r_eval**2
                end if
            else
                B_r = 0.0_dp
                dB_r = 0.0_dp
                d2B_r = 0.0_dp
            end if
        end if

        ! Interpolation of B_\vartheta and B_\varphi (flux functions)
        if (.not. bcovar_tp_batch_spline_ready) then
            error stop "splint_boozer_coord: Bcovar_tp batch spline not initialized"
        end if

        call evaluate_batch_splines_1d_der2(bcovar_tp_batch_spline, rho_tor, y1d, &
                                            dy1d, d2y1d)
        B_vartheta_B = y1d(1)
        dB_vartheta_B = dy1d(1)
        B_varphi_B = y1d(2)
        dB_varphi_B = dy1d(2)
        dB_vartheta_B = dB_vartheta_B*drhods
        dB_varphi_B = dB_varphi_B*drhods
        if (mode_secders > 0) then
            d2B_vartheta_B = d2y1d(1)*drhods2 - dy1d(1)*d2rhods2m
            d2B_varphi_B = d2y1d(2)*drhods2 - dy1d(2)*d2rhods2m
        else
            d2B_vartheta_B = 0.0_dp
            d2B_varphi_B = 0.0_dp
        end if

    end subroutine splint_boozer_coord

    subroutine compute_boozer_data
        ! Computes Boozer coordinate transformations and magnetic field data
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, &
                                          use_B_r
        use binsrc_sub, only: binsrc
        use plag_coeff_sub, only: plag_coeff
        use spline_vmec_sub

        implicit none

        real(dp), parameter :: s_min = 1.0e-6_dp, rho_min = sqrt(s_min)

        integer :: i, i_rho, i_theta, i_phi, npoilag, nder, nshift
        integer :: ibeg, iend, nqua
        real(dp) :: s, theta, varphi, A_theta, A_phi
        real(dp) :: dA_theta_ds, dA_phi_ds, aiota
        real(dp) :: sqg, alam, dl_ds, dl_dt, dl_dp
        real(dp) :: Bctrvr_vartheta, Bctrvr_varphi
        real(dp) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi
        real(dp) :: Bcovar_vartheta_B, Bcovar_varphi_B
        real(dp) :: denomjac, G00, Gbeg, aper
        real(dp) :: per_theta, per_phi, gridcellnum
        real(dp), allocatable :: wint_t(:), wint_p(:), theta_V(:), theta_B(:)
        real(dp), allocatable :: phi_V(:), phi_B(:), aiota_arr(:), rho_tor(:)
        real(dp), allocatable :: Bcovar_theta_V(:, :), Bcovar_varphi_V(:, :)
        real(dp), allocatable :: bmod_Vg(:, :), alam_2D(:, :)
        real(dp), allocatable :: sqrt_g_ss(:, :)
        real(dp), allocatable :: deltheta_BV_Vg(:, :), delphi_BV_Vg(:, :)
        real(dp), allocatable :: splcoe_t(:, :)
        real(dp), allocatable :: splcoe_p(:, :), coef(:, :)
        real(dp), allocatable :: perqua_t(:, :), perqua_p(:, :)
        real(dp), allocatable :: perqua_2D(:, :, :), Gfunc(:, :, :)
        real(dp), allocatable :: Bcovar_symfl(:, :, :, :)

        nqua = 7
        gridcellnum = real((n_theta_B - 1)*(n_phi_B - 1), dp)

        npoilag = ns_tp_B + 1
        nder = 0
        nshift = npoilag/2

        print *, 'Transforming to Boozer coordinates'

        if (use_B_r) then
            print *, 'B_r is computed'
        else
            print *, 'B_r is not computed'
        end if

        G00 = 0.0_dp

        allocate (rho_tor(ns_B))
        allocate (aiota_arr(1))
        allocate (Gfunc(1, 1, 1))
        allocate (Bcovar_symfl(1, 1, 1, 1))
        if (use_B_r) then
            deallocate (aiota_arr, Gfunc, Bcovar_symfl)
            allocate (aiota_arr(ns_B))
            allocate (Gfunc(ns_B, n_theta_B, n_phi_B))
            allocate (Bcovar_symfl(3, ns_B, n_theta_B, n_phi_B))
        end if

        allocate (Bcovar_theta_V(n_theta_B, n_phi_B))
        allocate (Bcovar_varphi_V(n_theta_B, n_phi_B))
        allocate (bmod_Vg(n_theta_B, n_phi_B))
        allocate (alam_2D(n_theta_B, n_phi_B))
        allocate (sqrt_g_ss(n_theta_B, n_phi_B))
        allocate (deltheta_BV_Vg(n_theta_B, n_phi_B))
        allocate (delphi_BV_Vg(n_theta_B, n_phi_B))
        allocate (wint_t(0:ns_tp_B), wint_p(0:ns_tp_B))
        allocate (coef(0:nder, npoilag))
        allocate (theta_V(2 - n_theta_B:2*n_theta_B - 1))
        allocate (theta_B(2 - n_theta_B:2*n_theta_B - 1))
        allocate (phi_V(2 - n_phi_B:2*n_phi_B - 1))
        allocate (phi_B(2 - n_phi_B:2*n_phi_B - 1))
        allocate (perqua_t(nqua, 2 - n_theta_B:2*n_theta_B - 1))
        allocate (perqua_p(nqua, 2 - n_phi_B:2*n_phi_B - 1))
        allocate (perqua_2D(nqua, n_theta_B, n_phi_B))

        allocate (splcoe_t(0:ns_tp_B, n_theta_B))
        allocate (splcoe_p(0:ns_tp_B, n_phi_B))

! allocate data arrays for Boozer data:
        if (.not. allocated(s_Bcovar_tp_B)) &
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, ns_B))

        ! Allocate module-level grids
        call ensure_grid_3d(bmod_grid, ns_B, n_theta_B, n_phi_B)
        call ensure_grid_3d(sqrt_g_ss_grid, ns_B, n_theta_B, n_phi_B)
        if (use_B_r) call ensure_grid_3d(br_grid, ns_B, n_theta_B, n_phi_B)

        do i = 0, ns_tp_B
            wint_t(i) = h_theta_B**(i + 1)/real(i + 1, dp)
            wint_p(i) = h_phi_B**(i + 1)/real(i + 1, dp)
        end do

        ! Set theta_V and phi_V linear, with value 0 at index 1 and stepsize h.
        ! Then expand this in both directions beyond 1:n_theta_B.
        do i_theta = 1, n_theta_B
            theta_V(i_theta) = real(i_theta - 1, dp)*h_theta_B
        end do
        per_theta = real(n_theta_B - 1, dp)*h_theta_B
        theta_V(2 - n_theta_B:0) = theta_V(1:n_theta_B - 1) - per_theta
        theta_V(n_theta_B + 1:2*n_theta_B - 1) = theta_V(2:n_theta_B) + per_theta

        do i_phi = 1, n_phi_B
            phi_V(i_phi) = real(i_phi - 1, dp)*h_phi_B
        end do
        per_phi = real(n_phi_B - 1, dp)*h_phi_B
        phi_V(2 - n_phi_B:0) = phi_V(1:n_phi_B - 1) - per_phi
        phi_V(n_phi_B + 1:2*n_phi_B - 1) = phi_V(2:n_phi_B) + per_phi

        do i_rho = 1, ns_B
            rho_tor(i_rho) = max(real(i_rho - 1, dp)*hs_B, rho_min)
            s = rho_tor(i_rho)**2

            do i_theta = 1, n_theta_B
                theta = real(i_theta - 1, dp)*h_theta_B
                do i_phi = 1, n_phi_B
                    varphi = real(i_phi - 1, dp)*h_phi_B

                    call vmec_field_evaluate(s, theta, varphi, &
                                             A_theta, A_phi, dA_theta_ds, &
                                             dA_phi_ds, aiota, &
                                             sqg, alam, dl_ds, &
                                             dl_dt, dl_dp, &
                                             Bctrvr_vartheta, &
                                             Bctrvr_varphi, &
                                             Bcovar_r, Bcovar_vartheta, &
                                             Bcovar_varphi)

                    alam_2D(i_theta, i_phi) = alam
                    bmod_Vg(i_theta, i_phi) = &
                        sqrt(Bctrvr_vartheta*Bcovar_vartheta &
                             + Bctrvr_varphi*Bcovar_varphi)
                    Bcovar_theta_V(i_theta, i_phi) = Bcovar_vartheta*(1.0_dp + dl_dt)
                    Bcovar_varphi_V(i_theta, i_phi) = &
                        Bcovar_varphi + Bcovar_vartheta*dl_dp
                    sqrt_g_ss(i_theta, i_phi) = get_sqrt_g_ss_contravariant(s, &
                                                                            theta, &
                                                                            varphi)

                    perqua_2D(4, i_theta, i_phi) = Bcovar_r
                    perqua_2D(5, i_theta, i_phi) = Bcovar_vartheta
                    perqua_2D(6, i_theta, i_phi) = Bcovar_varphi
                end do
            end do

! covariant components $B_\vartheta$ and $B_\varphi$ of Boozer coordinates:
            Bcovar_vartheta_B = sum(Bcovar_theta_V(2:n_theta_B, 2:n_phi_B))/gridcellnum
            Bcovar_varphi_B = sum(Bcovar_varphi_V(2:n_theta_B, 2:n_phi_B))/gridcellnum
            s_Bcovar_tp_B(1, 1, i_rho) = Bcovar_vartheta_B
            s_Bcovar_tp_B(2, 1, i_rho) = Bcovar_varphi_B

            denomjac = 1.0_dp/(aiota*Bcovar_vartheta_B + Bcovar_varphi_B)
            Gbeg = G00 + Bcovar_vartheta_B*denomjac*alam_2D(1, 1)

            splcoe_t(0, :) = Bcovar_theta_V(:, 1)

            call spl_per(ns_tp_B, n_theta_B, h_theta_B, splcoe_t)

            delphi_BV_Vg(1, 1) = 0.0_dp
            do i_theta = 1, n_theta_B - 1
                delphi_BV_Vg(i_theta + 1, 1) = &
                    delphi_BV_Vg(i_theta, 1) &
                    + sum(wint_t*splcoe_t(:, i_theta))
            end do
            ! Remove linear increasing component from delphi_BV_Vg
            aper = (delphi_BV_Vg(n_theta_B, 1) &
                    - delphi_BV_Vg(1, 1))/real(n_theta_B - 1, dp)
            do i_theta = 2, n_theta_B
                delphi_BV_Vg(i_theta, 1) = &
                    delphi_BV_Vg(i_theta, 1) - aper*real(i_theta - 1, dp)
            end do

            do i_theta = 1, n_theta_B
                splcoe_p(0, :) = Bcovar_varphi_V(i_theta, :)

                call spl_per(ns_tp_B, n_phi_B, h_phi_B, splcoe_p)

                do i_phi = 1, n_phi_B - 1
                    delphi_BV_Vg(i_theta, i_phi + 1) = &
                        delphi_BV_Vg(i_theta, i_phi) &
                        + sum(wint_p*splcoe_p(:, i_phi))
                end do
                aper = (delphi_BV_Vg(i_theta, n_phi_B) &
                        - delphi_BV_Vg(i_theta, 1))/real(n_phi_B - 1, dp)
                do i_phi = 2, n_phi_B
                    delphi_BV_Vg(i_theta, i_phi) = &
                        delphi_BV_Vg(i_theta, i_phi) &
                        - aper*real(i_phi - 1, dp)
                end do
            end do

! difference between Boozer and VMEC toroidal angle,
! $\Delta \varphi_{BV}=\varphi_B-\varphi=G$:
            delphi_BV_Vg = denomjac*delphi_BV_Vg + Gbeg
! difference between Boozer and VMEC poloidal angle,
! $\Delta \vartheta_{BV}=\vartheta_B-\theta$:
            deltheta_BV_Vg = aiota*delphi_BV_Vg + alam_2D

! At this point, all quantities are specified on
! equidistant grid in VMEC angles $(\theta,\varphi)$

! Re-interpolate to equidistant grid in $(\vartheta_B,\varphi)$:

            do i_phi = 1, n_phi_B
                perqua_t(1, 1:n_theta_B) = deltheta_BV_Vg(:, i_phi)
                perqua_t(2, 1:n_theta_B) = delphi_BV_Vg(:, i_phi)
                perqua_t(3, 1:n_theta_B) = bmod_Vg(:, i_phi)
                perqua_t(4:6, 1:n_theta_B) = perqua_2D(4:6, :, i_phi)
                perqua_t(7, 1:n_theta_B) = sqrt_g_ss(:, i_phi)
                ! Extend range of theta values
                perqua_t(:, 2 - n_theta_B:0) = perqua_t(:, 1:n_theta_B - 1)
                perqua_t(:, n_theta_B + 1:2*n_theta_B - 1) = perqua_t(:, 2:n_theta_B)
                theta_B = theta_V + perqua_t(1, :)
                do i_theta = 1, n_theta_B

                    call binsrc(theta_B, 2 - n_theta_B, 2*n_theta_B - 1, &
                                theta_V(i_theta), i)

                    ibeg = i - nshift
                    iend = ibeg + ns_tp_B

                    call plag_coeff(npoilag, nder, theta_V(i_theta), &
                                    theta_B(ibeg:iend), coef)

                    perqua_2D(:, i_theta, i_phi) = matmul(perqua_t(:, ibeg:iend), &
                                                          coef(0, :))
                end do
            end do

! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi)$

! Re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$:

            do i_theta = 1, n_theta_B
                perqua_p(:, 1:n_phi_B) = perqua_2D(:, i_theta, :)
                perqua_p(:, 2 - n_phi_B:0) = perqua_p(:, 1:n_phi_B - 1)
                ! Extend range of phi values
                perqua_p(:, n_phi_B + 1:2*n_phi_B - 1) = perqua_p(:, 2:n_phi_B)
                phi_B = phi_V + perqua_p(2, :)
                do i_phi = 1, n_phi_B

                    call binsrc(phi_B, 2 - n_phi_B, 2*n_phi_B - 1, phi_V(i_phi), i)

                    ibeg = i - nshift
                    iend = ibeg + ns_tp_B

                    call plag_coeff(npoilag, nder, phi_V(i_phi), phi_B(ibeg:iend), coef)

                    perqua_2D(:, i_theta, i_phi) = matmul(perqua_p(:, ibeg:iend), &
                                                          coef(0, :))
                end do
            end do

            bmod_grid(i_rho, :, :) = perqua_2D(3, :, :)
            sqrt_g_ss_grid(i_rho, :, :) = perqua_2D(7, :, :)

! End re-interpolate to equidistant grid in $(\vartheta_B,\varphi_B)$

            if (use_B_r) then
                aiota_arr(i_rho) = aiota
                Gfunc(i_rho, :, :) = perqua_2D(2, :, :)
! covariant components $B_k$ in symmetry flux coordinates on equidistant grid of
! Boozer coordinates:
                Bcovar_symfl(:, i_rho, :, :) = perqua_2D(4:6, :, :)
            end if

        end do

        if (use_B_r) then
            call compute_br_from_symflux(rho_tor, aiota_arr, Gfunc, Bcovar_symfl)
            deallocate (aiota_arr, Gfunc, Bcovar_symfl)
        end if

        deallocate (Bcovar_theta_V, Bcovar_varphi_V, bmod_Vg, alam_2D, &
                    sqrt_g_ss, deltheta_BV_Vg, delphi_BV_Vg, &
                    wint_t, wint_p, coef, theta_V, theta_B, phi_V, phi_B, &
                    perqua_t, perqua_p, perqua_2D)

        print *, 'done'

    end subroutine compute_boozer_data

    !> Original VMEC field evaluation using global splines (boozer_converter interface)
    subroutine vmec_field_evaluate(s, theta, varphi, &
                                   A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                   sqg, alam, dl_ds, dl_dt, dl_dp, &
                                   Bctrvr_vartheta, Bctrvr_varphi, &
                                   Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
        use spline_vmec_sub, only: vmec_field
        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds
        real(dp), intent(out) :: aiota, sqg, alam
        real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
        real(dp), intent(out) :: Bctrvr_vartheta, Bctrvr_varphi
        real(dp), intent(out) :: Bcovar_r, Bcovar_vartheta, Bcovar_varphi

        ! Call the existing VMEC routine
        call vmec_field(s, theta, varphi, &
                        A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, &
                        Bctrvr_vartheta, Bctrvr_varphi, &
                        Bcovar_r, Bcovar_vartheta, Bcovar_varphi)
    end subroutine vmec_field_evaluate

    !> Computes sqrt(g^{ss}) which is the same for Boozer (s, theta_B, varphi_B)
    !> and VMEC (s, theta, varphi) coordinates
    function get_sqrt_g_ss_contravariant(s, theta, varphi) result(sqrt_g_ss)
        use spline_vmec_sub, only: splint_vmec_data
        use spline_vmec_sub, only: metric_tensor_vmec
        real(dp), intent(in) :: s, theta, varphi
        real(dp) :: sqrt_g_ss

        real(dp) :: dummy(10)
        real(dp) :: R, dR_ds, dR_dtheta, dR_dphi
        real(dp) :: dZ_ds, dZ_dtheta, dZ_dphi

        real(dp) :: g_vmec(3, 3), sqrt_g_vmec

        call splint_vmec_data(s, theta, varphi, &
                              dummy(1), dummy(2), dummy(3), dummy(4), dummy(5), &
                              R, &
                              dummy(6), dummy(7), &
                              dR_ds, dR_dtheta, dR_dphi, &
                              dZ_ds, dZ_dtheta, dZ_dphi, &
                              dummy(8), dummy(9), dummy(10))
        call metric_tensor_vmec(R, dR_ds, dR_dtheta, dR_dphi, &
                                dZ_ds, dZ_dtheta, dZ_dphi, g_vmec, sqrt_g_vmec)

        !> contravariant metric component g^{ss} via cofactors of covariant components
        sqrt_g_ss = sqrt(g_vmec(2, 2)*g_vmec(3, 3) - g_vmec(2, 3)**2.0_dp) &
                    /abs(sqrt_g_vmec)
    end function get_sqrt_g_ss_contravariant

    !> Compute radial covariant magnetic field B_rho from symmetry flux coordinates
    subroutine compute_br_from_symflux(rho_tor, aiota_arr, Gfunc, Bcovar_symfl)
        use boozer_coordinates_mod, only: ns_B, n_theta_B, n_phi_B
        use plag_coeff_sub, only: plag_coeff

        real(dp), intent(in) :: rho_tor(:)
        real(dp), intent(in) :: aiota_arr(:)
        real(dp), intent(in) :: Gfunc(:, :, :)
        real(dp), intent(in) :: Bcovar_symfl(:, :, :, :)

        integer, parameter :: NPOILAG = 5
        integer, parameter :: NDER = 1

        integer :: i_rho, i_phi, ibeg, iend, nshift
        real(dp) :: coef(0:NDER, NPOILAG)

        nshift = NPOILAG/2

        do i_rho = 1, ns_B
            ibeg = i_rho - nshift
            iend = ibeg + NPOILAG - 1
            if (ibeg < 1) then
                ibeg = 1
                iend = ibeg + NPOILAG - 1
            else if (iend > ns_B) then
                iend = ns_B
                ibeg = iend - NPOILAG + 1
            end if

            call plag_coeff(NPOILAG, NDER, rho_tor(i_rho), rho_tor(ibeg:iend), coef)

            ! Compute B_rho (we spline covariant component B_rho instead of B_s)
            do i_phi = 1, n_phi_B
                br_grid(i_rho, :, i_phi) = &
                    2.0_dp*rho_tor(i_rho)*Bcovar_symfl(1, i_rho, :, i_phi) &
                    - matmul(coef(1, :)*aiota_arr(ibeg:iend), Gfunc(ibeg:iend, &
                                                                    :, i_phi)) &
                    *Bcovar_symfl(2, i_rho, :, i_phi) &
                    - matmul(coef(1, :), Gfunc(ibeg:iend, :, i_phi)) &
                    *Bcovar_symfl(3, i_rho, :, i_phi)
            end do
        end do

    end subroutine compute_br_from_symflux

    !> Ensure 3D grid is allocated with correct dimensions
    subroutine ensure_grid_3d(grid, n1, n2, n3)
        real(dp), allocatable, intent(inout) :: grid(:, :, :)
        integer, intent(in) :: n1, n2, n3

        if (.not. allocated(grid)) then
            allocate (grid(n1, n2, n3))
        else if (any(shape(grid) /= [n1, n2, n3])) then
            deallocate (grid)
            allocate (grid(n1, n2, n3))
        end if
    end subroutine ensure_grid_3d

    subroutine reset_boozer_batch_splines
        if (aphi_batch_spline_ready) then
            call destroy_batch_splines_1d(aphi_batch_spline)
            aphi_batch_spline_ready = .false.
        end if
        if (bcovar_tp_batch_spline_ready) then
            call destroy_batch_splines_1d(bcovar_tp_batch_spline)
            bcovar_tp_batch_spline_ready = .false.
        end if
        if (field3d_batch_spline_ready) then
            call destroy_batch_splines_3d(field3d_batch_spline)
            field3d_batch_spline_ready = .false.
            field3d_num_quantities = 0
        end if
        if (allocated(bmod_grid)) deallocate (bmod_grid)
        if (allocated(sqrt_g_ss_grid)) deallocate (sqrt_g_ss_grid)
        if (allocated(br_grid)) deallocate (br_grid)
    end subroutine reset_boozer_batch_splines

    subroutine build_boozer_aphi_batch_spline
        use vector_potentail_mod, only: ns, hs, sA_phi
        use new_vmec_stuff_mod, only: ns_A

        integer :: order

        if (aphi_batch_spline_ready) then
            call destroy_batch_splines_1d(aphi_batch_spline)
            aphi_batch_spline_ready = .false.
        end if

        order = ns_A
        if (order < 3 .or. order > 5) then
            error stop "build_boozer_aphi_batch_spline: spline order must be 3..5"
        end if

        aphi_batch_spline%order = order
        aphi_batch_spline%num_points = ns
        aphi_batch_spline%periodic = .false.
        aphi_batch_spline%x_min = 0.0_dp
        aphi_batch_spline%h_step = hs
        aphi_batch_spline%num_quantities = 1

        allocate (aphi_batch_spline%coeff(1, 0:order, ns))
        aphi_batch_spline%coeff(1, 0:order, :) = sA_phi(1:order + 1, :)

        aphi_batch_spline_ready = .true.
    end subroutine build_boozer_aphi_batch_spline

    subroutine build_boozer_bcovar_tp_batch_spline
        use boozer_coordinates_mod, only: ns_s_B, ns_B, hs_B, s_Bcovar_tp_B

        integer :: order
        real(dp) :: x_min, x_max
        real(dp), allocatable :: y_batch(:, :)

        if (bcovar_tp_batch_spline_ready) then
            call destroy_batch_splines_1d(bcovar_tp_batch_spline)
            bcovar_tp_batch_spline_ready = .false.
        end if

        order = ns_s_B
        if (order < 3 .or. order > 5) then
            error stop "build_boozer_bcovar_tp_batch_spline: spline order must be 3..5"
        end if

        x_min = 0.0_dp
        x_max = hs_B*real(ns_B - 1, dp)

        allocate (y_batch(ns_B, 2))
        y_batch(:, 1) = s_Bcovar_tp_B(1, 1, :)
        y_batch(:, 2) = s_Bcovar_tp_B(2, 1, :)

        call construct_batch_splines_1d(x_min, x_max, y_batch, order, .false., &
                                        bcovar_tp_batch_spline)
        bcovar_tp_batch_spline_ready = .true.
        deallocate (y_batch)
    end subroutine build_boozer_bcovar_tp_batch_spline

    subroutine build_boozer_field3d_batch_spline
        ! Combined 3D field batch spline: Bmod, sqrt_g_ss, optionally Br
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, use_B_r

        real(dp) :: x_min(3), x_max(3)
        real(dp), allocatable :: y_batch(:, :, :, :)
        integer :: order(3), nq, i_br
        logical :: periodic(3)

        if (.not. allocated(bmod_grid)) then
            error stop "build_boozer_field3d_batch_spline: bmod_grid not allocated"
        end if
        if (.not. allocated(sqrt_g_ss_grid)) then
            error stop "build_boozer_field3d_batch_spline: sqrt_g_ss_grid not allocated"
        end if
        if (use_B_r .and. .not. allocated(br_grid)) then
            error stop "build_boozer_field3d_batch_spline: br_grid not allocated"
        end if

        if (field3d_batch_spline_ready) then
            call destroy_batch_splines_3d(field3d_batch_spline)
            field3d_batch_spline_ready = .false.
            field3d_num_quantities = 0
        end if

        order = [ns_s_B, ns_tp_B, ns_tp_B]
        if (any(order < 3) .or. any(order > 5)) then
            error stop "build_boozer_field3d_batch_spline: spline order must be 3..5"
        end if

        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max(1) = hs_B*real(ns_B - 1, dp)
        x_max(2) = h_theta_B*real(n_theta_B - 1, dp)
        x_max(3) = h_phi_B*real(n_phi_B - 1, dp)

        periodic = [.false., .true., .true.]

        nq = 2  ! Bmod, sqrt_g_ss
        if (use_B_r) then
            nq = nq + 1
            i_br = nq
        end if

        allocate (y_batch(ns_B, n_theta_B, n_phi_B, nq))
        y_batch(:, :, :, 1) = bmod_grid(:, :, :)
        y_batch(:, :, :, 2) = sqrt_g_ss_grid(:, :, :)
        if (use_B_r) then
            y_batch(:, :, :, i_br) = br_grid(:, :, :)
        end if

        call construct_batch_splines_3d(x_min, x_max, y_batch, order, periodic, &
                                        field3d_batch_spline)
        field3d_batch_spline_ready = .true.
        field3d_num_quantities = nq
        deallocate (y_batch)
    end subroutine build_boozer_field3d_batch_spline

    !> Initialize Boozer splines from a libneo Boozer chartmap NetCDF.
    !>
    !> The chartmap stores Bmod on a (rho, theta, zeta) grid, plus
    !> A_phi(s), B_theta(rho), B_phi(rho) as radial surface functions.
    !> All quantities are in CGS-Gaussian units (G, cm, G*cm, G*cm^2).
    !>
    !> After this call the same batch-spline state is populated as after
    !> get_boozer_coordinates, so splint_boozer_coord can be called directly.
    subroutine get_boozer_coordinates_from_chartmap(chartmap_file)
        use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
        use netcdf, only: nf90_get_att, nf90_inq_varid, nf90_noerr, nf90_global
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper, rmajor, ns_s, ns_tp
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, use_B_r

        character(len=*), intent(in) :: chartmap_file

        integer :: ncid, varid, ncstat
        integer :: nrho, ntheta, nzeta
        integer :: nfp_int
        real(dp) :: torflux_val
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: Bmod_raw(:, :, :)
        real(dp), allocatable :: A_phi(:), B_theta(:), B_phi(:)
        real(dp), allocatable :: xyz_outer(:, :)
        real(dp), allocatable :: aphi_batch(:, :)
        real(dp) :: rho_min, rho_max, s_min, s_max, rmajor_cm
        real(dp) :: h_rho

        call nc_open(trim(chartmap_file), ncid)

        call nc_inq_dim(ncid, 'rho', nrho)
        call nc_inq_dim(ncid, 'theta', ntheta)
        call nc_inq_dim(ncid, 'zeta', nzeta)

        allocate (rho(nrho), theta(ntheta), zeta(nzeta))
        call nc_get(ncid, 'rho', rho)
        call nc_get(ncid, 'theta', theta)
        call nc_get(ncid, 'zeta', zeta)

        nfp_int = 1
        call nc_get(ncid, 'num_field_periods', nfp_int)

        torflux_val = 0.0_dp
        if (nf90_get_att(ncid, nf90_global, 'torflux', torflux_val) /= nf90_noerr) then
            torflux_val = 0.0_dp
        end if

        allocate (A_phi(nrho))
        allocate (B_theta(nrho))
        allocate (B_phi(nrho))
        call nc_get(ncid, 'A_phi', A_phi)
        call nc_get(ncid, 'B_theta', B_theta)
        call nc_get(ncid, 'B_phi', B_phi)

        allocate (Bmod_raw(nrho, ntheta, nzeta))
        call nc_get(ncid, 'Bmod', Bmod_raw)

        ! Estimate major radius from x,y coordinates on outermost surface (units: cm).
        ! Falls back to 1 cm if the variables are absent.
        allocate (xyz_outer(ntheta*nzeta, 2))
        ncstat = nf90_inq_varid(ncid, 'x', varid)
        if (ncstat == nf90_noerr) then
            block
                real(dp), allocatable :: x3(:, :, :), y3(:, :, :)
                integer :: it, iz, idx
                allocate (x3(nzeta, ntheta, nrho), y3(nzeta, ntheta, nrho))
                call nc_get(ncid, 'x', x3)
                call nc_get(ncid, 'y', y3)
                idx = 0
                do iz = 1, nzeta
                    do it = 1, ntheta
                        idx = idx + 1
                        xyz_outer(idx, 1) = x3(iz, it, nrho)
                        xyz_outer(idx, 2) = y3(iz, it, nrho)
                    end do
                end do
                deallocate (x3, y3)
                rmajor_cm = sum(sqrt(xyz_outer(:, 1)**2 + xyz_outer(:, 2)**2)) &
                            /real(ntheta*nzeta, dp)
            end block
        else
            rmajor_cm = 1.0_dp
        end if
        deallocate (xyz_outer)

        call nc_close(ncid)

        torflux = torflux_val
        nper = nfp_int
        rmajor = rmajor_cm

        ns_s = 3
        ns_tp = 3
        ns_s_B = ns_s
        ns_tp_B = ns_tp
        ns_B = nrho
        n_theta_B = ntheta
        n_phi_B = nzeta
        use_B_r = .false.

        rho_min = rho(1)
        rho_max = rho(nrho)
        h_rho = (rho_max - rho_min)/real(nrho - 1, dp)

        hs_B = h_rho
        h_theta_B = (theta(ntheta) - theta(1))/real(ntheta - 1, dp)
        h_phi_B = (zeta(nzeta) - zeta(1))/real(nzeta - 1, dp)

        call reset_boozer_batch_splines()

        call ensure_grid_3d(bmod_grid, nrho, ntheta, nzeta)
        call ensure_grid_3d(sqrt_g_ss_grid, nrho, ntheta, nzeta)
        bmod_grid = Bmod_raw
        sqrt_g_ss_grid = 0.0_dp

        if (.not. allocated(s_Bcovar_tp_B)) then
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, nrho))
        else if (any(shape(s_Bcovar_tp_B) /= [2, ns_s_B + 1, nrho])) then
            deallocate (s_Bcovar_tp_B)
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, nrho))
        end if
        s_Bcovar_tp_B(1, 1, :) = B_theta
        s_Bcovar_tp_B(2, 1, :) = B_phi

        s_min = rho(1)**2
        s_max = rho(nrho)**2

        allocate (aphi_batch(nrho, 1))
        aphi_batch(:, 1) = A_phi
        call construct_batch_splines_1d(s_min, s_max, aphi_batch, 3, .false., &
                                        aphi_batch_spline)
        aphi_batch_spline_ready = .true.
        deallocate (aphi_batch)

        block
            real(dp), allocatable :: bcovar_batch(:, :)
            allocate (bcovar_batch(nrho, 2))
            bcovar_batch(:, 1) = B_theta
            bcovar_batch(:, 2) = B_phi
            call construct_batch_splines_1d(rho_min, rho_max, bcovar_batch, 3, .false., &
                                            bcovar_tp_batch_spline)
            bcovar_tp_batch_spline_ready = .true.
            deallocate (bcovar_batch)
        end block

        block
            real(dp) :: x3_min(3), x3_max(3)
            real(dp), allocatable :: y3(:, :, :, :)
            integer :: order3(3)
            logical :: periodic3(3)

            x3_min = [rho_min, theta(1), zeta(1)]
            x3_max = [rho_max, theta(ntheta), zeta(nzeta)]
            order3 = [3, 3, 3]
            periodic3 = [.false., .true., .true.]

            if (field3d_batch_spline_ready) then
                call destroy_batch_splines_3d(field3d_batch_spline)
                field3d_batch_spline_ready = .false.
                field3d_num_quantities = 0
            end if

            allocate (y3(nrho, ntheta, nzeta, 2))
            y3(:, :, :, 1) = bmod_grid
            y3(:, :, :, 2) = sqrt_g_ss_grid
            call construct_batch_splines_3d(x3_min, x3_max, y3, order3, periodic3, &
                                            field3d_batch_spline)
            field3d_batch_spline_ready = .true.
            field3d_num_quantities = 2
            deallocate (y3)
        end block

        deallocate (rho, theta, zeta, Bmod_raw, A_phi, B_theta, B_phi)

    end subroutine get_boozer_coordinates_from_chartmap

    !> Initialize Boozer splines from a booz_xform boozmn NetCDF.
    !>
    !> Reads the Fourier harmonics bmnc_b, iota_b, buco_b, bvco_b, phi_b
    !> and the mode arrays ixm_b, ixn_b; then evaluates Bmod on a uniform
    !> (nrho, ntheta, nzeta) grid by Fourier summation.  The covariant
    !> components B_theta, B_phi and A_phi are recovered from the surface
    !> functions buco_b, bvco_b, phi_b on the full grid.
    !>
    !> nrho, ntheta, nzeta: output grid resolution (optional, defaults 30/48/96).
    subroutine get_boozer_coordinates_from_boozmn(boozmn_file, nrho_in, ntheta_in, nzeta_in)
        use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper, rmajor, ns_s, ns_tp
        use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B, &
                                          hs_B, h_theta_B, h_phi_B, &
                                          s_Bcovar_tp_B, use_B_r

        character(len=*), intent(in) :: boozmn_file
        integer, intent(in), optional :: nrho_in, ntheta_in, nzeta_in

        integer :: nrho, ntheta, nzeta
        integer :: ncid
        integer :: ns, nmn, nsurf_computed
        integer :: nfp_int
        integer, allocatable :: jlist(:), ixm(:), ixn(:)
        real(dp), allocatable :: iota_full(:), buco_full(:), bvco_full(:), phi_full(:)
        ! bmnc_h shape is (nmn, nsurf_computed): Fortran reads NetCDF (comput_surfs,mn_mode)
        ! with last NetCDF dim = first Fortran dim
        real(dp), allocatable :: bmnc_h(:, :)
        real(dp), allocatable :: rho_half(:), s_half(:)
        real(dp), allocatable :: rho_out(:), s_out(:)
        real(dp), allocatable :: B_theta(:), B_phi(:), A_phi_out(:)
        real(dp), allocatable :: bmnc_out(:, :)
        real(dp), allocatable :: rmnc_h(:, :)
        real(dp), allocatable :: Bmod(:, :, :)
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: aphi_batch(:, :)
        integer :: ir, k, it, iz, mn, mn00
        real(dp) :: angle, torflux_si, torflux_cgs, rmajor_m
        real(dp) :: rho_min, rho_max, s_min, s_max
        real(dp), parameter :: MU0 = 4.0e-7_dp*3.14159265358979_dp
        real(dp), parameter :: GAUSS_CM2_PER_TM2 = 1.0e8_dp
        real(dp), parameter :: GAUSS_CM_PER_TM = 1.0e6_dp
        real(dp), parameter :: GAUSS_PER_T = 1.0e4_dp

        nrho = 30
        ntheta = 48
        nzeta = 96
        if (present(nrho_in)) nrho = nrho_in
        if (present(ntheta_in)) ntheta = ntheta_in
        if (present(nzeta_in)) nzeta = nzeta_in

        call nc_open(trim(boozmn_file), ncid)
        call nc_get(ncid, 'ns_b', ns)
        call nc_inq_dim(ncid, 'ixm_b', nmn)
        call nc_inq_dim(ncid, 'jlist', nsurf_computed)

        allocate (jlist(nsurf_computed))
        allocate (ixm(nmn), ixn(nmn))
        allocate (iota_full(ns), buco_full(ns), bvco_full(ns), phi_full(ns))
        allocate (bmnc_h(nmn, nsurf_computed))
        allocate (rmnc_h(nmn, nsurf_computed))

        call nc_get(ncid, 'nfp_b', nfp_int)
        call nc_get(ncid, 'jlist', jlist)
        call nc_get(ncid, 'ixm_b', ixm)
        call nc_get(ncid, 'ixn_b', ixn)
        call nc_get(ncid, 'iota_b', iota_full)
        call nc_get(ncid, 'buco_b', buco_full)
        call nc_get(ncid, 'bvco_b', bvco_full)
        call nc_get(ncid, 'phi_b', phi_full)
        call nc_get(ncid, 'bmnc_b', bmnc_h)
        call nc_get(ncid, 'rmnc_b', rmnc_h)
        call nc_close(ncid)

        allocate (s_half(nsurf_computed))
        allocate (rho_half(nsurf_computed))
        do k = 1, nsurf_computed
            s_half(k) = (real(jlist(k), dp) - 1.5_dp)/real(ns - 1, dp)
        end do
        rho_half = sqrt(s_half)

        torflux_si = -phi_full(ns)/TWOPI
        torflux_cgs = torflux_si*GAUSS_CM2_PER_TM2

        rho_min = 1.0e-3_dp
        rho_max = 1.0_dp
        allocate (rho_out(nrho))
        allocate (s_out(nrho))
        allocate (theta(ntheta))
        allocate (zeta(nzeta))
        allocate (B_theta(nrho))
        allocate (B_phi(nrho))
        allocate (A_phi_out(nrho))
        allocate (bmnc_out(nrho, nmn))
        allocate (Bmod(nrho, ntheta, nzeta))

        do k = 1, nrho
            rho_out(k) = rho_min + (rho_max - rho_min)*real(k - 1, dp)/real(nrho - 1, dp)
        end do
        s_out = rho_out**2

        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = TWOPI/real(nfp_int, dp)*real(iz - 1, dp)/real(nzeta, dp)
        end do

        call boozmn_interp_modes(bmnc_h, rho_half, rho_out, ixm, nmn, nsurf_computed, &
                                 nrho, bmnc_out)
        call boozmn_interp_1d(buco_full(jlist), rho_half, rho_out, nsurf_computed, nrho, &
                              B_theta)
        call boozmn_interp_1d(bvco_full(jlist), rho_half, rho_out, nsurf_computed, nrho, &
                              B_phi)

        call boozmn_iota_integral(iota_full(jlist), rho_half, s_out, nsurf_computed, nrho, &
                                  torflux_si, A_phi_out)

        B_theta = B_theta*GAUSS_CM_PER_TM
        B_phi = B_phi*GAUSS_CM_PER_TM
        A_phi_out = A_phi_out*GAUSS_CM2_PER_TM2

        do ir = 1, nrho
            do it = 1, ntheta
                do iz = 1, nzeta
                    Bmod(ir, it, iz) = 0.0_dp
                    do mn = 1, nmn
                        angle = real(ixm(mn), dp)*theta(it) &
                                - real(ixn(mn), dp)*zeta(iz)
                        Bmod(ir, it, iz) = Bmod(ir, it, iz) + bmnc_out(ir, mn)*cos(angle)
                    end do
                    Bmod(ir, it, iz) = Bmod(ir, it, iz)*GAUSS_PER_T
                end do
            end do
        end do

        mn00 = 0
        do mn = 1, nmn
            if (ixm(mn) == 0 .and. ixn(mn) == 0) then
                mn00 = mn
                exit
            end if
        end do
        if (mn00 > 0) then
            rmajor_m = rmnc_h(mn00, nsurf_computed)
        else
            rmajor_m = 1.0_dp
        end if

        torflux = torflux_cgs
        nper = nfp_int
        rmajor = rmajor_m*1.0e2_dp

        ns_s = 3
        ns_tp = 3
        ns_s_B = ns_s
        ns_tp_B = ns_tp
        ns_B = nrho
        n_theta_B = ntheta
        n_phi_B = nzeta
        use_B_r = .false.

        hs_B = (rho_max - rho_min)/real(nrho - 1, dp)
        h_theta_B = theta(2) - theta(1)
        h_phi_B = zeta(2) - zeta(1)

        call reset_boozer_batch_splines()

        call ensure_grid_3d(bmod_grid, nrho, ntheta, nzeta)
        call ensure_grid_3d(sqrt_g_ss_grid, nrho, ntheta, nzeta)
        bmod_grid = Bmod
        sqrt_g_ss_grid = 0.0_dp

        if (.not. allocated(s_Bcovar_tp_B)) then
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, nrho))
        else if (any(shape(s_Bcovar_tp_B) /= [2, ns_s_B + 1, nrho])) then
            deallocate (s_Bcovar_tp_B)
            allocate (s_Bcovar_tp_B(2, ns_s_B + 1, nrho))
        end if
        s_Bcovar_tp_B(1, 1, :) = B_theta
        s_Bcovar_tp_B(2, 1, :) = B_phi

        s_min = rho_min**2
        s_max = rho_max**2
        allocate (aphi_batch(nrho, 1))
        aphi_batch(:, 1) = A_phi_out
        call construct_batch_splines_1d(s_min, s_max, aphi_batch, 3, .false., &
                                        aphi_batch_spline)
        aphi_batch_spline_ready = .true.
        deallocate (aphi_batch)

        block
            real(dp), allocatable :: bcovar_batch(:, :)
            allocate (bcovar_batch(nrho, 2))
            bcovar_batch(:, 1) = B_theta
            bcovar_batch(:, 2) = B_phi
            call construct_batch_splines_1d(rho_min, rho_max, bcovar_batch, 3, .false., &
                                            bcovar_tp_batch_spline)
            bcovar_tp_batch_spline_ready = .true.
            deallocate (bcovar_batch)
        end block

        block
            real(dp) :: x3_min(3), x3_max(3)
            real(dp), allocatable :: y3(:, :, :, :)
            integer :: order3(3)
            logical :: periodic3(3)

            x3_min = [rho_min, theta(1), zeta(1)]
            x3_max = [rho_max, theta(ntheta), zeta(nzeta)]
            order3 = [3, 3, 3]
            periodic3 = [.false., .true., .true.]

            if (field3d_batch_spline_ready) then
                call destroy_batch_splines_3d(field3d_batch_spline)
                field3d_batch_spline_ready = .false.
                field3d_num_quantities = 0
            end if

            allocate (y3(nrho, ntheta, nzeta, 2))
            y3(:, :, :, 1) = bmod_grid
            y3(:, :, :, 2) = sqrt_g_ss_grid
            call construct_batch_splines_3d(x3_min, x3_max, y3, order3, periodic3, &
                                            field3d_batch_spline)
            field3d_batch_spline_ready = .true.
            field3d_num_quantities = 2
            deallocate (y3)
        end block

        deallocate (jlist, ixm, ixn, iota_full, buco_full, bvco_full, phi_full)
        deallocate (bmnc_h, rmnc_h, s_half, rho_half, rho_out, s_out)
        deallocate (theta, zeta, B_theta, B_phi, A_phi_out, bmnc_out, Bmod)

    end subroutine get_boozer_coordinates_from_boozmn

    !> Interpolate boozmn Fourier coefficients from half grid to output rho grid.
    !> bmnc_h is shaped (nmn, nsurf): Fortran order from NetCDF (comput_surfs, mn_mode).
    subroutine boozmn_interp_modes(bmnc_h, rho_half, rho_out, ixm, nmn, nsurf, nrho_out, &
                                   bmnc_out)
        integer, intent(in) :: nmn, nsurf, nrho_out
        real(dp), intent(in) :: bmnc_h(nmn, nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: rho_out(nrho_out)
        integer, intent(in) :: ixm(nmn)
        real(dp), intent(out) :: bmnc_out(nrho_out, nmn)

        integer :: ir, mn, k
        real(dp) :: rho, frac
        real(dp) :: ratio

        do ir = 1, nrho_out
            rho = rho_out(ir)
            if (rho <= rho_half(1)) then
                do mn = 1, nmn
                    if (rho_half(1) > 0.0_dp) then
                        ratio = rho/rho_half(1)
                        bmnc_out(ir, mn) = bmnc_h(mn, 1)*ratio**min(ixm(mn), 50)
                    else
                        bmnc_out(ir, mn) = bmnc_h(mn, 1)
                    end if
                end do
            else if (rho >= rho_half(nsurf)) then
                do mn = 1, nmn
                    bmnc_out(ir, mn) = bmnc_h(mn, nsurf)
                end do
            else
                k = 1
                do while (k < nsurf .and. rho_half(k + 1) < rho)
                    k = k + 1
                end do
                if (rho_half(k + 1) > rho_half(k)) then
                    frac = (rho - rho_half(k))/(rho_half(k + 1) - rho_half(k))
                else
                    frac = 0.0_dp
                end if
                do mn = 1, nmn
                    bmnc_out(ir, mn) = (1.0_dp - frac)*bmnc_h(mn, k) &
                                       + frac*bmnc_h(mn, k + 1)
                end do
            end if
        end do
    end subroutine boozmn_interp_modes

    !> Interpolate a 1D radial profile from half grid to output rho grid.
    subroutine boozmn_interp_1d(vals_h, rho_half, rho_out, nsurf, nrho_out, vals_out)
        integer, intent(in) :: nsurf, nrho_out
        real(dp), intent(in) :: vals_h(nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: rho_out(nrho_out)
        real(dp), intent(out) :: vals_out(nrho_out)

        integer :: ir, k
        real(dp) :: rho, frac

        do ir = 1, nrho_out
            rho = rho_out(ir)
            if (rho <= rho_half(1)) then
                vals_out(ir) = vals_h(1)
            else if (rho >= rho_half(nsurf)) then
                vals_out(ir) = vals_h(nsurf)
            else
                k = 1
                do while (k < nsurf .and. rho_half(k + 1) < rho)
                    k = k + 1
                end do
                if (rho_half(k + 1) > rho_half(k)) then
                    frac = (rho - rho_half(k))/(rho_half(k + 1) - rho_half(k))
                else
                    frac = 0.0_dp
                end if
                vals_out(ir) = (1.0_dp - frac)*vals_h(k) + frac*vals_h(k + 1)
            end if
        end do
    end subroutine boozmn_interp_1d

    !> Compute A_phi(s) by integrating -torflux_si * iota(s) over s.
    !> A_phi(s) = -torflux_si * integral_0^s iota(s') ds'
    !> Units: T*m^2 (same as torflux_si).
    subroutine boozmn_iota_integral(iota_h, rho_half, s_out, nsurf, nrho_out, torflux_si, &
                                    A_phi_out)
        integer, intent(in) :: nsurf, nrho_out
        real(dp), intent(in) :: iota_h(nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: s_out(nrho_out)
        real(dp), intent(in) :: torflux_si
        real(dp), intent(out) :: A_phi_out(nrho_out)

        integer :: ir, k
        real(dp) :: iota_at_s, s, s_half_sq(nsurf), cum(nsurf)

        s_half_sq = rho_half**2

        cum(1) = 0.0_dp
        do k = 2, nsurf
            cum(k) = cum(k - 1) + 0.5_dp*(iota_h(k - 1) + iota_h(k)) &
                     *(s_half_sq(k) - s_half_sq(k - 1))
        end do

        do ir = 1, nrho_out
            s = s_out(ir)
            if (s <= s_half_sq(1)) then
                iota_at_s = iota_h(1)
                A_phi_out(ir) = -torflux_si*iota_at_s*s
            else if (s >= s_half_sq(nsurf)) then
                A_phi_out(ir) = -torflux_si*(cum(nsurf) + &
                                             iota_h(nsurf)*(s - s_half_sq(nsurf)))
            else
                k = 1
                do while (k < nsurf .and. s_half_sq(k + 1) < s)
                    k = k + 1
                end do
                if (s_half_sq(k + 1) > s_half_sq(k)) then
                    iota_at_s = iota_h(k) + (iota_h(k + 1) - iota_h(k))* &
                                (s - s_half_sq(k))/(s_half_sq(k + 1) - s_half_sq(k))
                else
                    iota_at_s = iota_h(k)
                end if
                A_phi_out(ir) = -torflux_si*(cum(k) + 0.5_dp*(iota_h(k) + iota_at_s)* &
                                             (s - s_half_sq(k)))
            end if
        end do
    end subroutine boozmn_iota_integral

end module boozer_sub
