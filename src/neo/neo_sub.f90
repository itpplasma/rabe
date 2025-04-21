module neo_sub_mod

contains

    subroutine neo_init(npsi)
        ! Initialization Routine
        ! **********************************************************************
        ! Modules
        ! **********************************************************************
        use neo_input
        use neo_work
        use neo_exchange
        use neo_units
        use neo_control
        use neo_spline
        ! **********************************************************************
        ! Local Definitions
        ! **********************************************************************
        implicit none
        integer, intent(out) :: npsi
        integer :: imn
        ! **********************************************************************
        ! Read input from data file and allocate necessary arrays
        ! **********************************************************************
        write (w_us, *) 'before neo_read'
        call neo_read
        write (w_us, *) 'after  neo_read'
        ! **********************************************************************
        npsi = ns
        call neo_prep_b00
        ! **********************************************************************
        ! Allocate and prepare necessary arrays
        ! **********************************************************************
        write (w_us, *) 'before neo_prep'
        call neo_prep
        write (w_us, *) 'after  neo_prep'
        ! **********************************************************************
        ! Allocate and prepare spline along s
        ! **********************************************************************
        write (w_us, *) 'before neo_init_spline'
        call neo_init_spline()
        write (w_us, *) 'after  neo_init_spline'
        ! **********************************************************************
        ! Calculation of rt0 and bmref (innermost flux surface)
        ! as reference values
        ! **********************************************************************

        bmref = 0.0_dp
        do imn = 1, mnmax
            if (ixm(imn) .eq. 0 .and. ixn(imn) .eq. 0) then
                bmref = bmnc(1, imn)
                bmref_g = bmref
            end if
        end do

        if (bmref .eq. 0.0_dp) then
            write (w_us, *) ' NEO_INIT: Fatal problem setting rt0 or bmref', rt0, bmref
            stop
        end if

    end subroutine neo_init

    subroutine neo_prep_b00
        ! Preparation of b00 with splines
        !***********************************************************************
        ! Modules
        use nrtype
        use neo_input
        use neo_spline_b00
        use inter_interfaces, only: splinecof1, splinecof3, tf
        use neo_spline_data, only: lsw_linear_boozer

        integer :: i, j
        integer(I4B) :: sw1, sw2
        real(dp) :: m0, c1, cn

        real(dp), dimension(:), allocatable :: lambda
        integer, dimension(:), allocatable :: index_i

        ! collect the b00 in a vector
        allocate (lambda(ns))
        allocate (index_i(ns))
        allocate (a_b00(ns), b_b00(ns))
        allocate (c_b00(ns), d_b00(ns))

        do i = 1, ns
            do j = 1, mnmax
                if (ixm(j) .eq. 0 .and. ixn(j) .eq. 0) then
                    b00(i) = bmnc(i, j)
                    exit
                end if
            end do
        end do
        ! spline b00 along s
        ! boundary types (natural spline)
        sw1 = 2
        sw2 = 4
        ! input for test function for spline
        m0 = 0.0_dp
        ! boundary condition for spline
        c1 = 0.0_dp
        cn = 0.0_dp
        ! we use no smoothing for spline
        lambda = 1.0d0
        index_i = (/(i, i=1, ns)/)
        if (lsw_linear_boozer) then
            call splinecof1(es, b00, c1, cn, lambda, index_i, sw1, sw2, &
                & a_b00, b_b00, c_b00, d_b00, m0, tf)
        else
            call splinecof3(es, b00, c1, cn, lambda, index_i, sw1, sw2, &
                & a_b00, b_b00, c_b00, d_b00, m0, tf)
        end if

        deallocate (lambda)
        deallocate (index_i)

    end subroutine neo_prep_b00

    subroutine neo_init_spline()
        ! Initialization for splines along s
        ! **********************************************************************
        ! Modules
        ! **********************************************************************
        use nrtype
        use neo_spline_data
        use neo_input
        use neo_exchange
        use inter_interfaces, only: splinecof3_hi_driv, splinecof3, tf, &
          & splinecof1_hi_driv, splinecof1
        use neo_spline_data, only: lsw_linear_boozer

        implicit none

        integer :: i
        integer(I4B) :: sw1, sw2
        real(dp) :: m0, c1, cn
        real(dp), dimension(:), allocatable :: lambda
        integer(I4B) :: m, maxpos(1)
        integer, parameter :: m_max_sp = 12

        allocate (a_bmnc(ns, mnmax), b_bmnc(ns, mnmax))
        allocate (c_bmnc(ns, mnmax), d_bmnc(ns, mnmax))

        allocate (a_iota(ns), b_iota(ns))
        allocate (c_iota(ns), d_iota(ns))
        allocate (a_curr_tor(ns), b_curr_tor(ns))
        allocate (c_curr_tor(ns), d_curr_tor(ns))
        allocate (a_curr_pol(ns), b_curr_pol(ns))
        allocate (c_curr_pol(ns), d_curr_pol(ns))

        allocate (r_m(mnmax), r_mhalf(mnmax))
        allocate (sp_index(ns))

        do i = 1, mnmax
            m = ixm(i)
            if (m .le. m_max_sp) then
                r_m(i) = dble(m)
            else
                if (modulo(m, 2) .eq. 1) then
                    r_m(i) = dble(m_max_sp + 1)
                else
                    r_m(i) = dble(m_max_sp)
                end if
            end if
            r_mhalf(i) = r_m(i)/2._dp
        end do
        sp_index = (/(i, i=1, ns)/)

        ! 1-d splines of 2-d arrays
        call splinecof3_hi_driv(es, bmnc, r_mhalf, &
            & a_bmnc, b_bmnc, c_bmnc, d_bmnc, sp_index, tf)

        ! boundary types (natural spline)
        sw1 = 2
        sw2 = 4
        ! input for test function for spline
        m0 = 0.0_dp
        ! boundary condition for spline
        c1 = 0.0_dp
        cn = 0.0_dp
        ! we use no smoothing for spline
        allocate (lambda(ns))
        lambda = 1.0d0

        ! 1-d splines of 1-d arrays
        call splinecof3(es, iota, c1, cn, lambda, sp_index, sw1, sw2, &
            & a_iota, b_iota, c_iota, d_iota, m0, tf)
        call splinecof3(es, curr_tor, c1, cn, lambda, sp_index, sw1, sw2, &
            & a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor, m0, tf)
        call splinecof3(es, curr_pol, c1, cn, lambda, sp_index, sw1, sw2, &
            & a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol, m0, tf)

        deallocate (lambda)

    end subroutine neo_init_spline

    subroutine neo_read_control
        ! Read Control File
        !***********************************************************************
        ! Modules
        !***********************************************************************
        use neo_units
        use neo_control
        use neo_input
        use neo_exchange
        use sizey_bo
        use sizey_cur
        use sizey_pla
        use neo_van
        !***********************************************************************
        ! Local definitions
        !***********************************************************************
        implicit none
        character(1) :: dummy
        integer :: i, n, stat
        integer, dimension(3) :: iarr
        !***********************************************************************
        ! Open input-unit and read data
        !***********************************************************************

        theta_n = 300
        phi_n = 300
        max_m_mode = 1000
        max_n_mode = 1000
        inp_swi = 6

        close (unit=r_u1)
        ! **********************************************************************
        return

    end subroutine neo_read_control
    ! **********************************************************************

    subroutine neo_read
        ! Read Boozer Files
        !***********************************************************************
        ! Modules
        !***********************************************************************
        use nrtype
        use neo_input
        use neo_units
        use neo_control
        use neo_work
        use neo_exchange
        !***********************************************************************
        ! Local definitions
        !***********************************************************************
        implicit none

        integer :: i, j, j_m, j_n
        integer :: m, n, num_m, num_n, m_found, n_found
        integer :: mm, nn
        integer :: i_alloc
        integer :: id1, id2, id3, id4, id5, id6, id7
        integer :: extra_count
        logical :: extra_zero
        character(5) :: dummy
        character(45) :: cdum
        real(kind=dp) :: xra, xrm
        real(kind=dp) :: r_small, r_big
        !***********************************************************************
        ! Open input-unit and read first quantities
        !***********************************************************************
        open (unit=r_u1, file=in_file, status='old', form='formatted')
        !***********************************************************************

        if (inp_swi .eq. 6) then        ! Stellerator symmetry
            read (r_u1, *) dummy
            read (r_u1, *) dummy
            read (r_u1, *) dummy
            read (r_u1, *) dummy
            read (r_u1, *) dummy
            read (r_u1, *) m0b, n0b, ns, nfp, flux, r_small, r_big
            m_max = m0b + 1
            n_max = 2*n0b + 1
            mnmax = m_max*n_max

            ! **********************************************************************
            ! Allocate storage arrays
            ! **********************************************************************
            allocate (ixm(mnmax), ixn(mnmax), stat=i_alloc)
            if (i_alloc /= 0) stop 'Allocation for integer arrays failed!'

            allocate (pixm(mnmax), pixn(mnmax), stat=i_alloc)
            if (i_alloc /= 0) stop 'Allocation for integer arrays pointers failed!'

            allocate (i_m(m_max), i_n(n_max), stat=i_alloc)
            if (i_alloc /= 0) stop 'Allocation for integer arrays failed!'

          allocate (es(ns), iota(ns), curr_pol(ns), curr_tor(ns), b00(ns), stat=i_alloc)
            if (i_alloc /= 0) stop 'Allocation for real arrays failed!'

            allocate (bmnc(ns, mnmax), stat=i_alloc)
            if (i_alloc /= 0) stop 'Allocation for fourier arrays (1) failed!'
            !***********************************************************************
            ! Read input arrays
            !***********************************************************************
            do i = 1, ns
                read (r_u1, *) dummy
                read (r_u1, *) dummy
                read (r_u1, *) es(i), iota(i), curr_pol(i), curr_tor(i), &
                    dummy, dummy
                read (r_u1, *) dummy

                extra_zero = .false.
                extra_count = 0
                do j = 1, mnmax
                    if (j .gt. 1) then
                        if (ixm(j - 1) .eq. 0 .and. ixn(j - 1) .eq. 0) then
                            extra_zero = .true.
                        end if
                    end if
                    if (extra_zero) then
                        extra_count = extra_count + 1
                        if (extra_count .eq. n0b) extra_zero = .false.
                        ixm(j) = 0
                        ixn(j) = -extra_count
                        bmnc(i, j) = 0.0d0
                    else
                        read (r_u1, *) ixm(j), ixn(j), &
                            dummy, dummy, dummy, &
                            bmnc(i, j)
                    end if
                end do
            end do
        else
            write (w_us, *) 'FATAL: There is yet no other input type defined'
            stop
        end if

        ! To silence a warning maybe used uninitialized (should be false positive).
        num_n = 0

        ! Filling of i_m and i_n
        ! and pointers pixm from ixm to i_m, and pixn from ixn to i_n
        do j = 1, mnmax
            m = ixm(j)
            n = ixn(j)
            if (j .eq. 1) then
                num_m = 1
                i_m(num_m) = m
                pixm(j) = num_m
                num_n = 1
                i_n(num_n) = n
                pixn(j) = num_n
            else
                m_found = 0
                do j_m = 1, num_m
                    if (m .eq. i_m(j_m)) then
                        pixm(j) = j_m
                        m_found = 1
                    end if
                end do
                if (m_found .eq. 0) then
                    num_m = num_m + 1
                    i_m(num_m) = m
                    pixm(j) = num_m
                end if
                n_found = 0
                do j_n = 1, num_n
                    if (n .eq. i_n(j_n)) then
                        pixn(j) = j_n
                        n_found = 1
                    end if
                end do
                if (n_found .eq. 0) then
                    num_n = num_n + 1
                    i_n(num_n) = n
                    pixn(j) = num_n
                end if
            end if
        end do

        curr_pol = curr_pol*2.d-7*nfp   ! = bcovar_phi [Tm]
        curr_tor = curr_tor*2.d-7         ! = bcovar_tht [Tm]
        max_n_mode = max_n_mode*nfp
        ixn = ixn*nfp
        i_n = i_n*nfp
        ixm = ixm
        i_m = i_m
        psi_pr = flux/twopi

        close (unit=r_u1)

        return
    end subroutine neo_read
    ! **********************************************************************

    subroutine neo_prep
        ! Preparation of Arrays
        !***********************************************************************
        ! Modules
        !***********************************************************************
        use nrtype
        use neo_input
        use neo_work
        use neo_exchange
        use neo_control
        use neo_units
        use neo_spline
        !***********************************************************************
        ! Local definitions
        !***********************************************************************
        implicit none

        integer :: i_alloc
        integer :: ip, it
        integer :: im, in
        integer :: m, n
        ! **********************************************************************
        ! Allocate Storage Arrays
        ! **********************************************************************
        allocate (cosmth(theta_n, m_max), &
                  sinmth(theta_n, m_max), &
                  cosnph(phi_n, n_max), &
                  sinnph(phi_n, n_max), &
                  stat=i_alloc)
        if (i_alloc /= 0) stop 'Allocation for cos/sin-arrays failed!'
        allocate (theta_arr(theta_n), &
                  phi_arr(phi_n), &
                  stat=i_alloc)
        if (i_alloc /= 0) stop 'Allocation for theta/phi-arrays failed!'
        ! **********************************************************************
        ! Some initial work
        ! **********************************************************************
        theta_start = 0.0
        theta_end = twopi
        theta_int = (theta_end - theta_start)/(theta_n - 1)
        phi_start = 0.0
        phi_end = twopi/nfp
        phi_int = (phi_end - phi_start)/(phi_n - 1)
        ! **********************************************************************
        ! Preparation of arrays
        ! **********************************************************************
        do it = 1, theta_n
            theta_arr(it) = theta_start + theta_int*(it - 1)
        end do

        do ip = 1, phi_n
            phi_arr(ip) = phi_start + phi_int*(ip - 1)
        end do

        do im = 1, m_max
            m = i_m(im)
            if (abs(m) .le. max_m_mode) then
                do it = 1, theta_n
                    sinmth(it, im) = sin(m*theta_arr(it))
                    cosmth(it, im) = cos(m*theta_arr(it))
                end do
            end if
        end do
        do in = 1, n_max
            n = i_n(in)
            if (abs(n) .le. max_n_mode) then
                do ip = 1, phi_n
                    sinnph(ip, in) = sin(n*phi_arr(ip))
                    cosnph(ip, in) = cos(n*phi_arr(ip))
                end do
            end if
        end do

    end subroutine neo_prep

end module neo_sub_mod
