module neo_spline_over_s

    implicit none

contains
!***********************************************************************
!
! routines for calculating spline coefficients
!              drivers
!
! Author:  Bernhard Seiwald
! Date:    16.12.2000
!          05.11.2001
!
!***********************************************************************

!***********************************************************************
!
! routines for third order spline
!
!***********************************************************************

! ------  third order spline: with testfunction, LSQ, smoothing
!
! AUTHOR: Bernhard Seiwald
!
! DATE:   05.07.2001

!> compute coefs for smoothing spline
!> positions of intervals are given by indx
!>
!> if dabs(c1) > 1e30 -> c1 = 0.0D0
!> if dabs(cn) > 1e30 -> cn = 0.0D0
!>
!> INPUT:
!>     INTEGER(I4B) ,       DIMENSION(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>                                             ATTENTION:
!>                                             x(1),y(1) and x(len_x),y(len_x)
!>                                             must be gridpoints!!!
!>     REAL (kind=dp), DIMENSION(len_x) :: x ...... x values
!>     REAL (kind=dp), DIMENSION(len_x) :: y ...... y values
!>     REAL (kind=dp)                :: c1, cn .... 1. and last 2. derivative
!>     REAL (kind=dp), DIMENSION(len_indx) :: lambda . weight for 3. derivative
!>     INTEGER(I4B)                        :: sw1 ....
!>                                               = 1 -> c1 = 1. deriv 1. point
!>                                               = 2 -> c1 = 2. deriv 1. point
!>                                               = 3 -> c1 = 1. deriv N. point
!>                                               = 4 -> c1 = 2. deriv N. point
!>     INTEGER(I4B)                         :: sw2 ....
!>                                               = 1 -> cn = 1. deriv 1. point
!>                                               = 2 -> cn = 2. deriv 1. point
!>                                               = 3 -> cn = 1. deriv N. point
!>                                               = 4 -> cn = 2. deriv N. point
!>     REAL (kind=dp)                :: m ...... powers of leading term
!>
!> OUTPUT:
!>     REAL (kind=dp), DIMENSION(len_indx) :: a, b, c, d ... spline coefs
!>
!> INTERNAL:
!>     INTEGER(I4B), PARAMETER :: VAR = 7 ... no of variables
!>
!> NEEDS:
!>     calc_opt_lambda3
    subroutine splinecof3(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                          a, b, c, d, m)
        !-----------------------------------------------------------------------
        ! Modules
        !-----------------------------------------------------------------------

        use nrtype, only: I4B, DP
 !! Modifications by Andreas F. Martitsch (06.08.2014)
        !Replace standard solver from Lapack with sparse solver
        !(Bad performance for more than 1000 flux surfaces ~ (3*nsurf)^2)
        use sparse_mod, only: sparse_solve
 !! End Modifications by Andreas F. Martitsch (06.08.2014)

        !---------------------------------------------------------------------

        implicit none

        real(DP), intent(INOUT) :: c1, cn
        real(DP), dimension(:), intent(IN) :: x
        real(DP), dimension(:), intent(IN) :: y
        real(DP), dimension(:), intent(IN) :: lambda1
        integer(I4B), dimension(:), intent(IN) :: indx
        real(DP), dimension(:), intent(OUT) :: a, b, c, d
        integer(I4B), intent(IN) :: sw1, sw2
        real(DP), intent(IN) :: m

        integer(I4B), parameter :: VAR = 7
        integer(I4B) :: size_dimension
        integer(I4B) :: i_alloc, info
        integer(I4B) :: len_x, len_indx
        integer(I4B) :: i, j, l, ii, ie
        integer(I4B) :: mu1, mu2, nu1, nu2
        integer(I4B) :: sig1, sig2, rho1, rho2
        integer(I4B), dimension(:), allocatable :: indx_lu
        real(DP) :: h, h_j, x_h, help_i, help_inh
        real(DP) :: help_a, help_b, help_c, help_d
        real(DP), dimension(:, :), allocatable :: MA
        real(DP), dimension(:), allocatable :: inh, simqa, lambda, omega
        character(200) :: error_message

        len_x = size(x)
        len_indx = size(indx)
        size_dimension = VAR*len_indx - 2

        if (.not. (size(x) == size(y))) then
            write (*, *) 'splinecof3: assertion 1 failed'
            stop 'program terminated'
        end if
        if (.not. (size(a) == size(b) .and. size(a) == size(c) &
                   .and. size(a) == size(d) .and. size(a) == size(indx) &
                   .and. size(a) == size(lambda1))) then
            write (*, *) 'splinecof3: assertion 2 failed'
            stop 'program terminated'
        end if

        ! check whether points are monotonously increasing or not
        do i = 1, len_x - 1
            if (x(i) >= x(i + 1)) then
                print *, 'SPLINECOF3: error i, x(i), x(i+1)', &
                    i, x(i), x(i + 1)
                stop 'SPLINECOF3: error  wrong order of x(i)'
            end if
        end do
        ! check indx
        do i = 1, len_indx - 1
            if (indx(i) < 1) then
                print *, 'SPLINECOF3: error i, indx(i)', i, indx(i)
                stop 'SPLINECOF3: error  indx(i) < 1'
            end if
            if (indx(i) >= indx(i + 1)) then
                print *, 'SPLINECOF3: error i, indx(i), indx(i+1)', &
                    i, indx(i), indx(i + 1)
                stop 'SPLINECOF3: error  wrong order of indx(i)'
            end if
            if (indx(i) > len_x) then
                print *, 'SPLINECOF3: error i, indx(i), indx(i+1)', &
                    i, indx(i), indx(i + 1)
                stop 'SPLINECOF3: error  indx(i) > len_x'
            end if
        end do
        if (indx(len_indx) < 1) then
            print *, 'SPLINECOF3: error len_indx, indx(len_indx)', &
                len_indx, indx(len_indx)
            stop 'SPLINECOF3: error  indx(max) < 1'
        end if
        if (indx(len_indx) > len_x) then
            print *, 'SPLINECOF3: error len_indx, indx(len_indx)', &
                len_indx, indx(len_indx)
            stop 'SPLINECOF3: error  indx(max) > len_x'
        end if

        if (sw1 == sw2) then
            stop 'SPLINECOF3: error  two identical boundary conditions'
        end if

        allocate (MA(size_dimension, size_dimension), stat=i_alloc, &
                  errmsg=error_message)
        if (i_alloc /= 0) then
            write (*, *) 'splinecof3: Allocation for array ma failed '// &
                'with error message:'
            write (*, *) trim(error_message)
            write (*, *) 'size should be ', size_dimension, ' x ', size_dimension
            stop
        end if
        allocate (inh(size_dimension), indx_lu(size_dimension), stat=i_alloc, &
                  errmsg=error_message)
        if (i_alloc /= 0) then
            write (*, *) 'splinecof3: Allocation for arrays inh and indx_lu failed '// &
                'with error message:'
            write (*, *) trim(error_message)
            write (*, *) 'size should be ', size_dimension
            stop
        end if
        allocate (simqa(size_dimension*size_dimension), stat=i_alloc, &
                  errmsg=error_message)
        if (i_alloc /= 0) then
            write (*, *) 'splinecof3: Allocation for array simqa failed '// &
                'with error message:'
            write (*, *) trim(error_message)
            write (*, *) 'size should be ', size_dimension*size_dimension
            stop
        end if
        allocate (lambda(size(lambda1)), stat=i_alloc, errmsg=error_message)
        if (i_alloc /= 0) then
            write (*, *) 'splinecof3: Allocation for array lambda failed '// &
                'with error message:'
            write (*, *) trim(error_message)
            write (*, *) 'size should be ', size(lambda1)
            stop
        end if
        allocate (omega(size(lambda1)), stat=i_alloc, errmsg=error_message)
        if (i_alloc /= 0) then
            write (*, *) 'splinecof3: Allocation for array omega failed with message:'
            write (*, *) trim(error_message)
            write (*, *) 'size should be ', size(lambda1)
            stop
        end if
        !---------------------------------------------------------------------

        if (DABS(c1) > 1.0e30) then
            c1 = 0.0d0; 
        end if
        if (DABS(cn) > 1.0e30) then
            cn = 0.0d0; 
        end if

        ! setting all to zero
        MA(:, :) = 0.0d0
        inh(:) = 0.0d0

        ! calculate optimal weights for smooting (lambda)
        if (maxval(lambda1) < 0.0d0) then
            call calc_opt_lambda3(x, y, omega)
        else
            omega = lambda1
        end if
        lambda = 1.0d0 - omega

        if (sw1 == 1) then
            mu1 = 1
            nu1 = 0
            sig1 = 0
            rho1 = 0
        else if (sw1 == 2) then
            mu1 = 0
            nu1 = 1
            sig1 = 0
            rho1 = 0
        else if (sw1 == 3) then
            mu1 = 0
            nu1 = 0
            sig1 = 1
            rho1 = 0
        else if (sw1 == 4) then
            mu1 = 0
            nu1 = 0
            sig1 = 0
            rho1 = 1
        else
            stop 'SPLINECOF3: error  in using boundary condition 1'
        end if

        if (sw2 == 1) then
            mu2 = 1
            nu2 = 0
            sig2 = 0
            rho2 = 0
        else if (sw2 == 2) then
            mu2 = 0
            nu2 = 1
            sig2 = 0
            rho2 = 0
        else if (sw2 == 3) then
            mu2 = 0
            nu2 = 0
            sig2 = 1
            rho2 = 0
        else if (sw2 == 4) then
            mu2 = 0
            nu2 = 0
            sig2 = 0
            rho2 = 1
        else
            stop 'SPLINECOF3: error  in using boundary condition 2'
        end if

        ! coefs for first point
        i = 0
        j = 1
        ii = indx((j - 1)/VAR + 1)
        ie = indx((j - 1)/VAR + 2) - 1
        h = x(indx((j - 1)/VAR + 2)) - x(ii)

        ! boundary condition 1
        i = i + 1
        MA(i, 2) = dble(mu1)
        MA(i, 3) = dble(nu1)
        MA(i, (len_indx - 1)*VAR + 2) = dble(sig1)
        MA(i, (len_indx - 1)*VAR + 3) = dble(rho1)
        inh(i) = c1

        ! A_i
        i = i + 1
        MA(i, j + 0 + 0) = 1.0d0
        MA(i, j + 0 + 1) = h
        MA(i, j + 0 + 2) = h*h
        MA(i, j + 0 + 3) = h*h*h
        MA(i, j + VAR + 0) = -1.0d0
        ! B_i
        i = i + 1
        MA(i, j + 0 + 1) = 1.0d0
        MA(i, j + 0 + 2) = 2.0d0*h
        MA(i, j + 0 + 3) = 3.0d0*h*h
        MA(i, j + VAR + 1) = -1.0d0
        ! C_i
        i = i + 1
        MA(i, j + 0 + 2) = 1.0d0
        MA(i, j + 0 + 3) = 3.0d0*h
        MA(i, j + VAR + 2) = -1.0d0
        ! delta a_i
        i = i + 1
        help_a = 0.0d0
        help_b = 0.0d0
        help_c = 0.0d0
        help_d = 0.0d0
        help_i = 0.0d0
        do l = ii, ie
            h_j = x(l) - x(ii)
            x_h = monom(x(l), m)*monom(x(l), m)
            help_a = help_a + x_h
            help_b = help_b + h_j*x_h
            help_c = help_c + h_j*h_j*x_h
            help_d = help_d + h_j*h_j*h_j*x_h
            help_i = help_i + monom(x(l), m)*y(l)
        end do  ! DO l = ii, ie
        MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
        MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
        MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
        MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d
        MA(i, j + 0 + 4) = 1.0d0
        inh(i) = omega((j - 1)/VAR + 1)*help_i
        ! delta b_i
        i = i + 1
        help_a = 0.0d0
        help_b = 0.0d0
        help_c = 0.0d0
        help_d = 0.0d0
        help_i = 0.0d0
        do l = ii, ie
            h_j = x(l) - x(ii)
            x_h = monom(x(l), m)*monom(x(l), m)
            help_a = help_a + h_j*x_h
            help_b = help_b + h_j*h_j*x_h
            help_c = help_c + h_j*h_j*h_j*x_h
            help_d = help_d + h_j*h_j*h_j*h_j*x_h
            help_i = help_i + h_j*monom(x(l), m)*y(l)
        end do  ! DO l = ii, ie
        MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
        MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
        MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
        MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d
        MA(i, j + 0 + 4) = h
        MA(i, j + 0 + 5) = 1.0d0
        MA(i, (len_indx - 1)*VAR + 4) = dble(mu1)
        MA(i, (len_indx - 1)*VAR + 5) = dble(mu2)
        inh(i) = omega((j - 1)/VAR + 1)*help_i
        ! delta c_i
        i = i + 1
        help_a = 0.0d0
        help_b = 0.0d0
        help_c = 0.0d0
        help_d = 0.0d0
        help_i = 0.0d0
        do l = ii, ie
            h_j = x(l) - x(ii)
            x_h = monom(x(l), m)*monom(x(l), m)
            help_a = help_a + h_j*h_j*x_h
            help_b = help_b + h_j*h_j*h_j*x_h
            help_c = help_c + h_j*h_j*h_j*h_j*x_h
            help_d = help_d + h_j*h_j*h_j*h_j*h_j*x_h
            help_i = help_i + h_j*h_j*monom(x(l), m)*y(l)
        end do  ! DO l = ii, ie
        MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
        MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
        MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
        MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d
        MA(i, j + 0 + 4) = h*h
        MA(i, j + 0 + 5) = 2.0d0*h
        MA(i, j + 0 + 6) = 1.0d0
        MA(i, (len_indx - 1)*VAR + 4) = dble(nu1)
        MA(i, (len_indx - 1)*VAR + 5) = dble(nu2)
        inh(i) = omega((j - 1)/VAR + 1)*help_i
        ! delta DELTA d_i
        i = i + 1
        help_a = 0.0d0
        help_b = 0.0d0
        help_c = 0.0d0
        help_d = 0.0d0
        help_i = 0.0d0
        do l = ii, ie
            h_j = x(l) - x(ii)
            x_h = monom(x(l), m)*monom(x(l), m)
            help_a = help_a + h_j*h_j*h_j*x_h
            help_b = help_b + h_j*h_j*h_j*h_j*x_h
            help_c = help_c + h_j*h_j*h_j*h_j*h_j*x_h
            help_d = help_d + h_j*h_j*h_j*h_j*h_j*h_j*x_h
            help_i = help_i + h_j*h_j*h_j*monom(x(l), m)*y(l)
        end do  ! DO l = ii, ie
        MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
        MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
        MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
        MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d + lambda((j - 1)/VAR + 1)
        MA(i, j + 0 + 4) = h*h*h
        MA(i, j + 0 + 5) = 3.0d0*h*h
        MA(i, j + 0 + 6) = 3.0d0*h
        inh(i) = omega((j - 1)/VAR + 1)*help_i

        ! coefs for point 2 to len_x_points-1
        do j = VAR + 1, VAR*(len_indx - 1) - 1, VAR
            ii = indx((j - 1)/VAR + 1)
            ie = indx((j - 1)/VAR + 2) - 1
            h = x(indx((j - 1)/VAR + 2)) - x(ii)
            ! A_i
            i = i + 1
            MA(i, j + 0 + 0) = 1.0d0
            MA(i, j + 0 + 1) = h
            MA(i, j + 0 + 2) = h*h
            MA(i, j + 0 + 3) = h*h*h
            MA(i, j + VAR + 0) = -1.0d0
            ! B_i
            i = i + 1
            MA(i, j + 0 + 1) = 1.0d0
            MA(i, j + 0 + 2) = 2.0d0*h
            MA(i, j + 0 + 3) = 3.0d0*h*h
            MA(i, j + VAR + 1) = -1.0d0
            ! C_i
            i = i + 1
            MA(i, j + 0 + 2) = 1.0d0
            MA(i, j + 0 + 3) = 3.0d0*h
            MA(i, j + VAR + 2) = -1.0d0
            ! delta a_i
            i = i + 1
            help_a = 0.0d0
            help_b = 0.0d0
            help_c = 0.0d0
            help_d = 0.0d0
            help_i = 0.0d0
            do l = ii, ie
                h_j = x(l) - x(ii)
                x_h = monom(x(l), m)*monom(x(l), m)
                help_a = help_a + x_h
                help_b = help_b + h_j*x_h
                help_c = help_c + h_j*h_j*x_h
                help_d = help_d + h_j*h_j*h_j*x_h
                help_i = help_i + monom(x(l), m)*y(l)
            end do   ! DO l = ii, ie
            MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
            MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
            MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
            MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d
            MA(i, j + 0 + 4) = 1.0d0
            MA(i, j - VAR + 4) = -1.0d0
            inh(i) = omega((j - 1)/VAR + 1)*help_i
            ! delta b_i
            i = i + 1
            help_a = 0.0d0
            help_b = 0.0d0
            help_c = 0.0d0
            help_d = 0.0d0
            help_i = 0.0d0
            do l = ii, ie
                h_j = x(l) - x(ii)
                x_h = monom(x(l), m)*monom(x(l), m)
                help_a = help_a + h_j*x_h
                help_b = help_b + h_j*h_j*x_h
                help_c = help_c + h_j*h_j*h_j*x_h
                help_d = help_d + h_j*h_j*h_j*h_j*x_h
                help_i = help_i + h_j*monom(x(l), m)*y(l)
            end do  ! DO l = ii, ie
            MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
            MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
            MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
            MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d
            MA(i, j + 0 + 4) = h
            MA(i, j + 0 + 5) = 1.0d0
            MA(i, j - VAR + 5) = -1.0d0
            inh(i) = omega((j - 1)/VAR + 1)*help_i
            ! delta c_i
            i = i + 1
            help_a = 0.0d0
            help_b = 0.0d0
            help_c = 0.0d0
            help_d = 0.0d0
            help_i = 0.0d0
            do l = ii, ie
                h_j = x(l) - x(ii)
                x_h = monom(x(l), m)*monom(x(l), m)
                help_a = help_a + h_j*h_j*x_h
                help_b = help_b + h_j*h_j*h_j*x_h
                help_c = help_c + h_j*h_j*h_j*h_j*x_h
                help_d = help_d + h_j*h_j*h_j*h_j*h_j*x_h
                help_i = help_i + h_j*h_j*monom(x(l), m)*y(l)
            end do  ! DO l = ii, ie
            MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
            MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
            MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
            MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d
            MA(i, j + 0 + 4) = h*h
            MA(i, j + 0 + 5) = 2.0d0*h
            MA(i, j + 0 + 6) = 1.0d0
            MA(i, j - VAR + 6) = -1.0d0
            inh(i) = omega((j - 1)/VAR + 1)*help_i
            ! delta DELTA d_i
            i = i + 1
            help_a = 0.0d0
            help_b = 0.0d0
            help_c = 0.0d0
            help_d = 0.0d0
            help_i = 0.0d0
            do l = ii, ie
                h_j = x(l) - x(ii)
                x_h = monom(x(l), m)*monom(x(l), m)
                help_a = help_a + h_j*h_j*h_j*x_h
                help_b = help_b + h_j*h_j*h_j*h_j*x_h
                help_c = help_c + h_j*h_j*h_j*h_j*h_j*x_h
                help_d = help_d + h_j*h_j*h_j*h_j*h_j*h_j*x_h
                help_i = help_i + h_j*h_j*h_j*monom(x(l), m)*y(l)
            end do  ! DO l = ii, ie
            MA(i, j + 0 + 0) = omega((j - 1)/VAR + 1)*help_a
            MA(i, j + 0 + 1) = omega((j - 1)/VAR + 1)*help_b
            MA(i, j + 0 + 2) = omega((j - 1)/VAR + 1)*help_c
            MA(i, j + 0 + 3) = omega((j - 1)/VAR + 1)*help_d + lambda((j - 1)/VAR + 1)
            MA(i, j + 0 + 4) = h*h*h
            MA(i, j + 0 + 5) = 3.0d0*h*h
            MA(i, j + 0 + 6) = 3.0d0*h
            inh(i) = omega((j - 1)/VAR + 1)*help_i
        end do  ! DO j = VAR+1, VAR*(len_indx-1)-1, VAR

        ! last point
        ! delta a_i
        i = i + 1
        ii = indx((j - 1)/VAR + 1)
        ie = ii
        help_a = 0.0d0
        help_inh = 0.0d0
        l = ii
        help_a = help_a + monom(x(l), m)*monom(x(l), m)
        help_inh = help_inh + monom(x(l), m)*y(l)

        MA(i, (len_indx - 1)*VAR + 1) = omega((j - 1)/VAR + 1)*help_a
        MA(i, (len_indx - 2)*VAR + 5) = omega((j - 1)/VAR + 1)*(-1.0d0)
        inh(i) = omega((j - 1)/VAR + 1)*help_inh
        ! delta b_i
        i = i + 1
        MA(i, (len_indx - 2)*VAR + 6) = -1.0d0
        MA(i, (len_indx - 1)*VAR + 4) = dble(sig1)
        MA(i, (len_indx - 1)*VAR + 5) = dble(sig2)
        ! delta c_i
        i = i + 1
        MA(i, (len_indx - 2)*VAR + 7) = -1.0d0
        MA(i, (len_indx - 1)*VAR + 4) = dble(rho1)
        MA(i, (len_indx - 1)*VAR + 5) = dble(rho2)

        ! boundary condition 2
        i = i + 1
        MA(i, 2) = dble(mu2)
        MA(i, 3) = dble(nu2)
        MA(i, (len_indx - 1)*VAR + 2) = dble(sig2)
        MA(i, (len_indx - 1)*VAR + 3) = dble(rho2)
        inh(i) = cn

! ---------------------------

        ! solve system
        call sparse_solve(MA, inh)

        ! take a(), b(), c(), d()
        do i = 1, len_indx
            a(i) = inh((i - 1)*VAR + 1)
            b(i) = inh((i - 1)*VAR + 2)
            c(i) = inh((i - 1)*VAR + 3)
            d(i) = inh((i - 1)*VAR + 4)
        end do

        deallocate (MA, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3: Deallocation for arrays 1 failed!'
        deallocate (inh, indx_lu, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3: Deallocation for arrays 2 failed!'
        deallocate (simqa, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3: Deallocation for arrays 3 failed!'
        deallocate (lambda, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3: Deallocation for lambda failed!'
        deallocate (omega, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3: Deallocation for omega failed!'

    end subroutine splinecof3

!> reconstruct spline coefficients (a, b, c, d) on x(i)
!>
!> h := (x - x_i)
!>
!> INPUT:
!>  REAL(DP)                :: ai, bi, ci, di ... old coefs
!>  REAL(DP)                :: h ................ h := x(i) - x(i-1)
!>
!> OUTPUT:
!>  REAL(DP)                :: a, b, c, d ....... new coefs
    subroutine reconstruction3(ai, bi, ci, di, h, a, b, c, d)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------

        use nrtype, only: DP

        !---------------------------------------------------------------------

        implicit none

        real(DP), intent(IN) :: ai, bi, ci, di
        real(DP), intent(IN) :: h
        real(DP), intent(OUT) :: a, b, c, d

        !---------------------------------------------------------------------

        d = di
        c = ci + 3.0d0*h*di
        b = bi + h*(2.0d0*ci + 3.0d0*h*di)
        a = ai + h*(bi + h*(ci + h*di))

    end subroutine reconstruction3

!> driver routine for splinecof3 ; used for Rmn, Zmn
!>
!> INPUT:
!>     INTEGER(I4B), DIMENSION(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>     REAL(DP),     DIMENSION(no) :: x ...... x values
!>     REAL(DP),     DIMENSION(no) :: y ...... y values
!>     REAL(DP)                    :: c1, cn . 1. and last 2. derivative
!>     REAL(DP),     DIMENSION(ns) :: lambda . weight for 3. derivative
!>     INTEGER(I4B), DIMENSION(ns) :: w ...... weight for point (0,1)
!>     INTEGER(I4B)                :: sw1 .... = 1 -> c1 = 1. deriv 1. point
!>                                             = 2 -> c1 = 2. deriv 1. point
!>                                             = 3 -> c1 = 1. deriv N. point
!>                                             = 4 -> c1 = 2. deriv N. point
!>     INTEGER(I4B)                :: sw2 .... = 1 -> cn = 1. deriv 1. point
!>                                             = 2 -> cn = 2. deriv 1. point
!>                                             = 3 -> cn = 1. deriv N. point
!>                                             = 4 -> cn = 2. deriv N. point
!>     REAL(DP)                :: m ...... powers of leading term
!>
!> OUTPUT:
!>     REAL(DP), DIMENSION(ns) :: a ...... spline coefs
!>     REAL(DP), DIMENSION(ns) :: b ...... spline coefs
!>     REAL(DP), DIMENSION(ns) :: c ...... spline coefs
!>     REAL(DP), DIMENSION(ns) :: d ...... spline coefs
!>
!> INTERNAL:
!>     INTEGER(I4B), PARAMETER :: VAR = 7 ... no of variables
    subroutine splinecof3_lo_driv(x, y, c1, cn, lambda, w, indx, &
                                  sw1, sw2, a, b, c, d, m)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------

        use nrtype, only: I4B, DP

        !---------------------------------------------------------------------

        implicit none

        integer(I4B), dimension(:), intent(IN) :: indx
        real(DP), intent(IN) :: m
        real(DP), intent(INOUT) :: c1, cn
        real(DP), dimension(:), intent(IN) :: x
        real(DP), dimension(:), intent(IN) :: y
        real(DP), dimension(:), intent(IN) :: lambda
        integer(I4B), dimension(:), intent(IN) :: w
        real(DP), dimension(:), intent(OUT) :: a, b, c, d
        integer(I4B), intent(IN) :: sw1, sw2

        integer(I4B) :: dim, no, ns, len_indx
        integer(I4B) :: i, j, ie, i_alloc
        integer(I4B) :: shift, shifti, shiftv
        integer(I4B), dimension(:), allocatable :: hi, indx1
        real(DP) :: h
        real(DP), dimension(:), allocatable :: xn, yn, lambda1
        real(DP), dimension(:), allocatable :: ai, bi, ci, di

        no = size(x)
        ns = size(a)
        len_indx = size(indx)

        !---------------------------------------------------------------------

        dim = sum(w)

        if (dim == 0) then
            stop 'error in splinecof3_lo_driv: w == 0'
        end if

        allocate (ai(dim), bi(dim), ci(dim), di(dim), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_lo_driv: allocation for arrays 1 failed!'
        allocate (indx1(dim), lambda1(dim), hi(no), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_lo_driv: allocation for arrays 2 failed!'

        hi = 1
        do i = 1, size(w)
            if ((w(i) /= 0) .and. (w(i) /= 1)) then
                stop 'splinecof3_lo_driv: wrong value for w  (0/1)'
            end if
            if (w(i) == 0) then
                if ((i + 1) <= size(w)) then
                    ie = indx(i + 1) - 1
                else
                    ie = size(hi)
                end if
                do j = indx(i), ie
                    hi(j) = 0
                end do
            end if
        end do

        dim = sum(hi)
        allocate (xn(dim), yn(dim), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_lo_driv: allocation for arrays 3 failed!'

        ! create new vectors for indx and lambda with respect to skipped points
        j = 1
        shifti = 0
        shiftv = 0
        do i = 1, size(indx)
            if (j <= size(indx1)) then
                indx1(j) = indx(i) - shiftv
                lambda1(j) = lambda(i - shifti)
            end if
            if (w(i) /= 0) then
                j = j + 1
            else
                shifti = shifti + 1
                if (i + 1 <= size(indx)) then
                    shiftv = shiftv + indx(i + 1) - indx(i)
                end if
            end if
        end do

        ! create new vectors for x and y with respect to skipped points
        j = indx1(1)
        do i = 1, size(hi)
            if (hi(i) /= 0) then
                xn(j) = x(i)
                yn(j) = y(i)
                j = j + 1
            end if
        end do

        call splinecof3(xn, yn, c1, cn, lambda1, indx1, sw1, sw2, &
                        ai, bi, ci, di, m)

        ! find first regular point
        shift = 1
        do while ((shift <= size(w)) .and. (w(shift) == 0))
            shift = shift + 1
        end do

        ! reconstruct spline coefficients from 0 to first calculated coeff.
        if ((shift > 1) .and. (shift < size(w))) then
            a(shift) = ai(1)
            b(shift) = bi(1)
            c(shift) = ci(1)
            d(shift) = di(1)
            do i = shift - 1, 1, -1
                h = x(indx(i)) - x(indx(i + 1))
                call reconstruction3(a(i + 1), b(i + 1), c(i + 1), d(i + 1), h, &
                                     a(i), b(i), c(i), d(i))
            end do
        end if

        ! reconstruct all other spline coefficients if needed
        j = 0
        do i = shift, ns
            if (w(i) == 1) then
                j = j + 1
                a(i) = ai(j)
                b(i) = bi(j)
                c(i) = ci(j)
                d(i) = di(j)
            else
                h = x(indx(i)) - x(indx(i - 1))
                call reconstruction3(a(i - 1), b(i - 1), c(i - 1), d(i - 1), h, &
                                     a(i), b(i), c(i), d(i))
            end if
        end do

        deallocate (ai, bi, ci, di, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_lo_driv: Deallocation for arrays 1 failed!'
        deallocate (indx1, lambda1, hi, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_lo_driv: Deallocation for arrays 2 failed!'
        deallocate (xn, yn, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_lo_driv: Deallocation for arrays 3 failed!'

    end subroutine splinecof3_lo_driv

!> driver routine for splinecof3_lo_driv
!>
!> INPUT:
!>     INTEGER(I4B) , DIMENSION(len_indx)  :: indx ... index vector
!>                                            contains index of grid points
!>     INTEGER(I4B),                       :: choose_rz  1: calc Rmn; 2: Zmn
!>     REAL(DP), DIMENSION(no)        :: x ...... x values
!>     REAL(DP), DIMENSION(no,no_cur) :: y ...... y values
!>     REAL(DP), DIMENSION(no_cur)    :: m ...... powers of leading term
!>
!> OUTPUT:
!>     REAL(DP), DIMENSION(ns,no_cur) :: a ...... spline coefs
!>     REAL(DP), DIMENSION(ns,no_cur) :: b ...... spline coefs
!>     REAL(DP), DIMENSION(ns,no_cur) :: c ...... spline coefs
!>     REAL(DP), DIMENSION(ns,no_cur) :: d ...... spline coefs
!> INTERNAL:
!>     REAL(DP),     DIMENSION(ns,no_cur) :: lambda3 . weight for 3. derivative
!>     INTEGER(I4B), DIMENSION(ns,no_cur) :: w ....... weight for point (0,1)
    subroutine splinecof3_hi_driv(x, y, m, a, b, c, d, indx)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------

        use nrtype, only: I4B, DP

        !---------------------------------------------------------------------

        implicit none

        integer(I4B), dimension(:), intent(IN) :: indx
        real(DP), dimension(:), intent(IN) :: m
        real(DP), dimension(:), intent(IN) :: x
        real(DP), dimension(:, :), intent(IN) :: y
        real(DP), dimension(:, :), intent(OUT) :: a, b, c, d

        real(DP), dimension(:, :), allocatable :: lambda3
        integer(I4B), dimension(:, :), allocatable :: w
        integer(I4B) :: ns, no_cur
        integer(I4B) :: i, sw1, sw2, i_alloc
        real(DP) :: c1, cn

        !---------------------------------------------------------------------

        ns = size(a, 1)
        no_cur = size(y, 2)

        allocate (lambda3(ns, size(y, 2)), w(ns, size(y, 2)), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_hi_driv: Allocation for arrays failed!'

        ! lambda3 = -1.0D0   !! automatic smoothing
        lambda3 = 1.0d0     !! no smoothing

        ! weights:  w(i)=0/1;  if(w(i)==0) ... do not use this point
        w = 1

        sw1 = 2
        sw2 = 4

        c1 = 0.0d0
        cn = 0.0d0

        do i = 1, no_cur
            if (m(i) /= 0.0d0) then
                w(1, i) = 0   ! system is not defined at y(0)=0
            end if
            call splinecof3_lo_driv(x, y(:, i), c1, cn, &
                                    lambda3(:, i), w(:, i), indx, sw1, sw2, &
                                    a(:, i), b(:, i), c(:, i), d(:, i), m(i))
        end do

        deallocate (lambda3, w, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof3_hi_driv: Deallocation for arrays failed!'

    end subroutine splinecof3_hi_driv

!> calculate optimal weights for smooting (lambda)
!>
!> \attention  NO FINAL VERSION NOW!!!!!
    subroutine calc_opt_lambda3(x, y, lambda)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------

        use nrtype, only: I4B, DP
        !---------------------------------------------------------------------

        implicit none

        real(DP), dimension(:), intent(IN) :: x, y
        real(DP), dimension(:), intent(OUT) :: lambda

        integer(I4B) :: i, no
        real(DP) :: av_a
        real(DP) :: ymax, xd(3), yd(3)

        !---------------------------------------------------------------------

        no = size(x)
        av_a = 0.0d0
        ymax = maxval(abs(y))
        if (ymax == 0.0d0) ymax = 1.0d0

        do i = 1, no
            if (i == 1) then
                xd(1) = x(2)
                xd(2) = x(1)
                xd(3) = x(3)
                yd(1) = y(2)
                yd(2) = y(1)
                yd(3) = y(3)
                call dist_lin(xd, yd, ymax, av_a)
            else if (i == no) then
                xd(1) = x(no - 2)
                xd(2) = x(no)
                xd(3) = x(no - 1)
                yd(1) = y(no - 2)
                yd(2) = y(no)
                yd(3) = y(no - 1)
                call dist_lin(xd, yd, ymax, av_a)
            else
                call dist_lin(x(i - 1:i + 1), y(i - 1:i + 1), ymax, av_a)
            end if
            lambda(i) = 1.0d0 - av_a**3
        end do
        av_a = sum(lambda)/dble(size(lambda))

        lambda = av_a
        lambda(1) = 1.0d0
        lambda(no) = 1.0d0

    end subroutine calc_opt_lambda3

    subroutine dist_lin(x, y, ymax, dist)

        use nrtype, only: DP

        implicit none

        real(DP), dimension(:), intent(IN) :: x, y
        real(DP), intent(IN) :: ymax
        real(DP), intent(OUT) :: dist

        real(DP) :: k, d
        ! --------------------------------------------------------------------

        k = (y(3) - y(1))/(x(3) - x(1))
        d = (y(1)*x(3) - y(3)*x(1))/(x(3) - x(1))

        dist = abs((y(2) - (k*x(2) + d))/ymax)

    end subroutine dist_lin

! ------  first order spline (linear interpolation)

!> compute coefs for smoothing spline
!> positions of intervals are given by indx
!>
!> if dabs(c1) > 1e30 -> c1 = 0.0D0
!> if dabs(cn) > 1e30 -> cn = 0.0D0
!>
!> INPUT:
!>     integer(I4B),   dimension(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>                                             ATTENTION:
!>                                             x(1),y(1) and x(len_x),y(len_x)
!>                                             must be gridpoints!!!
!>     real (kind=dp), dimension(len_x) :: x ...... x values
!>     real (kind=dp), dimension(len_x) :: y ...... y values
!>     real (kind=dp)                :: c1, cn .... ignored
!>     real (kind=dp), dimension(len_indx) :: lambda ignored
!>     integer(I4B)                        :: sw1 ignored
!>     integer(I4B)                         :: sw2 ignored
!>     real (kind=dp)                :: m ...... ignored
!>
!> OUTPUT:
!>     real (kind=dp), dimension(len_indx) :: a, b, c, d ... spline coefs
    subroutine splinecof1(x, y, c1, cn, lambda1, indx, sw1, sw2, &
       & a, b, c, d, m)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------
        use nrtype, only: I4B, DP

        implicit none

        real(DP), intent(inout) :: c1, cn
        real(DP), dimension(:), intent(in) :: x
        real(DP), dimension(:), intent(in) :: y
        real(DP), dimension(:), intent(in) :: lambda1
        integer(I4B), dimension(:), intent(in) :: indx
        real(DP), dimension(:), intent(out) :: a, b, c, d
        integer(I4B), intent(in) :: sw1, sw2
        real(DP), intent(in) :: m

        integer(I4B) :: len_x, len_indx
        integer(I4B) :: i

        len_x = size(x)
        len_indx = size(indx)

        if (.not. (size(x) == size(y))) then
            write (*, *) 'splinecof1: assertion 1 failed'
            stop 'program terminated'
        end if
        if (.not. (size(a) == size(b) .and. size(a) == size(c) &
                   .and. size(a) == size(d) .and. size(a) == size(indx) &
                   .and. size(a) == size(lambda1))) then
            write (*, *) 'splinecof1: assertion 2 failed'
            stop 'program terminated'
        end if

        ! check whether points are monotonously increasing or not
        do i = 1, len_x - 1
            if (x(i) >= x(i + 1)) then
                print *, 'SPLINECOF1: error i, x(i), x(i+1)', &
                    i, x(i), x(i + 1)
                stop 'SPLINECOF1: error  wrong order of x(i)'
            end if
        end do
        ! check indx
        do i = 1, len_indx - 1
            if (indx(i) < 1) then
                print *, 'SPLINECOF1: error i, indx(i)', i, indx(i)
                stop 'SPLINECOF1: error  indx(i) < 1'
            end if
            if (indx(i) >= indx(i + 1)) then
                print *, 'SPLINECOF1: error i, indx(i), indx(i+1)', &
                    i, indx(i), indx(i + 1)
                stop 'SPLINECOF1: error  wrong order of indx(i)'
            end if
            if (indx(i) > len_x) then
                print *, 'SPLINECOF1: error i, indx(i), indx(i+1)', &
                    i, indx(i), indx(i + 1)
                stop 'SPLINECOF1: error  indx(i) > len_x'
            end if
        end do
        if (indx(len_indx) < 1) then
            print *, 'SPLINECOF1: error len_indx, indx(len_indx)', &
                len_indx, indx(len_indx)
            stop 'SPLINECOF3: error  indx(max) < 1'
        end if
        if (indx(len_indx) > len_x) then
            print *, 'SPLINECOF1: error len_indx, indx(len_indx)', &
                len_indx, indx(len_indx)
            stop 'SPLINECOF1: error  indx(max) > len_x'
        end if

        if (sw1 == sw2) then
            stop 'SPLINECOF1: error  two identical boundary conditions'
        end if

        if (dabs(c1) > 1.0e30) then
            c1 = 0.0d0; 
        end if
        if (dabs(cn) > 1.0e30) then
            cn = 0.0d0; 
        end if

        ! ---------------------------

        do i = 1, len_indx - 1
            b(i) = (y(i + 1) - y(i))/(x(i + 1) - x(i))
            a(i) = y(i) ! - b(i) * x(i) ! this term cancels, because we assume coordinate system is centered at x(i), and thus x(i) = 0.
        end do

        a(len_indx) = a(len_indx - 1)
        b(len_indx) = b(len_indx - 1)

        c = 0.0
        d = 0.0

    end subroutine splinecof1

!> reconstruct spline coefficients (a, b, c, d) on x(i)
!>
!> h := (x - x_i)
!>
!> INPUT:
!>  rela(DP)                :: ai, bi, ci, di ... old coefs
!>  real(DP)                :: h ................ h := x(i) - x(i-1)
!>
!> OUTPUT:
!>  real(DP)                :: a, b, c, d ....... new coefs
    subroutine reconstruction1(ai, bi, ci, di, h, a, b, c, d)
        !-----------------------------------------------------------------------
        ! Modules
        !-----------------------------------------------------------------------
        use nrtype, only: DP

        implicit none

        real(DP), intent(in) :: ai, bi, ci, di
        real(DP), intent(in) :: h
        real(DP), intent(out) :: a, b, c, d

        d = 0.0
        c = 0.0
        b = bi
        a = ai + h*bi

    end subroutine reconstruction1

!> driver routine for splinecof1 ; used for Rmn, Zmn
!>
!> INPUT:
!>     integer(I4B), dimension(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>     real(DP),     dimension(no) :: x ...... x values
!>     real(DP),     dimension(no) :: y ...... y values
!>     real(DP)                    :: c1, cn . 1. and last 2. derivative
!>     real(DP),     dimension(ns) :: lambda . weight for 3. derivative
!>     integer(I4B), dimension(ns) :: w ...... weight for point (0,1)
!>     integer(I4B)                :: sw1 .... = 1 -> c1 = 1. deriv 1. point
!>                                             = 2 -> c1 = 2. deriv 1. point
!>                                             = 3 -> c1 = 1. deriv N. point
!>                                             = 4 -> c1 = 2. deriv N. point
!>     integer(I4B)                :: sw2 .... = 1 -> cn = 1. deriv 1. point
!>                                             = 2 -> cn = 2. deriv 1. point
!>                                             = 3 -> cn = 1. deriv N. point
!>                                             = 4 -> cn = 2. deriv N. point
!>     real(DP)                :: m ...... powers of leading term
!>
!> OUTPUT:
!>     real(DP), dimension(ns) :: a ...... spline coefs
!>     real(DP), dimension(ns) :: b ...... spline coefs
!>     real(DP), dimension(ns) :: c ...... spline coefs
!>     real(DP), dimension(ns) :: d ...... spline coefs
!>
!> INTERNAL:
!>     integer(I4B), parameter :: VAR = 7 ... no of variables
    subroutine splinecof1_lo_driv(x, y, c1, cn, lambda, w, indx, &
       & sw1, sw2, a, b, c, d, m)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------
        use nrtype, only: I4B, DP

        !-----------------------------------------------------------------------
        implicit none

        integer(I4B), dimension(:), intent(in) :: indx
        real(DP), intent(in) :: m
        real(DP), intent(inout) :: c1, cn
        real(DP), dimension(:), intent(in) :: x
        real(DP), dimension(:), intent(in) :: y
        real(DP), dimension(:), intent(in) :: lambda
        integer(I4B), dimension(:), intent(in) :: w
        real(DP), dimension(:), intent(out) :: a, b, c, d
        integer(I4B), intent(in) :: sw1, sw2

        integer(I4B) :: dim, no, ns, len_indx
        integer(I4B) :: i, j, ie, i_alloc
        integer(I4B) :: shift, shifti, shiftv
        integer(I4B), dimension(:), allocatable :: hi, indx1
        real(DP) :: h
        real(DP), dimension(:), allocatable :: xn, yn, lambda1
        real(DP), dimension(:), allocatable :: ai, bi, ci, di

        no = size(x)
        ns = size(a)
        len_indx = size(indx)

        !---------------------------------------------------------------------

        dim = sum(w)

        if (dim == 0) then
            stop 'error in splinecof1_lo_driv: w == 0'
        end if

        allocate (ai(dim), bi(dim), ci(dim), di(dim), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_lo_driv: allocation for arrays 1 failed!'
        allocate (indx1(dim), lambda1(dim), hi(no), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_lo_driv: allocation for arrays 2 failed!'

        hi = 1
        do i = 1, size(w)
            if ((w(i) /= 0) .and. (w(i) /= 1)) then
                stop 'splinecof1_lo_driv: wrong value for w  (0/1)'
            end if
            if (w(i) == 0) then
                if ((i + 1) <= size(w)) then
                    ie = indx(i + 1) - 1
                else
                    ie = size(hi)
                end if
                do j = indx(i), ie
                    hi(j) = 0
                end do
            end if
        end do

        dim = sum(hi)
        allocate (xn(dim), yn(dim), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_lo_driv: allocation for arrays 3 failed!'

        ! create new vectors for indx and lambda with respect to skipped points
        j = 1
        shifti = 0
        shiftv = 0
        do i = 1, size(indx)
            if (j <= size(indx1)) then
                indx1(j) = indx(i) - shiftv
                lambda1(j) = lambda(i - shifti)
            end if
            if (w(i) /= 0) then
                j = j + 1
            else
                shifti = shifti + 1
                if (i + 1 <= size(indx)) then
                    shiftv = shiftv + indx(i + 1) - indx(i)
                end if
            end if
        end do

        ! create new vectors for x and y with respect to skipped points
        j = indx1(1)
        do i = 1, size(hi)
            if (hi(i) /= 0) then
                xn(j) = x(i)
                yn(j) = y(i)
                j = j + 1
            end if
        end do

        call splinecof1(xn, yn, c1, cn, lambda1, indx1, sw1, sw2, &
            & ai, bi, ci, di, m)

        ! find first regular point
        shift = 1
        do while ((shift <= size(w)) .and. (w(shift) == 0))
            shift = shift + 1
        end do

        ! reconstruct spline coefficients from 0 to first calculated coeff.
        if ((shift > 1) .and. (shift < size(w))) then
            a(shift) = ai(1)
            b(shift) = bi(1)
            c(shift) = ci(1)
            d(shift) = di(1)
            do i = shift - 1, 1, -1
                h = x(indx(i)) - x(indx(i + 1))
                call reconstruction1(a(i + 1), b(i + 1), c(i + 1), d(i + 1), h, &
                    & a(i), b(i), c(i), d(i))
            end do
        end if

        ! reconstruct all other spline coefficients if needed
        j = 0
        do i = shift, ns
            if (w(i) == 1) then
                j = j + 1
                a(i) = ai(j)
                b(i) = bi(j)
                c(i) = ci(j)
                d(i) = di(j)
            else
                h = x(indx(i)) - x(indx(i - 1))
                call reconstruction1(a(i - 1), b(i - 1), c(i - 1), d(i - 1), h, &
                    & a(i), b(i), c(i), d(i))
            end if
        end do

        deallocate (ai, bi, ci, di, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_lo_driv: Deallocation for arrays 1 failed!'
        deallocate (indx1, lambda1, hi, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_lo_driv: Deallocation for arrays 2 failed!'
        deallocate (xn, yn, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_lo_driv: Deallocation for arrays 3 failed!'

    end subroutine splinecof1_lo_driv

!> driver routine for splinecof1_lo_driv
!>
!> INPUT:
!>     integer(I4B) , dimension(len_indx)  :: indx ... index vector
!>                                            contains index of grid points
!>     integer(I4B),                       :: choose_rz  1: calc Rmn; 2: Zmn
!>     real(DP), dimension(no)        :: x ...... x values
!>     real(DP), dimension(no,no_cur) :: y ...... y values
!>     real(DP), dimension(no_cur)    :: m ...... powers of leading term
!>
!> OUTPUT:
!>     real(DP), dimension(ns,no_cur) :: a ...... spline coefs
!>     real(DP), dimension(ns,no_cur) :: b ...... spline coefs
!>     real(DP), dimension(ns,no_cur) :: c ...... spline coefs
!>     real(DP), dimension(ns,no_cur) :: d ...... spline coefs
!> INTERNAL:
!>     real(DP),     dimension(ns,no_cur) :: lambda3 . weight for 3. derivative
!>     integer(I4B), dimension(ns,no_cur) :: w ....... weight for point (0,1)
    subroutine splinecof1_hi_driv(x, y, m, a, b, c, d, indx)
        !---------------------------------------------------------------------
        ! Modules
        !---------------------------------------------------------------------
        use nrtype, only: I4B, DP

        !---------------------------------------------------------------------

        implicit none

        integer(I4B), dimension(:), intent(in) :: indx
        real(DP), dimension(:), intent(in) :: m
        real(DP), dimension(:), intent(in) :: x
        real(DP), dimension(:, :), intent(in) :: y
        real(DP), dimension(:, :), intent(out) :: a, b, c, d

        real(DP), dimension(:, :), allocatable :: lambda3
        integer(I4B), dimension(:, :), allocatable :: w
        integer(I4B) :: ns, no_cur
        integer(I4B) :: i, sw1, sw2, i_alloc
        real(DP) :: c1, cn

        !---------------------------------------------------------------------

        ns = size(a, 1)
        no_cur = size(y, 2)

        allocate (lambda3(ns, size(y, 2)), w(ns, size(y, 2)), stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_hi_driv: Allocation for arrays failed!'

        lambda3 = 1.0d0     !! no smoothing

        ! weights:  w(i)=0/1;  if (w(i)==0) ... do not use this point
        w = 1

        sw1 = 2
        sw2 = 4

        c1 = 0.0d0
        cn = 0.0d0

        do i = 1, no_cur
            if (m(i) /= 0.0d0) then
                w(1, i) = 0   ! system is not defined at y(0)=0
            end if
            call splinecof1_lo_driv(x, y(:, i), c1, cn, &
                & lambda3(:, i), w(:, i), indx, sw1, sw2,&
                & a(:, i), b(:, i), c(:, i), d(:, i), m(i))
        end do

        deallocate (lambda3, w, stat=i_alloc)
        if (i_alloc /= 0) stop 'splinecof1_hi_driv: Deallocation for arrays failed!'

    end subroutine splinecof1_hi_driv

    subroutine spline_1d(xa, a, b, c, d, swd, order, x_in, y, yp, ypp, yppp)
        ! Computes value y(x) dy/dx(x) from cubic spline
        !
        ! Attention - fastest routine; no check at all
        !
        ! Input:
        !         xa(n)         x-values
        !         a(n),b(n),c(n),d(n)  coefs from spline
        !         swd           Switch for derivatives (0: no / 1: yes)
        !         order         powers of leading term
        !         x_in          x-value for y(x_in) and yp(x_in)
        ! Output:
        !         y             y-value at x_in
        !         yp            dy/dx-value at x_in
        !         ypp           d2y/dx2-value at x_in
        !         yppp          d3y/dx3-value at x_in
        !
        !-----------------------------------------------------------------------
        ! Modules
        !-----------------------------------------------------------------------

        use nrtype, only: I4B, DP, splinecof_compatibility

        !-----------------------------------------------------------------------

        implicit none

        integer(I4B), intent(IN) :: swd
        real(DP), intent(IN) :: order
        real(DP), dimension(:), intent(IN) :: xa, a, b, c, d
        real(DP), intent(IN) :: x_in
        real(DP), intent(OUT) :: y, yp, ypp, yppp

        integer(I4B) :: klo, khi, n
        integer(I4B) :: k
        real(DP) :: h, p, p1, p2, p3
        real(DP) :: x
        !-----------------------------------------------------------------------
        real(DP) :: delta
        !-----------------------------------------------------------------------

        n = size(a)

        if (.not. (n == size(b) .and. n == size(c) &
                   .and. n == size(d) .and. n == size(xa))) then
            write (*, *) 'splint_horner3: assertion 1 failed'
            print *, size(a), size(b), size(c), size(d), size(xa)
            stop 'program terminated'
        end if

        x = x_in

        klo = 1
        khi = n

        ! Bisection to find k value for which xa(klo) < x < xa(klo+1)
        do while ((khi - klo) .gt. 1)
            k = (khi + klo)/2
            if (xa(k) .gt. x) then
                khi = k
            else
                klo = k
            end if
        end do

        ! Checks to see if bisection was sucessfull.
        if ((klo < 0) .or. (klo > n)) then
            print *, 'splint_horner3: n, klo: ', n, klo
            stop
        end if
        if ((khi < 0) .or. (khi > n)) then
            print *, 'splint_horner3: n, khi: ', n, khi
            stop
        end if

        h = x - xa(klo)

        if (splinecof_compatibility) then
            !  Linear interpolation
            delta = xa(khi) - xa(klo)
            p = a(klo) + h*(b(klo) + delta*(c(klo) + delta*d(klo)))
        else
            ! Nonlinear interpolation
            p = a(klo) + h*(b(klo) + h*(c(klo) + h*d(klo)))
        end if

        y = monom(x, order)*p

        if (swd .ne. 0) then
            p1 = b(klo) + h*(2.0d0*c(klo) + 3.0d0*d(klo)*h)
            p2 = 2.0d0*c(klo) + 6.0d0*d(klo)*h
            p3 = 6.0d0*d(klo)
            yp = monom_deriv(x_in, order)*p + monom(x_in, order)*p1
            ypp = monom_2nd_deriv(x_in, order)*p &
                  + 2.0d0*monom_deriv(x_in, order)*p1 &
                  + monom(x_in, order)*p2
            yppp = monom_3rd_deriv(x_in, order)*p &
                   + 3.0d0*monom_2nd_deriv(x_in, order)*p1 &
                   + 3.0d0*monom_deriv(x_in, order)*p2 &
                   + monom(x_in, order)*p3
        else
            yp = 0.0d0
            ypp = 0.0d0
            yppp = 0.0d0
        end if

    end subroutine spline_1d

    function monom(x, order)
        use nrtype, only: DP

        implicit none

        real(DP), intent(IN) :: x
        real(DP), intent(IN) :: order
        real(DP) :: monom

        if (order .ne. 0.0d0) then
            if (x == 0.0d0) then
                monom = 0.0d0
            else
                monom = x**order
            end if
        else
            monom = 1.0d0
        end if

    end function monom

    function monom_deriv(x, order)
        use nrtype, only: DP

        implicit none

        real(DP), intent(IN) :: x
        real(DP), intent(IN) :: order
        real(DP) :: monom_deriv

        if ((order - 1.0d0) .ne. 0.0d0) then
            if (x == 0.0d0) then
                monom_deriv = 0.0d0
            else
                monom_deriv = order*x**(order - 1.0d0)
            end if
        else
            monom_deriv = 1.0d0
        end if

    end function monom_deriv

    function monom_2nd_deriv(x, order)
        use nrtype, only: DP

        implicit none

        real(DP), intent(IN) :: x
        real(DP), intent(IN) :: order
        real(DP) :: monom_2nd_deriv

        if ((order - 2.0d0) .ne. 0.0d0) then
            if (x == 0.0d0) then
                monom_2nd_deriv = 0.0d0
            else
                monom_2nd_deriv = order*(order - 1.0d0)*x**(order - 2.0d0)
            end if
        else
            monom_2nd_deriv = 1.0d0
        end if

    end function monom_2nd_deriv

    function monom_3rd_deriv(x, order)

        use nrtype, only: DP

        implicit none

        real(DP), intent(IN) :: x
        real(DP), intent(IN) :: order
        real(DP) :: monom_3rd_deriv

        if ((order - 3.0d0) .ne. 0.0d0) then
            if (x == 0.0d0) then
                monom_3rd_deriv = 0.0d0
            else
                monom_3rd_deriv = order*(order - 1.0d0)*(order - 2.0d0)* &
                                  x**(order - 3.0d0)
            end if
        else
            monom_3rd_deriv = 1.0d0
        end if

    end function monom_3rd_deriv

end module neo_spline_over_s
