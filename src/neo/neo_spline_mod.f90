module neo_spline_mod

contains

    subroutine spl2d(nx, ny, hx, hy, mx, my, f, spl)
        ! Makes a 2-dimensional cubic spline of function f(x,y)
        !
        ! Input:  nx, ny              number of values in x and y
        !         hx, hy              step size in x and y (aequidistant)
        !         mx, my              spline mode (0: standard, 1: periodic)
        !         f(nx,ny)            f(x,y)-values
        ! Output: spl                 Array with spline parameters

        use nrtype

        implicit none

        integer, intent(in) :: nx, ny, mx, my
        real(kind=dp), intent(in) :: hx, hy
        real(kind=dp), dimension(nx, ny), intent(in) :: f
        real(kind=dp), dimension(4, 4, nx, ny), intent(out) :: spl

        real(kind=dp), dimension(:), allocatable :: bi, ci, di, s
        integer :: i, j, k

        allocate (bi(nx), ci(nx), di(nx), s(nx))
        do j = 1, ny
            do i = 1, nx
                s(i) = f(i, j)
            end do
            if (mx .eq. 0) then
                call splreg(nx, hx, s, bi, ci, di)
            else
                call splper(nx, hx, s, bi, ci, di)
            end if
            do i = 1, nx
                spl(1, 1, i, j) = s(i)
                spl(2, 1, i, j) = bi(i)
                spl(3, 1, i, j) = ci(i)
                spl(4, 1, i, j) = di(i)
            end do
        end do
        deallocate (bi, ci, di, s)

        allocate (bi(ny), ci(ny), di(ny), s(ny))
        do k = 1, 4
            do i = 1, nx
                do j = 1, ny
                    s(j) = spl(k, 1, i, j)
                end do
                if (my .eq. 0) then
                    call splreg(ny, hy, s, bi, ci, di)
                else
                    call splper(ny, hy, s, bi, ci, di)
                end if
                do j = 1, ny
                    spl(k, 2, i, j) = bi(j)
                    spl(k, 3, i, j) = ci(j)
                    spl(k, 4, i, j) = di(j)
                end do
            end do
        end do
        deallocate (bi, ci, di, s)

    end subroutine spl2d

    !=====================================================
    subroutine eva2d(nx, ny, ix, iy, dx, dy, spl, spval)
        ! Evaluates a 2-dimensional cubic spline of function f(x,y)
        !
        ! Input:  nx, ny              number of values in x and y
        !         ix, iy              pointer into the spline array spl
        !         dx, dy              distance from x(ix) and y(iy)
        !         spl                 array with spline data
        ! Output: spval               evaluated function value

        use nrtype

        implicit none

        integer, intent(in) :: nx, ny, ix, iy
        real(kind=dp), intent(in) :: dx, dy
        real(kind=dp), dimension(4, 4, nx, ny), intent(in) :: spl
        real(kind=dp), intent(out) :: spval

        real(kind=dp), dimension(4) :: a
        integer :: l

        do l = 1, 4
            a(l) = spl(1, l, ix, iy) + dx*(spl(2, l, ix, iy) + &
                                          dx*(spl(3, l, ix, iy) + dx*spl(4, l, ix, iy)))
        end do
        spval = a(1) + dy*(a(2) + dy*(a(3) + dy*a(4)))

    end subroutine eva2d

    !=====================================================
    subroutine splreg(n, h, y, bi, ci, di)
        ! Makes a cubic spline of function y(x)
        !
        ! Input:  n                   number of values in y
        !         h                   step size in x (aequidistant)
        !         y(n)                y-values
        ! Output: bi(n),ci(n),di(n)   Spline parameters

        use nrtype

        implicit none

        integer, intent(in) :: n
        real(kind=dp), intent(in) :: h
        real(kind=dp), dimension(n), intent(in) :: y
        real(kind=dp), dimension(n), intent(out) :: bi, ci, di

        real(kind=dp) :: ak1, ak2, am1, am2, c, e, c1
        real(kind=dp), dimension(:), allocatable :: al, bt
        integer :: k, n2, i, i5

        allocate (al(n), bt(n))

        ak1 = 0.d0
        ak2 = 0.d0
        am1 = 0.d0
        am2 = 0.d0
        k = n - 1
        al(1) = ak1
        bt(1) = am1
        n2 = n - 2
        c = -4.d0*h
        do i = 1, n2
            e = -3.d0*((y(i + 2) - y(i + 1)) - (y(i + 1) - y(i)))/h
            c1 = c - al(i)*h
            al(i + 1) = h/c1
            bt(i + 1) = (h*bt(i) + e)/c1
        end do
        ci(n) = (am2 + ak2*bt(k))/(1.d0 - al(k)*ak2)
        do i = 1, k
            i5 = n - i
            ci(i5) = al(i5)*ci(i5 + 1) + bt(i5)
        end do
        n2 = n - 1
        do i = 1, n2
            bi(i) = (y(i + 1) - y(i))/h - h*(ci(i + 1) + 2.d0*ci(i))/3.d0
            di(i) = (ci(i + 1) - ci(i))/h/3.d0
        end do
        deallocate (al, bt)

    end subroutine splreg

    !=====================================================
    subroutine splper(n, h, y, bi, ci, di)
        ! Makes a cubic spline of periodic function y(x)
        !
        ! Input:  n                   number of values in y
        !         h                   step size in x (aequidistant)
        !         y(n)                y-values
        ! Output: bi(n),ci(n),di(n)   Spline parameters

        use nrtype

        implicit none

        integer, intent(in) :: n
        real(kind=dp), intent(in) :: h
        real(kind=dp), dimension(n), intent(in) :: y
        real(kind=dp), dimension(n), intent(out) :: bi, ci, di

        real(kind=dp) :: psi, ss
        real(kind=dp), dimension(:), allocatable :: bmx, yl
        real(kind=dp), dimension(:), allocatable :: amx1, amx2, amx3
        integer :: nmx, n1, n2, i, i1

        allocate (bmx(n), yl(n), amx1(n), amx2(n), amx3(n))

        bmx(1) = 1.d30

        nmx = n - 1
        n1 = nmx - 1
        n2 = nmx - 2
        psi = 3.d0/h/h

        call spfper(n, amx1, amx2, amx3)

        bmx(nmx) = (y(nmx + 1) - 2.d0*y(nmx) + y(nmx - 1))*psi
        bmx(1) = (y(2) - y(1) - y(nmx + 1) + y(nmx))*psi
        do i = 3, nmx
            bmx(i - 1) = (y(i) - 2.d0*y(i - 1) + y(i - 2))*psi
        end do
        yl(1) = bmx(1)/amx1(1)
        do i = 2, n1
            i1 = i - 1
            yl(i) = (bmx(i) - yl(i1)*amx2(i1))/amx1(i)
        end do
        ss = 0.d0
        do i = 1, n1
            ss = ss + yl(i)*amx3(i)
        end do
        yl(nmx) = (bmx(nmx) - ss)/amx1(nmx)
        bmx(nmx) = yl(nmx)/amx1(nmx)
        bmx(n1) = (yl(n1) - amx2(n1)*bmx(nmx))/amx1(n1)
        do i = n2, 1, -1
            bmx(i) = (yl(i) - amx3(i)*bmx(nmx) - amx2(i)*bmx(i + 1))/amx1(i)
        end do
        do i = 1, nmx
            ci(i) = bmx(i)
        end do

        do i = 1, n1
            bi(i) = (y(i + 1) - y(i))/h - h*(ci(i + 1) + 2.d0*ci(i))/3.d0
            di(i) = (ci(i + 1) - ci(i))/h/3.d0
        end do
        bi(nmx) = (y(n) - y(n - 1))/h - h*(ci(1) + 2.d0*ci(nmx))/3.d0
        di(nmx) = (ci(1) - ci(nmx))/h/3.d0

        ! Fix of problems at upper periodicity boundary
        bi(n) = bi(1)
        ci(n) = ci(1)
        di(n) = di(1)

        deallocate (bmx, yl, amx1, amx2, amx3)

    end subroutine splper

    !=====================================================
    subroutine spfper(np1, amx1, amx2, amx3)
        ! Helper routine for splfi

        use nrtype

        implicit none

        integer, intent(in) :: np1
        real(kind=dp), dimension(np1), intent(out) :: amx1, amx2, amx3
        real(kind=dp) :: beta, ss
        integer :: n, n1, i, i1

        n = np1 - 1

        n1 = n - 1
        amx1(1) = 2.d0
        amx2(1) = 0.5d0
        amx3(1) = 0.5d0
        amx1(2) = sqrt(15.d0)/2.d0
        amx2(2) = 1.d0/amx1(2)
        amx3(2) = -.25d0/amx1(2)
        beta = 3.75d0
        do i = 3, n1
            i1 = i - 1
            beta = 4.d0 - 1.d0/beta
            amx1(i) = sqrt(beta)
            amx2(i) = 1.d0/amx1(i)
            amx3(i) = -amx3(i1)/amx1(i)/amx1(i1)
        end do
        amx3(n1) = amx3(n1) + 1.d0/amx1(n1)
        amx2(n1) = amx3(n1)
        ss = 0.0d0
        do i = 1, n1
            ss = ss + amx3(i)*amx3(i)
        end do
        amx1(n) = sqrt(4.d0 - ss)

    end subroutine spfper

    !=====================================================
    subroutine poi2d(hx, hy, mx, my, &
                     xmin, xmax, ymin, ymax, &
                     x, y, ix, iy, dx, dy, ierr)
        ! Creates Pointers for eva2d
        !
        ! Input:  hx, hy              increment in x and y
        !         mx, my              standard (0) or periodic (1) spline
        !         xmin, xmax          Minimum and maximum x
        !         ymin, ymax          Minimum and maximum y
        !         x, y                x and y values for spline avaluation
        ! Output: spval               evaluated function value
        !         ix, iy              pointer into the spline array spl
        !         dx, dy              distance from x(ix) and y(iy)
        !         ierr                error (> 0)

        use nrtype

        implicit none

        real(kind=dp), intent(in) :: hx, hy
        integer, intent(in) :: mx, my
        real(kind=dp), intent(in) :: xmin, xmax, ymin, ymax
        real(kind=dp), intent(in) :: x, y

        integer, intent(out) :: ix, iy
        real(kind=dp), intent(out) :: dx, dy
        integer, intent(out) :: ierr

        real(kind=dp) :: dxx, x1, dyy, y1
        real(kind=dp) :: dxmax, dymax

        ierr = 0

        dxx = x - xmin
        if (mx .eq. 0) then
            if (dxx .lt. 0.d0) then
                ierr = 1
                return
            end if
            if (x .gt. xmax) then
                ierr = 2
                return
            end if
        else
            dxmax = xmax - xmin
            if (dxx .lt. 0.d0) then
                dxx = dxx + dble(float(1 + int(abs(dxx/dxmax))))*dxmax
            else if (dxx .gt. dxmax) then
                dxx = dxx - dble(float(int(abs(dxx/dxmax))))*dxmax
            end if
        end if
        x1 = dxx/hx
        ix = int(x1)
        dx = hx*(x1 - dble(float(ix)))
        ix = ix + 1

        dyy = y - ymin
        if (my .eq. 0) then
            if (dyy .lt. 0.d0) then
                ierr = 3
                return
            end if
            if (y .gt. ymax) then
                ierr = 4
                return
            end if
        else
            dymax = ymax - ymin
            if (dyy .lt. 0.d0) then
                dyy = dyy + dble(float(1 + int(abs(dyy/dymax))))*dymax
            else if (dyy .gt. dymax) then
                dyy = dyy - dble(float(int(abs(dyy/dymax))))*dymax
            end if
        end if
        y1 = dyy/hy
        iy = int(y1)
        dy = hy*(y1 - dble(float(iy)))
        iy = iy + 1

    end subroutine poi2d

end module neo_spline_mod
