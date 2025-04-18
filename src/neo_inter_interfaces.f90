
! --------------------------------------------------------------------
!
!  inter_interfaces.f90
!
! --------------------------------------------------------------------

module inter_interfaces

    interface lubksb
        subroutine lubksb_a(a, indx, b)
            use nrtype
            real(DP), dimension(:, :), intent(IN) :: a
            integer(I4B), dimension(:), intent(IN) :: indx
            real(DP), dimension(:), intent(INOUT) :: b
        end subroutine lubksb_a
    end interface

    interface ludcmp
        subroutine ludcmp_a(a, indx, d)
            use nrtype, only: I4B, DP
            real(DP), dimension(:, :), intent(INOUT) :: a
            integer(I4B), dimension(:), intent(OUT) :: indx
            real(DP), intent(OUT) :: d
        end subroutine ludcmp_a
    end interface

    interface splinecof3
        subroutine splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                a, b, c, d, m, f)
            use nrtype, only: I4B, DP
            real(DP), intent(INOUT) :: c1, cn
            real(DP), dimension(:), intent(IN) :: x
            real(DP), dimension(:), intent(IN) :: y
            real(DP), dimension(:), intent(IN) :: lambda1
            integer(I4B), dimension(:), intent(IN) :: indx
            real(DP), dimension(:), intent(OUT) :: a, b, c, d
            integer(I4B), intent(IN) :: sw1, sw2
            real(DP), intent(IN) :: m
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof3_a
    end interface

    interface splinecof1
        subroutine splinecof1_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
            & a, b, c, d, m, f)
            use nrtype, only: I4B, DP

            real(DP), intent(inout) :: c1, cn
            real(DP), dimension(:), intent(in) :: x
            real(DP), dimension(:), intent(in) :: y
            real(DP), dimension(:), intent(in) :: lambda1
            integer(I4B), dimension(:), intent(in) :: indx
            real(DP), dimension(:), intent(out) :: a, b, c, d
            integer(I4B), intent(in) :: sw1, sw2
            real(DP), intent(in) :: m
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof1_a
    end interface

    interface reconstruction3
        subroutine reconstruction3_a(ai, bi, ci, di, h, a, b, c, d)
            use nrtype, only: DP
            real(DP), intent(IN) :: ai, bi, ci, di
            real(DP), intent(IN) :: h
            real(DP), intent(OUT) :: a, b, c, d
        end subroutine reconstruction3_a
    end interface

    interface reconstruction1
        subroutine reconstruction1_a(ai, bi, ci, di, h, a, b, c, d)
            use nrtype, only: DP
            real(DP), intent(in) :: ai, bi, ci, di
            real(DP), intent(in) :: h
            real(DP), intent(out) :: a, b, c, d
        end subroutine reconstruction1_a
    end interface

    interface calc_opt_lambda3
        subroutine calc_opt_lambda3_a(x, y, lambda)
            use nrtype, only: DP
            real(DP), dimension(:), intent(IN) :: x, y
            real(DP), dimension(:), intent(OUT) :: lambda
        end subroutine calc_opt_lambda3_a
    end interface

    interface dist_lin
        subroutine dist_lin_a(x, y, ymax, dist)
            use nrtype, only: DP
            real(DP), dimension(:), intent(IN) :: x, y
            real(DP), intent(IN) :: ymax
            real(DP), intent(OUT) :: dist
        end subroutine dist_lin_a
    end interface

    interface splinecof3_lo_driv
        subroutine splinecof3_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
                                        sw1, sw2, a, b, c, d, m, f)
            use nrtype, only: I4B, DP
            integer(I4B), dimension(:), intent(IN) :: indx
            real(DP), intent(IN) :: m
            real(DP), intent(INOUT) :: c1, cn
            real(DP), dimension(:), intent(IN) :: x
            real(DP), dimension(:), intent(IN) :: y
            real(DP), dimension(:), intent(IN) :: lambda
            integer(I4B), dimension(:), intent(IN) :: w
            real(DP), dimension(:), intent(OUT) :: a, b, c, d
            integer(I4B), intent(IN) :: sw1, sw2
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof3_lo_driv_a
    end interface

    interface splinecof1_lo_driv
        subroutine splinecof1_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
            & sw1, sw2, a, b, c, d, m, f)
            use nrtype, only: I4B, DP
            integer(I4B), dimension(:), intent(in) :: indx
            real(DP), intent(in) :: m
            real(DP), intent(inout) :: c1, cn
            real(DP), dimension(:), intent(in) :: x
            real(DP), dimension(:), intent(in) :: y
            real(DP), dimension(:), intent(in) :: lambda
            integer(I4B), dimension(:), intent(in) :: w
            real(DP), dimension(:), intent(out) :: a, b, c, d
            integer(I4B), intent(in) :: sw1, sw2
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof1_lo_driv_a
    end interface

    interface splinecof3_hi_driv
        subroutine splinecof3_hi_driv_a(x, y, m, a, b, c, d, indx, f)
            use nrtype, only: I4B, DP
            integer(I4B), dimension(:), intent(IN) :: indx
            real(DP), dimension(:), intent(IN) :: m
            real(DP), dimension(:), intent(IN) :: x
            real(DP), dimension(:, :), intent(IN) :: y
            real(DP), dimension(:, :), intent(OUT) :: a, b, c, d
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof3_hi_driv_a
    end interface

    interface splinecof1_hi_driv
        subroutine splinecof1_hi_driv_a(x, y, m, a, b, c, d, indx, f)
            use nrtype, only: I4B, DP
            integer(I4B), dimension(:), intent(in) :: indx
            real(DP), dimension(:), intent(in) :: m
            real(DP), dimension(:), intent(in) :: x
            real(DP), dimension(:, :), intent(in) :: y
            real(DP), dimension(:, :), intent(out) :: a, b, c, d
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: f
                end function f
            end interface
        end subroutine splinecof1_hi_driv_a
    end interface

    interface splint_horner3
        subroutine splint_horner3_a(xa, a, b, c, d, swd, m, x_in, f, fp, fpp, fppp, &
                                    y, yp, ypp, yppp)
            use nrtype, only: I4B, DP
            integer(I4B), intent(IN) :: swd
            real(DP), intent(IN) :: m
            real(DP), dimension(:), intent(IN) :: xa, a, b, c, d
            real(DP), intent(IN) :: x_in
            real(DP), intent(OUT) :: y, yp, ypp, yppp
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: f
                end function f
            end interface
            interface
                function fp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fp
                end function fp
            end interface
            interface
                function fpp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fpp
                end function fpp
            end interface
            interface
                function fppp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fppp
                end function fppp
            end interface
        end subroutine splint_horner3_a
    end interface

    interface splint_horner1
        subroutine splint_horner1_a(xa, a, b, c, d, swd, m, x_in, &
            & f, fp, fpp, fppp, y, yp, ypp, yppp)
            use nrtype, only: I4B, DP

            integer(I4B), intent(in) :: swd
            real(DP), intent(in) :: m
            real(DP), dimension(:), intent(in) :: xa, a, b, c, d
            real(DP), intent(in) :: x_in
            real(DP), intent(out) :: y, yp, ypp, yppp
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: f
                end function f
            end interface
            interface
                function fp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: fp
                end function fp
            end interface
            interface
                function fpp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: fpp
                end function fpp
            end interface
            interface
                function fppp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(in) :: x
                    real(DP), intent(in) :: m
                    real(DP) :: fppp
                end function fppp
            end interface
        end subroutine splint_horner1_a
    end interface

    interface splint_horner3_driv_s
    subroutine splint_horner3_driv_s_a(svec, a, b, c, d, swd, ixm, ixn, s, theta, phi, &
                                           f, fp, fpp, fppp, y, ys, yt, yp)
            use nrtype, only: I4B, DP
            real(DP), dimension(:), intent(IN) :: svec
            real(DP), dimension(:, :), intent(IN) :: a, b, c, d
            integer(I4B), intent(IN) :: swd
            real(DP), dimension(:), intent(IN) :: ixm
            integer(I4B), dimension(:), intent(IN) :: ixn
            real(DP), intent(IN) :: s, theta, phi
            real(DP), intent(OUT) :: y, ys, yt, yp
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: f
                end function f
            end interface
            interface
                function fp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fp
                end function fp
            end interface
            interface
                function fpp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fpp
                end function fpp
            end interface
            interface
                function fppp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fppp
                end function fppp
            end interface
        end subroutine splint_horner3_driv_s_a
    end interface

    interface splint_horner3_driv_c
    subroutine splint_horner3_driv_c_a(svec, a, b, c, d, swd, ixm, ixn, s, theta, phi, &
                                           f, fp, fpp, fppp, y, ys, yt, yp)
            use nrtype, only: I4B, DP
            real(DP), dimension(:), intent(IN) :: svec
            real(DP), dimension(:, :), intent(IN) :: a, b, c, d
            integer(I4B), intent(IN) :: swd
            real(DP), dimension(:), intent(IN) :: ixm
            integer(I4B), dimension(:), intent(IN) :: ixn
            real(DP), intent(IN) :: s, theta, phi
            real(DP), intent(OUT) :: y, ys, yt, yp
            interface
                function f(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: f
                end function f
            end interface
            interface
                function fp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fp
                end function fp
            end interface
            interface
                function fpp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fpp
                end function fpp
            end interface
            interface
                function fppp(x, m)
                    use nrtype, only: DP
                    implicit none
                    real(DP), intent(IN) :: x
                    real(DP), intent(IN) :: m
                    real(DP) :: fppp
                end function fppp
            end interface
        end subroutine splint_horner3_driv_c_a
    end interface

    interface
        function tf(x, m)
            use nrtype, only: DP
            implicit none
            real(DP), intent(IN) :: x
            real(DP), intent(IN) :: m
            real(DP) :: tf
        end function tf
    end interface
    interface
        function tfp(x, m)
            use nrtype, only: DP
            implicit none
            real(DP), intent(IN) :: x
            real(DP), intent(IN) :: m
            real(DP) :: tfp
        end function tfp
    end interface
    interface
        function tfpp(x, m)
            use nrtype, only: DP
            implicit none
            real(DP), intent(IN) :: x
            real(DP), intent(IN) :: m
            real(DP) :: tfpp
        end function tfpp
    end interface
    interface
        function tfppp(x, m)
            use nrtype, only: DP
            implicit none
            real(DP), intent(IN) :: x
            real(DP), intent(IN) :: m
            real(DP) :: tfppp
        end function tfppp
    end interface

    interface
        function tfone(x, m)
            use nrtype, only: DP
            implicit none
            real(DP), intent(IN) :: x
            real(DP), intent(IN) :: m
            real(DP) :: tfone
        end function tfone
    end interface
    interface
        function tfzero(x, m)
            use nrtype, only: DP
            implicit none
            real(DP), intent(IN) :: x
            real(DP), intent(IN) :: m
            real(DP) :: tfzero
        end function tfzero
    end interface

end module inter_interfaces
