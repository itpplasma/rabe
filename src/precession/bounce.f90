module bounce
    use omp_lib
    use constants, only: dp, pi
    use params, only: ntimstep, nplagr, nder, npl_half, nfp, dtaumin, relerr, ntau

    implicit none

contains

    subroutine trace_orbit_till_bounce(z_start, z_end)
        use plag_coeff_sub, only: plag_coeff
        use orbit, only: orbit_timestep_axis
        use field_instance, only: magfie

        real(dp), dimension(5), intent(in) :: z_start
        real(dp), dimension(5), intent(out) :: z_end

        integer :: ierr, ierr_coll
        real(dp), dimension(5) :: z
        real(dp) :: bmod, sqrtg
        real(dp), dimension(3) :: bder, hcovar, hctrvr, hcurl
        integer :: it, ktau
        integer(8) :: kt
        logical :: passing

        integer :: ifp_tip, ifp_per
        integer, dimension(:), allocatable :: ipoi
        real(dp), dimension(:), allocatable :: xp
        real(dp), dimension(:, :), allocatable :: coef, orb_sten
        integer, parameter :: n_tip_vars = 6
        real(dp), parameter :: zerolam = 0.0_dp
        real(dp), dimension(n_tip_vars) :: var_tip
        real(dp) :: phiper, alam_prev, par_inv
        integer :: iper, itip, kper, nfp_tip, nfp_per

        integer, parameter :: nfp_dim = 3
        integer :: nfp_cot, ideal, ijpar, ierr_cot, iangvar
        real(dp), dimension(nfp_dim) :: fpr_in
        logical :: did_bounce

        did_bounce = .false.
        z_end = 0d0
        !
        iangvar = 2

        z = z_start

        call magfie(z(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        !$omp critical
        if (.not. allocated(ipoi)) &
          allocate (ipoi(nplagr), coef(0:nder, nplagr), orb_sten(6, nplagr), xp(nplagr))
        !$omp end critical
        do it = 1, nplagr
            ipoi(it) = it
        end do

        nfp_tip = nfp             !<= initial array dimension for tips
        nfp_per = nfp             !<= initial array dimension for periods

        ifp_tip = 0               !<= initialize footprint counter on tips
        ifp_per = 0               !<= initialize footprint counter on periods

        phiper = 0.0d0

        kt = 0

        itip = npl_half + 1
        alam_prev = z(5)

        ! End initialize tip detector
        !--------------------------------
        !
        par_inv = 0d0
        TIMELOOP: do ktau = 1, ntau*(ntimstep - 1)
            call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
            if (ierr .ne. 0) then
                error stop "Orbit lost during integration!"
            end if
            kt = kt + 1
            par_inv = par_inv + z(5)**2*dtaumin ! parallel adiabatic invariant
            if (kt .le. nplagr) then          !<=first nplagr points to initialize stencil
                orb_sten(1:5, kt) = z
                orb_sten(6, kt) = par_inv
            else                          !<=normal case, shift stencil
                orb_sten(1:5, ipoi(1)) = z
                orb_sten(6, ipoi(1)) = par_inv
                ipoi = cshift(ipoi, 1)
            end if

            ! Tip detection and interpolation
            if (alam_prev .lt. 0.d0 .and. z(5) .gt. 0.d0) itip = 0   !<=tip has been passed
            itip = itip + 1
            alam_prev = z(5)
            if (kt .gt. nplagr) then          !<=use only initialized stencil
                if (itip .eq. npl_half) then   !<=stencil around tip is complete, interpolate
                    xp = orb_sten(5, ipoi)

                    call plag_coeff(nplagr, nder, zerolam, xp, coef)

                    var_tip = matmul(orb_sten(:, ipoi), coef(0, :))
                    var_tip(2) = modulo(var_tip(2), 2.0_dp*pi)
                    var_tip(3) = modulo(var_tip(3), 2.0_dp*pi)

                    z = var_tip(1:5) !<= update z to interpolated tip value
                    did_bounce = .true.
                    exit TIMELOOP
                end if
            end if
            ! End tip detection and interpolation
        end do TIMELOOP

        !$omp critical
        z_end = z
        !$omp end critical
        if (.not. did_bounce) then
            print *, "Warning: trace_orbit_till_bounce reached end of time loop without detecting bounce."
        end if
    end subroutine trace_orbit_till_bounce

end module bounce
