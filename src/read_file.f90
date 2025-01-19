module read_file
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    abstract interface
        subroutine read_field_file(file_name, B)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            character(len=256), intent(in) :: file_name
            real(dp), intent(out) :: B
        end subroutine read_field_file
    end interface

contains

    subroutine read_boozer_file(file_name, B)
        use constants, only: pi

        character(len=*), intent(in) :: file_name
        real(dp), intent(out) :: B

        integer :: m_max, n_max, n_s, nfp
        integer :: num_pol_modes, num_tor_modes, num_modes
        real(dp) :: tor_flux_separatrix, psi_tor_separatrix
        integer :: file_unit

        integer, dimension(:), allocatable :: m, n
        real(dp), dimension(:), allocatable :: s_tor, iota, curr_pol, curr_tor
        real(dp), dimension(:), allocatable :: dp_ds, sqrtg_00
        real(dp), dimension(:, :), allocatable :: rmnc, zmns, vmns, bmnc
        character(len=256) :: dummy
        integer :: extra_count, i, j
        logical :: extra_zero
        integer :: ios, i_alloc

        open (newunit=file_unit, file=file_name)
        call skip_lines(file_unit, 5)
    read (file_unit,"(I6, I6, I6, I6, F8.8)", iostat=ios) m_max, n_max, n_s, nfp, tor_flux_separatrix
        num_pol_modes = m_max + 1
        num_tor_modes = 2*n_max + 1
        num_modes = num_pol_modes*num_tor_modes

        allocate (m(num_modes), n(num_modes), stat=i_alloc)
        if (i_alloc /= 0) stop 'Allocation for integer arrays failed!'

    allocate(s_tor(n_s), iota(n_s), curr_pol(n_s), curr_tor(n_s), dp_ds(n_s), sqrtg_00(n_s), stat = i_alloc)
        if (i_alloc /= 0) stop 'Allocation for real arrays failed!'

    allocate(rmnc(n_s,num_modes), zmns(n_s,num_modes), vmns(n_s,num_modes), bmnc(n_s,num_modes), stat = i_alloc)
        if (i_alloc /= 0) stop 'Allocation for fourier arrays (1) failed!'
        do i = 1, n_s
            read (file_unit, *) dummy
            read (file_unit, *) dummy
            read (file_unit, *) s_tor(i), iota(i), curr_pol(i), curr_tor(i), &
                dp_ds(i), sqrtg_00(i)
            read (file_unit, *) dummy

            extra_zero = .false.
            extra_count = 0
            do j = 1, num_modes
                if (j .gt. 1) then
                    if (m(j - 1) .eq. 0 .and. n(j - 1) .eq. 0) then
                        extra_zero = .true.
                    end if
                end if
                if (extra_zero) then
                    extra_count = extra_count + 1
                    if (extra_count .eq. n_max) extra_zero = .false.
                    m(j) = 0
                    n(j) = -extra_count
                    rmnc(i, j) = 0.0d0
                    zmns(i, j) = 0.0d0
                    vmns(i, j) = 0.0d0
                    bmnc(i, j) = 0.0d0
                else
                    read (file_unit, *) m(j), n(j), &
                        rmnc(i, j), zmns(i, j), vmns(i, j), &
                        bmnc(i, j)
                end if
            end do
        end do

        !**********************************************************
        ! Change from Gernot Kapper - 02.12.2015
        ! This corrects the direction of the poloidal current to
        ! match the Boozer file w7x-m24li.bc
        !**********************************************************
        ! curr_pol = - curr_pol * 2.d-7 * nfp   ! ? -   ! Before patch
        curr_pol = curr_pol*2.d-7*nfp               ! After patch

        !**********************************************************
        ! Patch from Gernot Kapper - 20.11.2014
        ! See mail from Winfried Kernbichler archived at
        ! /proj/plasma/doCUMENTS/Neo2/Archive/
        !**********************************************************
        ! curr_tor = curr_tor * 2.d-7 * nfp   ! Henning   ! Before patch
        curr_tor = -curr_tor*2.d-7         ! Henning   ! After patch

        n_max = n_max*nfp
        n = n*nfp
        m = m
        psi_tor_separatrix = abs(tor_flux_separatrix)/2.0_dp*pi
    end subroutine read_boozer_file

    subroutine skip_lines(file_unit, n_lines)
        integer, intent(in) :: file_unit, n_lines

        integer :: ios, line

        do line = 1, n_lines
            read (file_unit, '(A)', iostat=ios)
            if (ios /= 0) then
                print *, "Error when skipping lines."
                error stop
            end if
        end do
    end subroutine

end module read_file
