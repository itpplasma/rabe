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
    character(len=256), intent(in) :: file_name
    real(dp), intent(out) :: B 
    
    integer :: m_max, n_max, n_s, nfp
    integer :: num_pol_modes, num_tor_modes, num_modes
    real(dp) :: tor_flux_separatrix
    integer :: file_unit

    integer :: ios

    open(newunit=file_unit, file=file_name)
    call skip_lines(file_unit, 5)
    read (file_unit,"(I6, I6, I6, I6, F8.8)", iostat=ios) m_max, n_max, n_s, nfp, tor_flux_separatrix
    num_pol_modes = m_max+1
    num_tor_modes = 2*n_max+1
    num_modes = num_pol_modes*num_tor_modes

      ! **********************************************************************
      ! Allocate storage arrays
      ! **********************************************************************
      ALLOCATE(ixm(num_modes), ixn(num_modes), stat = i_alloc)
      if(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(pixm(num_modes), pixn(num_modes), stat = i_alloc)
      if(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

      ALLOCATE(i_m(num_pol_modes), i_n(num_tor_modes), stat = i_alloc)
      if(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

      ALLOCATE(es(n_s), iota(n_s), curr_pol(n_s), curr_tor(n_s),               &
          pprime(n_s), sqrtg00(n_s), b00(n_s), stat = i_alloc)
      if(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

      ALLOCATE(rmnc(n_s,num_modes), zmnc(n_s,num_modes), lmnc(n_s,num_modes),             &
          bmnc(n_s,num_modes),                                                 &
          stat = i_alloc)
      if(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
      !***********************************************************************
      ! Read input arrays
      !***********************************************************************
      do i =1, n_s
        read(file_unit,*) dummy
        read(file_unit,*) dummy
        read(file_unit,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
             pprime(i),sqrtg00(i)
        read(file_unit,*) dummy

        extra_zero = .FALSE.
        extra_count = 0
        do j=1,num_modes
          if (j .GT. 1) THEN
            if (ixm(j-1) .EQ. 0 .AND. ixn(j-1) .EQ. 0) THEN
              extra_zero = .TRUE.
            ENDIF
          END if
          if (extra_zero) THEN
            extra_count =  extra_count + 1
            if (extra_count .EQ. n_max) extra_zero = .FALSE.
            ixm(j) = 0
            ixn(j) = -extra_count
            rmnc(i,j) = 0.0d0
            zmnc(i,j) = 0.0d0
            lmnc(i,j) = 0.0d0
            bmnc(i,j) = 0.0d0
          ELSE
            read(file_unit,*) ixm(j),ixn(j),                                    &
                 rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                 bmnc(i,j)
          ENDIF
        END do
      END do

        !**********************************************************
      ! Change from Gernot Kapper - 02.12.2015
      ! This corrects the direction of the poloidal current to
      ! match the Boozer file w7x-m24li.bc
      !**********************************************************
      ! curr_pol = - curr_pol * 2.d-7 * nfp   ! ? -   ! Before patch
      curr_pol = curr_pol * 2.d-7 * nfp               ! After patch

      !**********************************************************
      ! Patch from Gernot Kapper - 20.11.2014
      ! See mail from Winfried Kernbichler archived at
      ! /proj/plasma/doCUMENTS/Neo2/Archive/
      !**********************************************************
      ! curr_tor = curr_tor * 2.d-7 * nfp   ! Henning   ! Before patch
      curr_tor = - curr_tor * 2.d-7         ! Henning   ! After patch

      max_n_mode = max_n_mode * nfp
      ixn =  ixn * nfp
      i_n =  i_n * nfp
      ixm =  ixm
      i_m =  i_m
      psi_pr = ABS(tor_flux_separatrix) / twopi
end subroutine read_boozer_file

subroutine skip_lines(file_unit, n_lines)
    integer, intent(in) :: file_unit, n_lines

    integer :: ios, line

    do line = 1, n_lines
        read(file_unit, '(A)', iostat=ios)
        if (ios /= 0) THEN
            print *, "Error when skipping lines."
            error stop
        END if
    END do
end subroutine

end module read_file
