! AUTHOR: Bernhard Seiwald
!
! DATE:   18.07.2001

! simple wrapper for solvers for real system of linear
! equations  A * X = B

module solve_systems

    use nrtype, only: I4B, DP

    implicit none

    public :: solve_eqsys

! --------------------------------------------------------------------

contains

    subroutine solve_eqsys(a, b, info)

        implicit none

        real(DP), dimension(:, :), intent(INOUT) :: a
        real(DP), dimension(:), intent(INOUT) :: b
        integer(I4B), intent(OUT) :: info
        integer(I4B) :: i_alloc
        integer(I4B) :: n, nrhs, lda, ldb
        integer(I4B), dimension(:), allocatable :: ipiv
! --------------------------------------------------------------------

        lda = size(a, 1)
        n = size(a, 2)
        ldb = size(b, 1)
        nrhs = 1

        allocate (ipiv(n), stat=i_alloc)
        if (i_alloc /= 0) stop 'solve_eqsys: Allocation for array failed!'

        call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

        info = 0

        deallocate (ipiv, stat=i_alloc)
        if (i_alloc /= 0) stop 'solve_eqsys: Deallocation for array failed!'

    end subroutine solve_eqsys

end module solve_systems
