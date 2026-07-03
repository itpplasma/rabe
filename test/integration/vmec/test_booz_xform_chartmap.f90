!> Chartmap loader test for boozer_field_t.
!>
!> Loads input/circ_chartmap.nc via init_from_chartmap and checks that Bmod
!> evaluated at grid nodes reproduces the Bmod values stored in the same file.
!> The reference values are read from the file itself, so the test pins the
!> field interface, not the provenance of the fixture.
program test_booz_xform_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: chartmap_file = 'input/circ_chartmap.nc'

    integer, parameter :: n_points = 3
    ! Grid indices (rho, theta, zeta), 1-based, spread over the volume
    integer, parameter :: ir(n_points) = [10, 5, 15]
    integer, parameter :: it(n_points) = [13, 7, 19]
    integer, parameter :: iz(n_points) = [13, 7, 7]

    ! Interpolating splines reproduce grid-node values exactly (measured 0.0
    ! reldiff locally); 1e-9 leaves headroom for compiler/FMA variation.
    real(dp), parameter :: reltol = 1.0e-9_dp
    real(dp), parameter :: abstol = 1.0e-12_dp
    ! Chartmap files are CGS (Bmod in Gauss); boozer_field_t evaluates in SI.
    real(dp), parameter :: gauss2tesla = 1.0e-4_dp

    type(boozer_field_t) :: field
    real(dp), allocatable :: rho_grid(:), theta_grid(:), zeta_grid(:)
    real(dp), allocatable :: Bmod_grid(:, :, :)
    real(dp) :: Bmod
    integer :: point
    logical :: test_failed

    call read_chartmap_reference(chartmap_file, rho_grid, theta_grid, &
                                 zeta_grid, Bmod_grid)
    Bmod_grid = Bmod_grid*gauss2tesla

    call field%init_from_chartmap(chartmap_file)
    if (.not. field%initialized) error stop 'chartmap field not initialized'

    test_failed = .false.
    do point = 1, n_points
        call field%fix_to_surface(rho_grid(ir(point))**2)
        call field%compute_B_mod(theta_grid(it(point)), &
                                 zeta_grid(iz(point)), Bmod)
        if (not_same(Bmod, Bmod_grid(ir(point), it(point), iz(point)), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, 'FAIL: chartmap Bmod mismatch at grid node', &
                ir(point), it(point), iz(point)
            print *, '  computed: ', Bmod
            print *, '  file:     ', Bmod_grid(ir(point), it(point), iz(point))
            print *, '  reldiff:  ', abs(Bmod - Bmod_grid(ir(point), it(point), &
                                                          iz(point))) &
                /Bmod_grid(ir(point), it(point), iz(point))
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

contains

    subroutine read_chartmap_reference(path, rho, theta, zeta, Bmod_ref)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: rho(:), theta(:), zeta(:)
        real(dp), allocatable, intent(out) :: Bmod_ref(:, :, :)

        integer :: ncid, dimid, varid
        integer :: nrho, ntheta, nzeta

        call check_nc(nf90_open(path, nf90_nowrite, ncid), 'open chartmap')
        call check_nc(nf90_inq_dimid(ncid, 'rho', dimid), 'inq rho dim')
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=nrho), 'rho len')
        call check_nc(nf90_inq_dimid(ncid, 'theta', dimid), 'inq theta dim')
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=ntheta), &
                      'theta len')
        call check_nc(nf90_inq_dimid(ncid, 'zeta', dimid), 'inq zeta dim')
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=nzeta), 'zeta len')

        allocate (rho(nrho), theta(ntheta), zeta(nzeta))
        allocate (Bmod_ref(nrho, ntheta, nzeta))

        call check_nc(nf90_inq_varid(ncid, 'rho', varid), 'inq rho')
        call check_nc(nf90_get_var(ncid, varid, rho), 'read rho')
        call check_nc(nf90_inq_varid(ncid, 'theta', varid), 'inq theta')
        call check_nc(nf90_get_var(ncid, varid, theta), 'read theta')
        call check_nc(nf90_inq_varid(ncid, 'zeta', varid), 'inq zeta')
        call check_nc(nf90_get_var(ncid, varid, zeta), 'read zeta')
        call check_nc(nf90_inq_varid(ncid, 'Bmod', varid), 'inq Bmod')
        call check_nc(nf90_get_var(ncid, varid, Bmod_ref), 'read Bmod')
        call check_nc(nf90_close(ncid), 'close chartmap')
    end subroutine read_chartmap_reference

    subroutine check_nc(status, context)
        integer, intent(in) :: status
        character(len=*), intent(in) :: context

        if (status /= nf90_noerr) then
            print *, trim(context), ': ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check_nc

end program test_booz_xform_chartmap
