program compare_rabe_deriv
    !! rabe_deriv.x run on a field and rerun on its own output as spectrum file
    !! must give the same values, Jacobian and lambda_off (NetCDF round trip of
    !! the documented layout).
    use constants, only: dp
    use netcdf_mod, only: netcdf_t

    implicit none

    real(dp), parameter :: tol = 1e-10_dp
    type(netcdf_t) :: a, b
    integer :: n_surf, n_param, n_out, n_nu
    real(dp), allocatable :: va(:, :), vb(:, :), ja(:, :), jb(:, :)
    real(dp), allocatable :: la(:, :), lb(:, :)
    real(dp) :: err_v, err_j, err_l

    call a%open("from_field/rabe_deriv.nc")
    call b%open("from_spectrum/rabe_deriv.nc")
    call a%read_dim("surface", n_surf)
    call a%read_dim("param", n_param)
    call a%read_dim("output", n_out)
    call a%read_dim("nu", n_nu)
    allocate (va(n_out, n_surf), vb(n_out, n_surf))
    allocate (ja(n_param, n_out*n_surf), jb(n_param, n_out*n_surf))
    allocate (la(n_nu, n_surf), lb(n_nu, n_surf))
    call a%read_real_2d("values", va)
    call b%read_real_2d("values", vb)
    call read_jacobian(a, ja)
    call read_jacobian(b, jb)
    call a%read_real_2d("lambda_off", la)
    call b%read_real_2d("lambda_off", lb)
    call a%close()
    call b%close()

    err_v = maxval(abs(va - vb)/abs(va))
    err_j = maxval(abs(ja - jb))/maxval(abs(ja))
    err_l = maxval(abs(la - lb)/abs(la))
    print *, "round trip rel errors: values ", err_v, " jacobian ", err_j, &
        " lambda_off ", err_l
    if (err_v > tol .or. err_j > tol .or. err_l > tol .or. n_param < 10) then
        print *, "FAILED: rabe_deriv.nc round trip"
        error stop 1
    end if

contains

    subroutine read_jacobian(nc, jac)
        type(netcdf_t), intent(inout) :: nc
        real(dp), intent(out) :: jac(:, :)
        real(dp), allocatable :: jac3(:, :, :)

        allocate (jac3(size(jac, 1), n_out, n_surf))
        call nc%read_real_3d("jacobian", jac3)
        jac = reshape(jac3, shape(jac))
    end subroutine read_jacobian

end program compare_rabe_deriv
