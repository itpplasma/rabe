program test_against_neo2_anti_sigma
    use neo_field, only: neo_field_t
    use constants, only: dp, pi
    use myplot_module, only: myplot
    use utils, only: linspace

    implicit none

    real(dp) :: retol = 1e-11
    character(len=*), parameter :: bc_filename = "input/helical_anti.bc"

    type(neo_field_t) :: field
    integer, parameter :: n_cases = 4
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: bmod
    real(dp) :: bmod_neo2(n_cases)
    real(dp), parameter :: B_max = 1.05_dp
    integer :: case

    logical :: test_failed

    integer, parameter :: n_stor = 20
    real(dp), dimension(n_stor) :: stor_plot, B_at_chi_pi
    type(myplot) :: plt

    test_failed = .false.

    stor = (/0.50_dp, 0.70_dp, 0.80_dp, 0.9999_dp/)
    theta = (/0.0_dp, 0.25_dp*pi, 2.0_dp*pi, 0.0_dp/)
    phi = (/0.0_dp, 0.15_dp*pi, 0.0_dp, -0.1_dp*pi/)

    bmod_neo2 = (/0.86660052969026313_dp, &
                  1.0488269770770471_dp, &
                  0.87920244657041735_dp, &
                  1.0499999999555609_dp/)

    call field%neo_field_init(bc_filename, stor(1))
    do case = 1, n_cases
        call field%neo_change_stor(stor(case))
        call field%compute_B_mod(theta(case), phi(case), bmod)
        if (abs(bmod/bmod_neo2(case) - 1) > retol) then
            print *, "-------------------------------------------------------------"
            print *, "test_against_neo2_field failed: B NEO-2"
            print *, "B: ", bmod
            print *, "NEO-2: ", bmod_neo2(case)
            print *, "relative error: ", 1.0_dp - bmod/bmod_neo2(case)
            test_failed = .true.
        end if
    end do

    call plt%initialize(xlabel="$s_{tor}$", ylabel="$B$ [T]", legend=.true.)
    call linspace(0.1_dp, 0.9999_dp, n_stor, stor_plot)
    do case = 1, n_stor
        call field%neo_change_stor(stor_plot(case))
        call field%compute_B_mod(0.0_dp, -0.1_dp*pi, B_at_chi_pi(case))
    end do
    call plt%add_plot(stor_plot, B_at_chi_pi, "$B(\chi = \pi)$", "r-")
    call plt%add_plot((/stor_plot(1), stor_plot(n_stor)/), &
                      (/B_max, B_max/), &
                      "analytic $B_{max}$", &
                      "k-")
    call plt%add_plot((/stor(2), stor(4)/), &
                      (/bmod_neo2(2), bmod_neo2(4)/), &
                      "NEO-2", &
                      "bx", &
                      markersize=6)
    call plt%show()

    if (test_failed) error stop

end program test_against_neo2_anti_sigma
