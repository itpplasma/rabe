program test_against_neo2

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: retol = 1e-12

    call test_against_neo2_field()

contains

    subroutine test_against_neo2_field
        use neo_field, only: neo_field_t

        type(neo_field_t) :: field
        real(dp) :: stor(3), theta(3), phi(3)
        real(dp) :: bmod, sqrtg, dB_dx(3)
        real(dp) :: bmod_neo2(3), sqrtg_neo2(3), dB_dx_neo2(3, 3)
        character(len=*), parameter :: bc_filename = "input/quasi_helical.bc"
        integer :: idx

        stor = (/0.02_dp, 0.50_dp, 0.98_dp/)
        theta = (/1.00_dp, -1.00_dp, 0.00_dp/)
        phi = (/-1.00_dp, 0.00_dp, 1.00_dp/)

     bmod_neo2 = (/5.8461732541782538_dp, 5.9754087905985811_dp, 5.8111625647178791_dp/)
 sqrtg_neo2 = (/-19930779.196453653_dp, -19077980.387761779_dp, -20171657.785802342_dp/)
        dB_dx_neo2(1,:) = (/-0.92227365739728728_dp,      -0.52436170385830649_dp,       0.13109042596457662_dp/)
        dB_dx_neo2(2,:) = (/2.3191833894393170_dp,      0.42134787926653139_dp,      -0.10533696981663285_dp/)
        dB_dx_neo2(3,:) = (/-1.7916431811050300_dp,      -0.45822445110624527_dp,       0.11455611277656132_dp/)

        call field%neo_field_init(bc_filename, stor(1))
        do idx = 1, size(stor)
            call field%compute_B_sqrtg_dB_dx(theta(idx), phi(idx), bmod, sqrtg, dB_dx)
            if (abs(bmod/bmod_neo2(idx) - 1) > retol) then
                print *, "-------------------------------------------------------------"
                print *, "test_against_neo2_field failed: B"
                print *, "B: ", bmod
                print *, "NEO-2: ", bmod_neo2(idx)
                error stop
            end if
            if (abs(sqrtg/sqrtg_neo2(idx) - 1) > retol) then
                print *, "-------------------------------------------------------------"
                print *, "test_against_neo2_field failed: sqrtg"
                print *, "sqrtg: ", sqrtg
                print *, "NEO-2: ", sqrtg_neo2(idx)
                error stop
            end if
            if (any(abs(dB_dx/dB_dx_neo2(idx, :) - 1) > retol)) then
                print *, "-------------------------------------------------------------"
                print *, "test_against_neo2_field failed: dB_dx"
                print *, "dB_dx: ", dB_dx
                print *, "NEO-2: ", dB_dx_neo2(idx, :)
                error stop
            end if
        end do
    end subroutine test_against_neo2_field

end program test_against_neo2
