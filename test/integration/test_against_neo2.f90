program test_against_neo2
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_magfie, only: neo_magfie_a

    implicit none

    real(dp) :: retol = 1e-12

    call test_against_neo2_field()

contains

    subroutine test_against_neo2_field
        real(dp) :: x(3, 3)
        real(dp) :: bmod, sqrtg, dB_dx(3)
        real(dp) :: bmod_neo2(3), sqrtg_neo2(3), dB_dx_neo2(3, 3)
        integer :: idx

        x(1, :) = (/0.02_dp, 1.00_dp, -1.00_dp/)
        bmod_neo2(1) = 5.8461732541782538_dp
        sqrtg_neo2(1) = -19930779.196453653_dp
        dB_dx_neo2(1,:) = (/-0.92227365739728728_dp,      -0.52436170385830649_dp,       0.13109042596457662_dp/)

        x(2, :) = (/0.50_dp, -1.00_dp, 0.00_dp/)
        bmod_neo2(3) = 5.9754087905985811_dp
        sqrtg_neo2(3) = -19077980.387761779_dp
        dB_dx_neo2(3,:) = (/2.3191833894393170_dp,      0.42134787926653139_dp,      -0.10533696981663285_dp/)

        x(3, :) = (/0.98_dp, 0.00_dp, 1.00_dp/)
        bmod_neo2(3) = 5.8111625647178791_dp
        sqrtg_neo2(3) = -20171657.785802342_dp
        dB_dx_neo2(3,:) = (/-1.7916431811050300_dp,      -0.45822445110624527_dp,       0.11455611277656132_dp/)
        do idx = 1, size(x, 1)
            call neo_magfie_a(x(idx, :), bmod, sqrtg, dB_dx)
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
