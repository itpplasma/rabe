program test_against_neo2
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_magfie, only: neo_magfie_a

    implicit none

    real(dp) :: retol = 1e-12, retol_dB_dx = 1e-11

    call test_against_neo2_field()

contains

    subroutine test_against_neo2_field
        real(dp) :: x(3, 3)
        real(dp) :: bmod, sqrtg, dB_dx(3)
        real(dp) :: bmod_neo2(3), sqrtg_neo2(3), dB_dx_neo2(3, 3)
        integer :: idx

        x(1, :) = (/0.02_dp, 1.00_dp, -1.00_dp/)
        bmod_neo2(1) = 5.8458951223088489_dp
        sqrtg_neo2(1) = -19932756.361788746_dp
        dB_dx_neo2(1,:) = (/-0.93545545807886743_dp, -0.52324687626235267_dp, 0.13081171906558817_dp/)

        x(2, :) = (/0.50_dp, -1.00_dp, 0.00_dp/)
        bmod_neo2(3) = 5.9751856841370872_dp
        sqrtg_neo2(3) = -19079482.277618125_dp
        dB_dx_neo2(3,:) = (/2.3020250573447201_dp,       0.42393972670620234_dp,      -0.10598493167655058_dp/)

        x(3, :) = (/0.98_dp, 0.00_dp, 1.00_dp/)
        bmod_neo2(3) = 5.8109715984855734_dp
        sqrtg_neo2(3) = -20173065.198901489_dp
        dB_dx_neo2(3,:) = (/-1.8031078794845030_dp,      -0.45690311437961667_dp,       0.11422577859490417_dp/)
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
            if (any(abs(dB_dx/dB_dx_neo2(idx, :) - 1) > retol_dB_dx)) then
                print *, "-------------------------------------------------------------"
                print *, "test_against_neo2_field failed: dB_dx"
                print *, "dB_dx: ", dB_dx
                print *, "NEO-2: ", dB_dx_neo2(idx, :)
                error stop
            end if
        end do
    end subroutine test_against_neo2_field

end program test_against_neo2
