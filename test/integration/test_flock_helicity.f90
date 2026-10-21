program test_flock_helicity
    use constants, only: dp
    use fourier_field, only: fourier_field_t, fourier_field_init
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    implicit none
    type(fourier_field_t) :: field
    type(flock_of_fieldlines_t) :: flock
    real(dp) :: expected

    ! B = 1.5 - 0.3 cos(theta) + 0.3 cos(phi). Along the check line
    ! theta=pi*phi, the QA well |theta|<=pi has its minimum at the origin.
    ! An incorrectly enlarged well includes a deeper, unrelated minimum.
    call fourier_field_init(field, [0, 1, 0], [0, 0, 1], &
                            [1.5_dp, -0.3_dp, 0.3_dp], 0.0_dp, 1.0_dp)
    call make_flock_from_labels(flock, [0.0_dp], 2.0_dp, field, &
                                M_pol=1.0_dp, N_tor=0.0_dp, nfp=1.0_dp)
    ! Along theta=2*phi, dB/dphi=0 at cos(phi)=1/4, both maxima.
    expected = acos(0.25_dp)
    if (maxval(abs(flock%fieldlines(1)%phi_max - [-expected, expected])) &
        > 1.0e-5_dp) error stop 'wrong analytic field-line maxima'
end program test_flock_helicity
