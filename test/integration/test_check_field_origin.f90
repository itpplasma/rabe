program test_check_field_origin
    use constants, only: dp
    use mock_field, only: mock_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use field_checks, only: suspect_omnigenous_origin_not_minimum

    implicit none

    real(dp), parameter :: N_tor = 1.0_dp, M_pol = 1.0_dp
    real(dp), parameter :: B_0 = 3.0_dp
    real(dp) :: eps
    type(mock_field_t) :: base_field
    type(mock_perturbed_field_t) :: field
    logical :: test_failed

    test_failed = .false.

    ! B = B_0 - cos(theta - phi) + eps*cos(2*theta), origin is minimum
    call base_field%mock_field_init(1.0_dp, 1.0_dp, B_0, -1.0_dp)

    ! eps = 1/5: origin should still be a minimum -> false
    eps = 0.2_dp
    call field%mock_perturbed_field_init(base_field, 2.0_dp, 0.0_dp, eps)
    if (suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_check_field_origin failed: eps=1/5 should not be suspect"
        test_failed = .true.
    end if

    ! eps = 1/2: origin should still be a small maximum bump, still -> false
    eps = 0.5_dp
    call field%mock_perturbed_field_init(base_field, 2.0_dp, 0.0_dp, eps)
    if (suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_check_field_origin failed: eps=1/2 should not be suspect"
        test_failed = .true.
    end if

    ! eps = 1: origin is a significant maximum, violation too strong -> true
    eps = 1.0_dp
    call field%mock_perturbed_field_init(base_field, 2.0_dp, 0.0_dp, eps)
    if (.not. suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_check_field_origin failed: eps=1 should be suspect"
        test_failed = .true.
    end if

    ! B = B_0 + cos(theta - phi) + eps*cos(2*theta), origin is maximum
    call base_field%mock_field_init(1.0_dp, 1.0_dp, B_0, 1.0_dp)

    ! eps = 1/5: origin is already a maximum -> true
    eps = 0.2_dp
    call field%mock_perturbed_field_init(base_field, 2.0_dp, 0.0_dp, eps)
    if (.not. suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_check_field_origin failed: maximum at origin should be suspect"
        test_failed = .true.
    end if

    eps = -0.2_dp
    call field%mock_perturbed_field_init(base_field, 2.0_dp, 0.0_dp, eps)
    if (.not. suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_check_field_origin failed: "
        print *, "small local minimum on top of maximum at origin should be suspect"
        test_failed = .true.
    end if

    if (test_failed) error stop

end program test_check_field_origin
