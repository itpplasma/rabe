module rabe
    use boozer_field,        only: boozer_field_t
    use fieldline_mod,       only: fieldline_t
    use fieldline_labels,    only: get_labels
    use make_fieldline,      only: make_flock_of_fieldlines
    use deviation,           only: calc_deviation
    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients,        only: calc_nu_star_crit, &
                                   calc_finite_boundary_layer_correction
    use shaing_callen_mod,   only: calc_trapped_fraction, &
                                   get_non_omnigenous_remainder
    use error_handling,      only: set_unsafe_mode, &
                                   reset_failed_check_counter, &
                                   did_fail_any_sanity_check
    use git_version,         only: git_hash

    implicit none
end module rabe
