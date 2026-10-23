module wout_chain
    !! Chain of the derivative mode with libneo's VMEC -> Boozer tangent and
    !! adjoint (module `boozer_vmec_ad`), so that the offset coefficients are
    !! differentiated with respect to the VMEC wout harmonics.
    !!
    !! The composed map is wout -> (bmnc, iota, G, I) -> rabe outputs. The first
    !! factor is libneo's discrete VMEC -> Boozer algorithm with its exact JVP
    !! and VJP with respect to `rmnc`, `zmns`, `lmns` (wout shapes, half-mesh
    !! `lmns`), the full-mesh `iotaf` and the scalars `phiedge` (scales the
    !! whole toroidal flux profile) and `rmajor` (`Rmajor_p`); the second
    !! factor is the dense spectrum Jacobian of `spectrum_derivative`. libneo
    !! provides every spectrum parameter: `bmnc`, `nabla_s_mnc` (its |grad s|
    !! spectrum in 1/m), `iota`, `G`, `I`, `psi_edge` = phi_edge/(2 pi) and
    !! `R` = `Rmajor_p`, so all columns of the Jacobian enter the chain. `R` is
    !! read from the file, as rabe's field does, and depends only on `rmajor`.
    !! The `iota` column of the spectrum
    !! Jacobian is the frozen-rational surrogate (see README), and chaining it
    !! with libneo's iota tangent makes the field-line iota move with the true
    !! iota at fixed `p/q` and labels.
    !!
    !! Mode convention: rabe's `mboz` is the largest poloidal mode number
    !! (`m = 0..mboz`), booz_xform and libneo count modes (`m = 0..mboz-1`), so
    !! libneo is initialised with `mboz + 1`. Both use the same ordering.
    use constants, only: dp
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_wout_tangent_t, &
                              boozer_spectrum_t, boozer_ad_init, &
                              boozer_ad_primal, boozer_ad_jvp, boozer_ad_vjp, &
                              wout_tangent_alloc, spectrum_alloc
    use spectrum_derivative, only: surface_spectrum_t, n_outputs, n_params, &
                                   off_iota, off_G, off_I, off_psi_edge, off_R, &
                                   off_dbmnc_ds, off_diota_ds, off_dG_ds, off_dI_ds

    implicit none
    private

    public :: wout_chain_t, boozer_wout_tangent_t, wout_tangent_alloc
    public :: wout_chain_init, wout_chain_spectrum, wout_chain_jvp, wout_chain_vjp

    !> Differentiated VMEC -> Boozer transform of one wout file. The context
    !! holds the wout data and does not use libneo's VMEC spline state, but
    !! `boozer_ad_init` sets libneo's module settings (file name, spline
    !! orders, grid refinement), so a rabe field loaded afterwards uses them.
    type :: wout_chain_t
        type(boozer_ad_context_t) :: ctx
    end type wout_chain_t

contains

    subroutine wout_chain_init(chain, wout_file, mboz, nboz, grid_refinement)
        type(wout_chain_t), intent(out) :: chain
        character(len=*), intent(in) :: wout_file
        integer, intent(in) :: mboz, nboz
            !! rabe convention: m = 0..mboz, |n| <= nboz field-period harmonics
        integer, intent(in) :: grid_refinement

        call boozer_ad_init(chain%ctx, wout_file, mboz + 1, nboz, &
                            grid_refinement=grid_refinement)
    end subroutine wout_chain_init

    !> Set all fields of `spec` (`nfp`, `xm`, `xn`, `bmnc`, `nabla_s_mnc`,
    !! `iota`, `G`, `I`, `psi_edge`, `R`) from libneo's transform at `s_tor`.
    !! If `spec` already holds a mode set, it must be the same. With
    !! `ds_radial > 0` the radial block (`dbmnc_ds`, `diota_ds`, `dG_ds`,
    !! `dI_ds`) is the central difference of the transforms at `s_tor +-
    !! ds_radial` (alpha_criteria.md, Gamma_s).
    subroutine wout_chain_spectrum(chain, s_tor, spec, ds_radial)
        type(wout_chain_t), intent(in) :: chain
        real(dp), intent(in) :: s_tor
        type(surface_spectrum_t), intent(inout) :: spec
        real(dp), intent(in), optional :: ds_radial

        type(boozer_spectrum_t) :: y
        real(dp) :: h

        h = radial_step(ds_radial)
        call boozer_ad_primal(chain%ctx, surface_list([s_tor], h), y)
        if (allocated(spec%xm)) then
            if (size(spec%xm) /= chain%ctx%mnboz) &
                error stop "wout_chain_spectrum: mode count differs"
            if (any(nint(spec%xm) /= chain%ctx%xm_b) .or. &
                any(nint(spec%xn) /= chain%ctx%xn_b)) &
                error stop "wout_chain_spectrum: mode set differs"
        end if
        spec%nfp = real(chain%ctx%nper, dp)
        spec%xm = real(chain%ctx%xm_b, dp)
        spec%xn = real(chain%ctx%xn_b, dp)
        spec%bmnc = y%bmnc(:, 1)
        spec%nabla_s_mnc = y%gradsmnc(:, 1)
        spec%iota = y%iota(1)
        spec%G = y%g(1)
        spec%I = y%i(1)
        spec%psi_edge = y%psi_edge(1)
        spec%R = y%rmajor(1)
        spec%has_radial = h > 0.0_dp
        if (spec%has_radial) then
            spec%dbmnc_ds = (y%bmnc(:, 3) - y%bmnc(:, 2))/(2.0_dp*h)
            spec%diota_ds = (y%iota(3) - y%iota(2))/(2.0_dp*h)
            spec%dG_ds = (y%g(3) - y%g(2))/(2.0_dp*h)
            spec%dI_ds = (y%i(3) - y%i(2))/(2.0_dp*h)
        end if
    end subroutine wout_chain_spectrum

    !> Tangent of the outputs for a wout tangent `t`:
    !! d values(:, j) = J_j^T dparams_j on each surface j, with dparams the
    !! tangent of (bmnc, nabla_s_mnc, iota, G, I, psi_edge, R) and, with
    !! `ds_radial > 0`, of the central differences of the radial block.
    subroutine wout_chain_jvp(chain, s_list, jacobians, t, dvalues, ds_radial)
        type(wout_chain_t), intent(in) :: chain
        real(dp), intent(in) :: s_list(:)
        real(dp), intent(in) :: jacobians(:, :, :)
            !! (n_params, n_outputs, surface) from `evaluate` at s_list
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(out) :: dvalues(n_outputs, size(s_list))
        real(dp), intent(in), optional :: ds_radial

        type(boozer_spectrum_t) :: y, yd
        integer :: j, K, n, jp, jm
        real(dp) :: h, c

        K = chain%ctx%mnboz
        n = size(s_list)
        h = radial_step(ds_radial)
        call check_shape(jacobians, K, n)
        call boozer_ad_jvp(chain%ctx, surface_list(s_list, h), t, y, yd)
        do j = 1, n
            dvalues(:, j) = matmul(yd%bmnc(:, j), jacobians(1:K, :, j)) &
                            + matmul(yd%gradsmnc(:, j), jacobians(K + 1:2*K, :, j)) &
                            + yd%iota(j)*jacobians(2*K + off_iota, :, j) &
                            + yd%g(j)*jacobians(2*K + off_G, :, j) &
                            + yd%i(j)*jacobians(2*K + off_I, :, j) &
                            + yd%psi_edge(j)*jacobians(2*K + off_psi_edge, :, j) &
                            + yd%rmajor(j)*jacobians(2*K + off_R, :, j)
            if (h > 0.0_dp) then
                jm = n + 2*j - 1
                jp = n + 2*j
                c = 0.5_dp/h
                dvalues(:, j) = dvalues(:, j) &
                    + c*matmul(yd%bmnc(:, jp) - yd%bmnc(:, jm), &
                               jacobians(2*K + off_dbmnc_ds + 1:3*K + off_dbmnc_ds, :, j)) &
                    + c*(yd%iota(jp) - yd%iota(jm))*jacobians(3*K + off_diota_ds, :, j) &
                    + c*(yd%g(jp) - yd%g(jm))*jacobians(3*K + off_dG_ds, :, j) &
                    + c*(yd%i(jp) - yd%i(jm))*jacobians(3*K + off_dI_ds, :, j)
            end if
        end do
    end subroutine wout_chain_jvp

    !> Adjoint of the chain, the wout cotangent `t` (overwritten) of the
    !! output cotangent `cotangent(output, surface)`, summed over all surfaces.
    subroutine wout_chain_vjp(chain, s_list, jacobians, cotangent, t, ds_radial)
        type(wout_chain_t), intent(in) :: chain
        real(dp), intent(in) :: s_list(:)
        real(dp), intent(in) :: jacobians(:, :, :)
        real(dp), intent(in) :: cotangent(:, :)
        type(boozer_wout_tangent_t), intent(out) :: t
        real(dp), intent(in), optional :: ds_radial

        type(boozer_spectrum_t) :: yb
        integer :: j, K, n, jp, jm
        real(dp) :: h, c
        real(dp), allocatable :: s_all(:)

        K = chain%ctx%mnboz
        n = size(s_list)
        h = radial_step(ds_radial)
        call check_shape(jacobians, K, n)
        if (any(shape(cotangent) /= [n_outputs, n])) &
            error stop "wout_chain_vjp: cotangent must be (n_outputs, surface)"
        s_all = surface_list(s_list, h)
        call spectrum_alloc(chain%ctx, size(s_all), yb)
        do j = 1, n
            yb%bmnc(:, j) = matmul(jacobians(1:K, :, j), cotangent(:, j))
            yb%gradsmnc(:, j) = matmul(jacobians(K + 1:2*K, :, j), cotangent(:, j))
            yb%iota(j) = dot_product(jacobians(2*K + off_iota, :, j), cotangent(:, j))
            yb%g(j) = dot_product(jacobians(2*K + off_G, :, j), cotangent(:, j))
            yb%i(j) = dot_product(jacobians(2*K + off_I, :, j), cotangent(:, j))
            yb%psi_edge(j) = dot_product(jacobians(2*K + off_psi_edge, :, j), &
                                         cotangent(:, j))
            yb%rmajor(j) = dot_product(jacobians(2*K + off_R, :, j), cotangent(:, j))
            if (h > 0.0_dp) then
                jm = n + 2*j - 1
                jp = n + 2*j
                c = 0.5_dp/h
                yb%bmnc(:, jp) = c*matmul(jacobians(2*K + off_dbmnc_ds + 1: &
                                                    3*K + off_dbmnc_ds, :, j), &
                                          cotangent(:, j))
                yb%iota(jp) = c*dot_product(jacobians(3*K + off_diota_ds, :, j), &
                                            cotangent(:, j))
                yb%g(jp) = c*dot_product(jacobians(3*K + off_dG_ds, :, j), cotangent(:, j))
                yb%i(jp) = c*dot_product(jacobians(3*K + off_dI_ds, :, j), cotangent(:, j))
                yb%bmnc(:, jm) = -yb%bmnc(:, jp)
                yb%iota(jm) = -yb%iota(jp)
                yb%g(jm) = -yb%g(jp)
                yb%i(jm) = -yb%i(jp)
            end if
        end do
        call boozer_ad_vjp(chain%ctx, s_all, yb, t)
    end subroutine wout_chain_vjp

    pure real(dp) function radial_step(ds_radial)
        real(dp), intent(in), optional :: ds_radial

        radial_step = 0.0_dp
        if (present(ds_radial)) radial_step = max(ds_radial, 0.0_dp)
    end function radial_step

    !> The surfaces of the transform: s_list, then s_j - h, s_j + h per surface
    !! when h > 0.
    pure function surface_list(s_list, h) result(s_all)
        real(dp), intent(in) :: s_list(:), h
        real(dp), allocatable :: s_all(:)

        integer :: j, n

        n = size(s_list)
        if (h > 0.0_dp) then
            allocate (s_all(3*n))
            s_all(1:n) = s_list
            do j = 1, n
                s_all(n + 2*j - 1) = s_list(j) - h
                s_all(n + 2*j) = s_list(j) + h
            end do
        else
            s_all = s_list
        end if
    end function surface_list

    subroutine check_shape(jacobians, K, n_surf)
        real(dp), intent(in) :: jacobians(:, :, :)
        integer, intent(in) :: K, n_surf

        if (size(jacobians, 1) /= n_params(K) .or. &
            size(jacobians, 2) /= n_outputs .or. size(jacobians, 3) /= n_surf) &
            error stop "wout_chain: jacobians must be (n_params(K), n_outputs, surface)"
    end subroutine check_shape

end module wout_chain
