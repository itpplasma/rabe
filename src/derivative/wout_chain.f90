module wout_chain
    !! Chain of the derivative mode with libneo's VMEC -> Boozer tangent and
    !! adjoint (module `boozer_vmec_ad`), so that the offset coefficients are
    !! differentiated with respect to the VMEC wout harmonics.
    !!
    !! The composed map is wout -> (bmnc, iota, G, I) -> rabe outputs. The first
    !! factor is libneo's discrete VMEC -> Boozer algorithm with its exact JVP
    !! and VJP with respect to `rmnc`, `zmns`, `lmns` (wout shapes, half-mesh
    !! `lmns`) and the full-mesh `iotaf`; the second factor is the dense
    !! spectrum Jacobian of `spectrum_derivative`. `nabla_s_mnc`, `psi_edge` and
    !! `R` are not produced by libneo's transform and are held fixed, so their
    !! columns do not enter the chain. The `iota` column of the spectrum
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
                                   off_iota, off_G, off_I

    implicit none
    private

    public :: wout_chain_t, boozer_wout_tangent_t, wout_tangent_alloc
    public :: wout_chain_init, wout_chain_spectrum, wout_chain_jvp, wout_chain_vjp

    !> Differentiated VMEC -> Boozer transform of one wout file. libneo keeps the
    !! VMEC splines as module state, so only one chain can be active, and no
    !! other VMEC field may be loaded after `wout_chain_init`.
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

    !> Overwrite `nfp`, `xm`, `xn`, `bmnc`, `iota`, `G` and `I` of `spec` by
    !! libneo's transform at `s_tor`. If `spec` already holds a mode set (for
    !! `nabla_s_mnc`), it must be the same.
    subroutine wout_chain_spectrum(chain, s_tor, spec)
        type(wout_chain_t), intent(in) :: chain
        real(dp), intent(in) :: s_tor
        type(surface_spectrum_t), intent(inout) :: spec

        type(boozer_spectrum_t) :: y

        call boozer_ad_primal(chain%ctx, [s_tor], y)
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
        if (.not. allocated(spec%nabla_s_mnc)) then
            allocate (spec%nabla_s_mnc(chain%ctx%mnboz))
            spec%nabla_s_mnc = 0.0_dp
        end if
        spec%iota = y%iota(1)
        spec%G = y%g(1)
        spec%I = y%i(1)
    end subroutine wout_chain_spectrum

    !> Tangent of the outputs for a wout tangent `t`:
    !! d values(:, j) = J_j^T (dbmnc, diota, dG, dI)_j on each surface j.
    subroutine wout_chain_jvp(chain, s_list, jacobians, t, dvalues)
        type(wout_chain_t), intent(in) :: chain
        real(dp), intent(in) :: s_list(:)
        real(dp), intent(in) :: jacobians(:, :, :)
            !! (n_params, n_outputs, surface) from `evaluate` at s_list
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(out) :: dvalues(n_outputs, size(s_list))

        type(boozer_spectrum_t) :: y, yd
        integer :: j, K

        K = chain%ctx%mnboz
        call check_shape(jacobians, K, size(s_list))
        call boozer_ad_jvp(chain%ctx, s_list, t, y, yd)
        do j = 1, size(s_list)
            dvalues(:, j) = matmul(yd%bmnc(:, j), jacobians(1:K, :, j)) &
                            + yd%iota(j)*jacobians(2*K + off_iota, :, j) &
                            + yd%g(j)*jacobians(2*K + off_G, :, j) &
                            + yd%i(j)*jacobians(2*K + off_I, :, j)
        end do
    end subroutine wout_chain_jvp

    !> Adjoint: wout cotangent `t` (overwritten) of the output cotangent
    !! `cotangent(output, surface)`, summed over all surfaces.
    subroutine wout_chain_vjp(chain, s_list, jacobians, cotangent, t)
        type(wout_chain_t), intent(in) :: chain
        real(dp), intent(in) :: s_list(:)
        real(dp), intent(in) :: jacobians(:, :, :)
        real(dp), intent(in) :: cotangent(:, :)
        type(boozer_wout_tangent_t), intent(out) :: t

        type(boozer_spectrum_t) :: yb
        integer :: j, K

        K = chain%ctx%mnboz
        call check_shape(jacobians, K, size(s_list))
        if (any(shape(cotangent) /= [n_outputs, size(s_list)])) &
            error stop "wout_chain_vjp: cotangent must be (n_outputs, surface)"
        call spectrum_alloc(chain%ctx, size(s_list), yb)
        do j = 1, size(s_list)
            yb%bmnc(:, j) = matmul(jacobians(1:K, :, j), cotangent(:, j))
            yb%iota(j) = dot_product(jacobians(2*K + off_iota, :, j), cotangent(:, j))
            yb%g(j) = dot_product(jacobians(2*K + off_G, :, j), cotangent(:, j))
            yb%i(j) = dot_product(jacobians(2*K + off_I, :, j), cotangent(:, j))
        end do
        call boozer_ad_vjp(chain%ctx, s_list, yb, t)
    end subroutine wout_chain_vjp

    subroutine check_shape(jacobians, K, n_surf)
        real(dp), intent(in) :: jacobians(:, :, :)
        integer, intent(in) :: K, n_surf

        if (size(jacobians, 1) /= n_params(K) .or. &
            size(jacobians, 2) /= n_outputs .or. size(jacobians, 3) /= n_surf) &
            error stop "wout_chain: jacobians must be (n_params(K), n_outputs, surface)"
    end subroutine check_shape

end module wout_chain
