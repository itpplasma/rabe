module bc_field

    use constants, only: dp, pi
    use field_base, only: field_t
    use interpolate, only: BatchSplineData1D, construct_batch_splines_1d, &
                           evaluate_batch_splines_1d, evaluate_batch_splines_1d_der
    use bc_file, only: bc_data_t, read_bc_file

    implicit none
    private

    public :: bc_field_t

    !> Covariant components from the .bc current columns (NEO-2 convention):
    !! B_phi_covariant = -mu0/(2 pi) * nper * (Jpol/nper) [T*m]
    !! B_theta_covariant = -mu0/(2 pi) * Itor [T*m]
    real(dp), parameter :: mu0_over_two_pi = 2.0e-7_dp
    real(dp), parameter :: two_pi = 2.0_dp*pi

    type, extends(field_t) :: bc_field_t
        !>
        !! \brief Magnetic field read directly from a Boozer .bc file.
        !!
        !! \details Evaluates B by direct Fourier sum of the bmnc harmonics,
        !! `B(theta, phi) = sum_k bmnc(k) * cos(m(k)*theta - n(k)*nper*phi)`,
        !! with the harmonics and flux functions splined radially across the
        !! file's surfaces. Only stellarator-symmetric files are supported.
        !<
        logical :: initialized = .false.
        logical :: fixed_to_surface = .false.
        !> Use the legacy NEO-2 v-sign convention (phi_cyl = phi_B + 2pi/nper*v)
        !! of the PoP-2026 .bc data in compute_nabla_s. Default .false. uses
        !! phi_cyl = phi_B - 2pi/nper*v.
        logical :: neo2_nabla_s = .false.
        real(dp) :: fixed_stor
        real(dp) :: nfp
        !> psi_tor_edge = flux/(2 pi) where flux is the header value [T*m^2].
        real(dp) :: psi_tor_edge
        real(dp) :: R
        real(dp) :: a
        integer :: nsurf, nmode
        real(dp), allocatable :: m(:)
        real(dp), allocatable :: n_angle(:) !< file n times nper
        real(dp), allocatable :: s(:)
        !> (nsurf, 3): iota, B_theta_covariant, B_phi_covariant
        real(dp), allocatable :: flux_funcs(:, :)
        real(dp), allocatable :: bmnc(:, :) !< (nsurf, nmode) [T]
        real(dp), allocatable :: geom(:, :) !< (nsurf, 3*nmode): rmnc, zmns, vmns
        type(BatchSplineData1D) :: spl_flux, spl_bmnc, spl_geom
        ! Surface cache filled by fix_to_surface
        real(dp) :: iota_s, B_theta_cov_s, B_phi_cov_s
        real(dp), allocatable :: bmnc_s(:), dbmnc_ds_s(:)
        real(dp), allocatable :: rmnc_s(:), zmns_s(:), vmns_s(:)
    contains
        procedure :: bc_field_init
        procedure :: fix_to_surface
        procedure :: get_iota
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
        procedure :: get_covariant_components
    end type bc_field_t

contains

    !>
    !! \brief Read a Boozer .bc file and prepare radial splines.
    !!
    !! \param[in] bc_filename path to the .bc file
    !! \param[in] neo2_nabla_s opt-in legacy NEO-2 v-sign convention (see type doc)
    !<
    subroutine bc_field_init(self, bc_filename, neo2_nabla_s)
        class(bc_field_t), intent(inout) :: self
        character(len=*), intent(in) :: bc_filename
        logical, intent(in), optional :: neo2_nabla_s

        type(bc_data_t) :: d
        integer :: order

        if (present(neo2_nabla_s)) self%neo2_nabla_s = neo2_nabla_s

        call read_bc_file(bc_filename, d)
        call populate_from_bc_data(self, d)

        if (self%nsurf > 1) then
            call check_uniform_grid(self%s, bc_filename)
            if (self%nsurf < 4) then
                error stop "bc_field_init: need 1 or at least 4 surfaces "// &
                    "for the radial spline"
            end if
            order = min(5, self%nsurf - 1)
            call construct_batch_splines_1d(self%s(1), self%s(self%nsurf), &
                                            self%flux_funcs, order, .false., &
                                            self%spl_flux)
            call construct_batch_splines_1d(self%s(1), self%s(self%nsurf), &
                                            self%bmnc, order, .false., &
                                            self%spl_bmnc)
            call construct_batch_splines_1d(self%s(1), self%s(self%nsurf), &
                                            self%geom, order, .false., &
                                            self%spl_geom)
        end if

        allocate (self%bmnc_s(self%nmode), self%dbmnc_ds_s(self%nmode))
        allocate (self%rmnc_s(self%nmode), self%zmns_s(self%nmode))
        allocate (self%vmns_s(self%nmode))

        self%initialized = .true.
    end subroutine bc_field_init

    !>
    !! \brief Fix the field to the flux surface at normalized toroidal flux stor.
    !!
    !! \details Evaluates the radial splines once and caches the surface
    !! harmonics, their radial derivatives, and the flux functions. Must be
    !! called before any compute_* call or get_covariant_components. For a
    !! single-surface file the stored surface is used and d/ds is zero.
    !<
    subroutine fix_to_surface(self, stor)
        class(bc_field_t), intent(inout) :: self
        real(dp), intent(in) :: stor

        real(dp) :: s_eval, flux_vals(3), geom_vals(3*self%nmode)

        if (.not. self%initialized) &
            error stop "bc_field fix_to_surface: call bc_field_init first"

        if (self%nsurf == 1) then
            self%iota_s = self%flux_funcs(1, 1)
            self%B_theta_cov_s = self%flux_funcs(1, 2)
            self%B_phi_cov_s = self%flux_funcs(1, 3)
            self%bmnc_s = self%bmnc(1, :)
            self%dbmnc_ds_s = 0.0_dp
            self%rmnc_s = self%geom(1, 1:self%nmode)
            self%zmns_s = self%geom(1, self%nmode + 1:2*self%nmode)
            self%vmns_s = self%geom(1, 2*self%nmode + 1:3*self%nmode)
        else
            s_eval = min(max(stor, self%s(1)), self%s(self%nsurf))
            call evaluate_batch_splines_1d(self%spl_flux, s_eval, flux_vals)
            self%iota_s = flux_vals(1)
            self%B_theta_cov_s = flux_vals(2)
            self%B_phi_cov_s = flux_vals(3)
            call evaluate_batch_splines_1d_der(self%spl_bmnc, s_eval, &
                                               self%bmnc_s, self%dbmnc_ds_s)
            call evaluate_batch_splines_1d(self%spl_geom, s_eval, geom_vals)
            self%rmnc_s = geom_vals(1:self%nmode)
            self%zmns_s = geom_vals(self%nmode + 1:2*self%nmode)
            self%vmns_s = geom_vals(2*self%nmode + 1:3*self%nmode)
        end if

        self%fixed_stor = stor
        self%fixed_to_surface = .true.
    end subroutine fix_to_surface

    !>
    !! \brief Rotational transform iota at flux surface stor.
    !<
    subroutine get_iota(self, stor, iota)
        class(bc_field_t), intent(in) :: self
        real(dp), intent(in) :: stor
        real(dp), intent(out) :: iota

        real(dp) :: s_eval, flux_vals(3)

        if (.not. self%initialized) &
            error stop "bc_field get_iota: call bc_field_init first"
        if (self%nsurf == 1) then
            iota = self%flux_funcs(1, 1)
        else
            s_eval = min(max(stor, self%s(1)), self%s(self%nsurf))
            call evaluate_batch_splines_1d(self%spl_flux, s_eval, flux_vals)
            iota = flux_vals(1)
        end if
    end subroutine get_iota

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(bc_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: cosarg(self%nmode)

        if (.not. self%fixed_to_surface) &
            error stop "bc_field compute_B_mod: call fix_to_surface first"
        cosarg = cos(self%m*theta - self%n_angle*phi)
        B_mod = sum(self%bmnc_s*cosarg)
    end subroutine compute_B_mod

    !>
    !! \brief B and its derivatives w.r.t. (s, theta, phi), all SI.
    !<
    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(bc_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: arg(self%nmode), cosarg(self%nmode), sinarg(self%nmode)

        if (.not. self%fixed_to_surface) &
            error stop "bc_field compute_B_and_dB_dx: call fix_to_surface first"
        arg = self%m*theta - self%n_angle*phi
        cosarg = cos(arg)
        sinarg = sin(arg)
        B_mod = sum(self%bmnc_s*cosarg)
        dB_dx(1) = sum(self%dbmnc_ds_s*cosarg)
        dB_dx(2) = -sum(self%m*self%bmnc_s*sinarg)
        dB_dx(3) = sum(self%n_angle*self%bmnc_s*sinarg)
    end subroutine compute_B_and_dB_dx

    !>
    !! \brief B, the Boozer Jacobian sqrt(g), and dB/d(s, theta, phi).
    !!
    !! \details sqrt(g) = -psi_tor_edge*(iota*B_theta + B_phi)/B^2, matching
    !! boozer_field_t (libneo torflux = -psi_tor_edge).
    !<
    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(bc_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        call self%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)
        sqrtg = -self%psi_tor_edge &
                *(self%iota_s*self%B_theta_cov_s + self%B_phi_cov_s)/B_mod**2
    end subroutine compute_B_sqrtg_dB_dx

    !>
    !! \brief |nabla s| from the angular metric of the (R, Z, v) harmonics.
    !!
    !! \details |nabla s| = sqrt(g_tt*g_pp - g_tp^2)*B^2
    !!                      /(|B_phi + iota*B_theta|*|psi_tor_edge|)
    !! with the angular metric built from the cylindrical map
    !! phi_cyl = phi_B - 2pi/nper * v (default) or + (neo2_nabla_s mode).
    !<
    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(bc_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        real(dp) :: arg(self%nmode), cosarg(self%nmode), sinarg(self%nmode)
        real(dp) :: B_mod, R_val, R_t, R_p, Z_t, Z_p, v_t, v_p, p_t, p_p
        real(dp) :: g_tt, g_pp, g_tp, v_factor

        if (.not. self%fixed_to_surface) &
            error stop "bc_field compute_nabla_s: call fix_to_surface first"

        arg = self%m*theta - self%n_angle*phi
        cosarg = cos(arg)
        sinarg = sin(arg)

        B_mod = sum(self%bmnc_s*cosarg)
        R_val = sum(self%rmnc_s*cosarg)
        R_t = -sum(self%m*self%rmnc_s*sinarg)
        R_p = sum(self%n_angle*self%rmnc_s*sinarg)
        Z_t = sum(self%m*self%zmns_s*cosarg)
        Z_p = -sum(self%n_angle*self%zmns_s*cosarg)
        v_t = sum(self%m*self%vmns_s*cosarg)
        v_p = -sum(self%n_angle*self%vmns_s*cosarg)

        v_factor = -two_pi/self%nfp
        if (self%neo2_nabla_s) v_factor = -v_factor
        p_t = v_factor*v_t
        p_p = 1.0_dp + v_factor*v_p

        g_tt = R_t**2 + Z_t**2 + R_val**2*p_t**2
        g_pp = R_p**2 + Z_p**2 + R_val**2*p_p**2
        g_tp = R_t*R_p + Z_t*Z_p + R_val**2*p_t*p_p

        nabla_s = sqrt(abs(g_tt*g_pp - g_tp**2))*B_mod**2 &
                  /(abs(self%iota_s*self%B_theta_cov_s + self%B_phi_cov_s) &
                    *abs(self%psi_tor_edge))
    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(bc_field_t), intent(in) :: self

        rel_accuracy_B = 1e-12_dp
    end function rel_accuracy_B

    !>
    !! \brief Covariant B_theta and B_phi [T*m] on the fixed surface.
    !<
    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(bc_field_t), intent(in) :: self
        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        if (.not. self%fixed_to_surface) &
            error stop "bc_field get_covariant_components: call fix_to_surface first"
        B_theta_covariant = self%B_theta_cov_s
        B_phi_covariant = self%B_phi_cov_s
    end subroutine get_covariant_components

    !> Map the raw libneo bc_data_t to rabe's conventions: currents to
    !! covariant components via -mu0/(2 pi), file n to n*nper, flux to
    !! psi_tor_edge = flux/(2 pi).
    subroutine populate_from_bc_data(self, d)
        class(bc_field_t), intent(inout) :: self
        type(bc_data_t), intent(in) :: d

        self%nsurf = d%nsurf
        self%nmode = d%nmode
        self%nfp = real(d%nper, dp)
        self%psi_tor_edge = d%flux/two_pi
        self%a = d%a
        self%m = real(d%m, dp)
        self%n_angle = real(d%n, dp)*self%nfp
        self%s = d%s
        allocate (self%flux_funcs(d%nsurf, 3))
        self%flux_funcs(:, 1) = d%iota
        self%flux_funcs(:, 2) = -mu0_over_two_pi*d%Itor
        self%flux_funcs(:, 3) = -mu0_over_two_pi*self%nfp*d%Jpol_over_nper
        self%bmnc = d%bmnc
        allocate (self%geom(d%nsurf, 3*d%nmode))
        self%geom(:, 1:d%nmode) = d%rmnc
        self%geom(:, d%nmode + 1:2*d%nmode) = d%zmns
        self%geom(:, 2*d%nmode + 1:3*d%nmode) = d%vmns
        call set_major_radius(self, d)
    end subroutine populate_from_bc_data

    !> R from the (0,0) harmonic of R on the innermost surface;
    !! the header R is the fallback.
    subroutine set_major_radius(self, d)
        class(bc_field_t), intent(inout) :: self
        type(bc_data_t), intent(in) :: d

        integer :: imode

        self%R = d%R
        do imode = 1, d%nmode
            if (d%m(imode) == 0 .and. d%n(imode) == 0) then
                self%R = d%rmnc(1, imode)
                return
            end if
        end do
    end subroutine set_major_radius

    subroutine check_uniform_grid(s, bc_filename)
        real(dp), intent(in) :: s(:)
        character(len=*), intent(in) :: bc_filename

        real(dp) :: h
        integer :: i

        h = (s(size(s)) - s(1))/real(size(s) - 1, dp)
        if (h <= 0.0_dp) error stop "bc_field: s grid must be increasing"
        do i = 2, size(s)
            if (abs(s(i) - s(1) - real(i - 1, dp)*h) > 1e-6_dp*abs(h)) then
                print *, "bc_field: non-uniform s grid in ", trim(bc_filename)
                error stop
            end if
        end do
    end subroutine check_uniform_grid

    subroutine grow(m_tmp, n_tmp, table)
        real(dp), allocatable, intent(inout) :: m_tmp(:), n_tmp(:), table(:, :)

        real(dp), allocatable :: m_new(:), n_new(:), table_new(:, :)
        integer :: n_old

        n_old = size(m_tmp)
        allocate (m_new(2*n_old), n_new(2*n_old), table_new(2*n_old, 4))
        m_new(1:n_old) = m_tmp
        n_new(1:n_old) = n_tmp
        table_new(1:n_old, :) = table
        call move_alloc(m_new, m_tmp)
        call move_alloc(n_new, n_tmp)
        call move_alloc(table_new, table)
    end subroutine grow

end module bc_field
