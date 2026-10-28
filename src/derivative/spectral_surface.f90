module spectral_surface
    !! Flux-surface field given exactly by a finite Boozer Fourier series.
    !!
    !! `B(theta, phi) = sum_k bmnc(k)*cos(xm(k)*theta - xn(k)*phi)` and
    !! `|grad s|(theta, phi) = sum_k nabla_s_mnc(k)*cos(xm(k)*theta - xn(k)*phi)`
    !! with the booz_xform convention that `xn` already contains the number of
    !! field periods. Unlike `fourier_field_t`, nothing is splined, so the field
    !! is an exact (linear) function of its coefficients. This is the field used
    !! by the derivative mode in `spectrum_derivative`, and it can also be fed to
    !! the standard rabe pipeline, because it extends `field_t`.
    !!
    !! The angles passed to the field may be measured from a shifted origin
    !! (`theta_origin`, `phi_origin`): the series is evaluated at
    !! `(theta + theta_origin, phi + phi_origin)`. rabe needs a minimum of the
    !! ideal omnigenous field at the angular origin, and a quasi-isodynamic
    !! spectrum in booz_xform convention has a maximum there, so the derivative
    !! mode moves the origin by half a field period. The coefficients stay
    !! those of the given spectrum, so derivatives with respect to them do not
    !! change.
    use constants, only: dp
    use field_base, only: field_t

    implicit none
    private

    public :: spectral_surface_t, spectral_surface_init

    type, extends(field_t) :: spectral_surface_t
        logical :: initialized = .false.
        integer :: n_modes = 0
        real(dp), allocatable :: xm(:), xn(:)
        real(dp), allocatable :: bmnc(:), nabla_s_mnc(:)
        real(dp) :: B_theta_covariant = 0.0_dp
        real(dp) :: B_phi_covariant = 0.0_dp
        real(dp) :: theta_origin = 0.0_dp
        real(dp) :: phi_origin = 0.0_dp
    contains
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
        procedure :: get_covariant_components
        procedure :: basis
        procedure :: B_derivatives
        procedure :: nabla_s_derivatives
    end type spectral_surface_t

contains

    subroutine spectral_surface_init(self, xm, xn, bmnc, nabla_s_mnc, &
                                     B_theta_covariant, B_phi_covariant)
        type(spectral_surface_t), intent(out) :: self
        real(dp), intent(in) :: xm(:)
            !! poloidal mode numbers
        real(dp), intent(in) :: xn(:)
            !! toroidal mode numbers including the number of field periods
        real(dp), intent(in) :: bmnc(:)
            !! cosine coefficients of B in Tesla
        real(dp), intent(in) :: nabla_s_mnc(:)
            !! cosine coefficients of |grad s| in 1/m
        real(dp), intent(in) :: B_theta_covariant
            !! I, covariant poloidal component of B in T*m
        real(dp), intent(in) :: B_phi_covariant
            !! G, covariant toroidal component of B in T*m

        if (size(xn) /= size(xm) .or. size(bmnc) /= size(xm) .or. &
            size(nabla_s_mnc) /= size(xm)) then
            error stop "spectral_surface_init: mode arrays differ in size"
        end if
        self%n_modes = size(xm)
        self%xm = xm
        self%xn = xn
        self%bmnc = bmnc
        self%nabla_s_mnc = nabla_s_mnc
        self%B_theta_covariant = B_theta_covariant
        self%B_phi_covariant = B_phi_covariant
        self%initialized = .true.
    end subroutine spectral_surface_init

    !> Basis functions cos(xm*theta - xn*phi) and sin(xm*theta - xn*phi).
    subroutine basis(self, theta, phi, cos_k, sin_k)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: cos_k(:), sin_k(:)

        integer :: k
        real(dp) :: arg

        do k = 1, self%n_modes
            arg = self%xm(k)*(theta + self%theta_origin) &
                  - self%xn(k)*(phi + self%phi_origin)
            cos_k(k) = cos(arg)
            sin_k(k) = sin(arg)
        end do
    end subroutine basis

    !> B and its first and second angular derivatives from given basis values.
    pure subroutine B_derivatives(self, cos_k, sin_k, B, B_t, B_p, &
                                  B_tt, B_tp, B_pp)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: cos_k(:), sin_k(:)
        real(dp), intent(out) :: B, B_t, B_p, B_tt, B_tp, B_pp

        B = sum(self%bmnc*cos_k)
        B_t = -sum(self%bmnc*self%xm*sin_k)
        B_p = sum(self%bmnc*self%xn*sin_k)
        B_tt = -sum(self%bmnc*self%xm**2*cos_k)
        B_tp = sum(self%bmnc*self%xm*self%xn*cos_k)
        B_pp = -sum(self%bmnc*self%xn**2*cos_k)
    end subroutine B_derivatives

    !> |grad s| and its first angular derivatives from given basis values.
    pure subroutine nabla_s_derivatives(self, cos_k, sin_k, s, s_t, s_p)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: cos_k(:), sin_k(:)
        real(dp), intent(out) :: s, s_t, s_p

        s = sum(self%nabla_s_mnc*cos_k)
        s_t = -sum(self%nabla_s_mnc*self%xm*sin_k)
        s_p = sum(self%nabla_s_mnc*self%xn*sin_k)
    end subroutine nabla_s_derivatives

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        B_mod = 0.0_dp
        sqrtg = 0.0_dp
        dB_dx = 0.0_dp
        error stop "spectral_surface_t does not provide sqrtg"
    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: cos_k(self%n_modes), sin_k(self%n_modes)

        if (.not. self%initialized) &
            error stop "compute_B_and_dB_dx: spectral_surface_t not initialized"
        call self%basis(theta, phi, cos_k, sin_k)
        B_mod = sum(self%bmnc*cos_k)
        dB_dx(1) = 0.0_dp
        dB_dx(2) = -sum(self%bmnc*self%xm*sin_k)
        dB_dx(3) = sum(self%bmnc*self%xn*sin_k)
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        integer :: k

        if (.not. self%initialized) &
            error stop "compute_B_mod: spectral_surface_t not initialized"
        B_mod = 0.0_dp
        do k = 1, self%n_modes
            B_mod = B_mod + self%bmnc(k)*cos(self%xm(k)*(theta + self%theta_origin) &
                                             - self%xn(k)*(phi + self%phi_origin))
        end do
    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        integer :: k

        if (.not. self%initialized) &
            error stop "compute_nabla_s: spectral_surface_t not initialized"
        nabla_s = 0.0_dp
        do k = 1, self%n_modes
            nabla_s = nabla_s + &
                      self%nabla_s_mnc(k)*cos(self%xm(k)*(theta + self%theta_origin) &
                                              - self%xn(k)*(phi + self%phi_origin))
        end do
    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(spectral_surface_t), intent(in) :: self

        ! The series is summed exactly, so only roundoff remains.
        rel_accuracy_B = 1.0e-12_dp
    end function rel_accuracy_B

    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(spectral_surface_t), intent(in) :: self
        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        if (.not. self%initialized) &
            error stop "get_covariant_components: spectral_surface_t not initialized"
        B_theta_covariant = self%B_theta_covariant
        B_phi_covariant = self%B_phi_covariant
    end subroutine get_covariant_components

end module spectral_surface
