module bmod_spectral
    !> Trigonometric surface representation of |B| on a fixed flux surface.
    !>
    !> Built from |B| values on the equidistant Boozer angular grid, where
    !> tensor-product interpolating splines reproduce the data exactly, so
    !> the representation carries no angular interpolation error between
    !> grid nodes. Evaluation is a full 2D Fourier sum; intended for the
    !> field-line maxima search, where between-node spline error directly
    !> corrupts the maxima-height differences that drive the bootstrap
    !> offset coefficients.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: bmod_spectral_t

    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)

    type :: bmod_spectral_t
        logical :: ready = .false.
        integer :: n_theta = 0
        integer :: n_phi = 0
        integer :: nper = 1
        integer, allocatable :: m_freq(:)
        integer, allocatable :: n_freq(:)
        complex(dp), allocatable :: coeff(:, :)
    contains
        procedure :: build
        procedure :: evaluate
        procedure :: destroy
    end type bmod_spectral_t

contains

    subroutine build(self, values, nper)
        !> values(j, k) = |B| at theta_j = (j-1)*2pi/M, phi_k = (k-1)*2pi/(nper*N)
        !> with M = size(values, 1), N = size(values, 2) distinct grid nodes
        !> (periodic endpoints excluded).
        class(bmod_spectral_t), intent(inout) :: self
        real(dp), intent(in) :: values(:, :)
        integer, intent(in) :: nper

        complex(dp), parameter :: imag = (0.0_dp, 1.0_dp)
        complex(dp), allocatable :: e_theta(:, :), e_phi(:, :)
        integer :: m_nodes, n_nodes, row, col, j, k

        call self%destroy()

        m_nodes = size(values, 1)
        n_nodes = size(values, 2)
        self%n_theta = m_nodes
        self%n_phi = n_nodes
        self%nper = nper

        allocate (self%m_freq(m_nodes), self%n_freq(n_nodes))
        do row = 1, m_nodes
            self%m_freq(row) = signed_frequency(row, m_nodes)
        end do
        do col = 1, n_nodes
            self%n_freq(col) = signed_frequency(col, n_nodes)
        end do

        allocate (e_theta(m_nodes, m_nodes), e_phi(n_nodes, n_nodes))
        do j = 1, m_nodes
            do row = 1, m_nodes
                e_theta(row, j) = exp(-imag*TWOPI*real(self%m_freq(row)*(j - 1), dp) &
                                      /real(m_nodes, dp))
            end do
        end do
        do k = 1, n_nodes
            do col = 1, n_nodes
                e_phi(col, k) = exp(-imag*TWOPI*real(self%n_freq(col)*(k - 1), dp) &
                                    /real(n_nodes, dp))
            end do
        end do

        allocate (self%coeff(m_nodes, n_nodes))
        self%coeff = matmul(e_theta, matmul(cmplx(values, kind=dp), &
                                            transpose(e_phi))) &
                     /real(m_nodes*n_nodes, dp)

        self%ready = .true.
    end subroutine build

    integer function signed_frequency(index, n_nodes)
        integer, intent(in) :: index, n_nodes

        signed_frequency = index - 1
        if (signed_frequency > n_nodes/2) signed_frequency = signed_frequency - n_nodes
    end function signed_frequency

    real(dp) function evaluate(self, theta, phi) result(B)
        class(bmod_spectral_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi

        complex(dp), parameter :: imag = (0.0_dp, 1.0_dp)
        complex(dp) :: e_theta(self%n_theta), e_phi(self%n_phi)
        integer :: row, col

        if (.not. self%ready) error stop "bmod_spectral_t%evaluate: not built"

        do row = 1, self%n_theta
            e_theta(row) = exp(imag*real(self%m_freq(row), dp)*theta)
        end do
        do col = 1, self%n_phi
            e_phi(col) = exp(imag*real(self%n_freq(col)*self%nper, dp)*phi)
        end do

        B = real(sum(e_theta*matmul(self%coeff, e_phi)), dp)
    end function evaluate

    subroutine destroy(self)
        class(bmod_spectral_t), intent(inout) :: self

        if (allocated(self%coeff)) deallocate (self%coeff)
        if (allocated(self%m_freq)) deallocate (self%m_freq)
        if (allocated(self%n_freq)) deallocate (self%n_freq)
        self%ready = .false.
        self%n_theta = 0
        self%n_phi = 0
        self%nper = 1
    end subroutine destroy

end module bmod_spectral
