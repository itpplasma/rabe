<!-- When adding/removing modules or changing interfaces in this directory, update this README. -->

# Fieldline Library (`fieldline_lib`)

Modules for representing field lines on a flux surface, computing integrals along
them, and deriving surface-averaged quantities.

## Modules

### `fieldline_mod` (`fieldline.f90`)

Core type `fieldline_t` representing a single field line between two B maxima.

**Fields set by `make_flock_of_fieldlines`:**

| Field | Description |
|-------|-------------|
| `xi_0` | Fieldline label in [0, 2pi) |
| `theta_0`, `phi_0` | Starting angles on the chi=0 line |
| `iota` | Rotational transform (rational approximation) |
| `iota_p` | Effective rotational transform in omnigenous frame |
| `phi_max(2)` | Toroidal positions of left/right B maxima |
| `B_max(2)` | B values at left/right maxima |
| `eta_b` | 1/B_max_global (bounce parameter) |
| `delta_eta` | 1/B_max_local - eta_b (well depth variation) |
| `integral_lambda_b_over_B_squared` | Boundary layer width integral I_j (up to covariant factor) |
| `I_ref` | Reference I = (n/Σ1/√(I_j))², chosen so that average delta_aspect = 0 |
| `delta_aspect_ratio` | sqrt(I_ref/I_j) - 1 |
| `integral_one_over_B_squared` | Used for surface average normalization |
| `radial_drift` | Integrated radial drift |
| `integral_nabla_s_over_B_squared` | Integrated nabla_s/B^2 |
| `M_pol`, `N_tor`, `nfp` | Mode numbers and field periods |

**Methods:**
- `get_theta(phi)` — returns theta along the fieldline (scalar or array)

### `make_fieldline` (`make_fieldline.f90`)

**Key subroutine:**
```fortran
call make_flock_of_fieldlines(fieldlines, xi_0, iota, field, M_pol, N_tor, nfp, phi_tol [, split_maxima])
```
- `fieldlines`: pre-allocated `fieldline_t(:)` array
- `xi_0`: fieldline labels from `get_labels`
- `iota`: rational approximation of iota from `get_labels`
- `field`: any `field_t` subtype
- `phi_tol`: tolerance for finding B maxima
- `split_maxima` (optional): set to 1 if more than 2 maxima found per period

Workflow: places fieldlines on chi=0 line, finds B maxima per fieldline,
computes all integrals, then derives `eta_b`, `delta_eta`, `I_ref`, and
`delta_aspect_ratio` across the flock.

### `fieldline_labels` (`fieldline_labels.f90`)

**Key subroutine:**
```fortran
call get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfp, xi_0, approx_iota)
```
- Returns allocatable `xi_0` array and rational `approx_iota`
- Grid is constructed so symmetry points xi_0=pi and xi_0=pi-iota_p
  are either on the grid or symmetric between two grid points

**Helper:** `suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol, phi_tol)` —
checks whether the origin is a B minimum (required assumption).

**Fourier analysis:** `fourier_transform_over_label(fieldlines, fieldline_modes)` —
Fourier analysis of radial drift, delta_eta, and delta_aspect_ratio over the label.

**Types:** `modes_t` (cos/sin coefficients + mode numbers), `fieldline_modes_t`
(contains modes for radial_drift, delta_eta, delta_aspect_ratio).

### `fieldline_integrals` (`fieldline_integrals.f90`)

Computes per-fieldline integrals between B maxima:
- `integral_lambda_b_over_B_squared` (boundary layer width)
- `integral_one_over_B_squared` (normalization)
- `radial_drift`
- `integral_nabla_s_over_B_squared`

### `fieldline_integrands` (`fieldline_integrands.f90`)

Pure integrand functions used by `fieldline_integrals`:
- `local_radial_drift(field, theta, phi, eta)`
- `lambda_over_B_squared(field, theta, phi, eta)`
- `B_squared(field, theta, phi)`
- `nabla_s_over_B_squared(field, theta, phi)`

### `surface_average_mod` (`surface_average.f90`)

Type `surface_average_t` with fields: `normalization`, `B_squared`, `lambda_b`, `nabla_s`.

```fortran
call calc_surface_averages(fieldlines, surface_average)
```
Computes surface averages by taking the mean over fieldline label.

## Typical usage

```fortran
use fieldline_labels, only: get_labels
use make_fieldline, only: make_flock_of_fieldlines
use fieldline_mod, only: fieldline_t

call field%fix_to_surface(stor)
call get_labels(max_n, iota, M_pol, N_tor, nfp, xi_0, approx_iota)
allocate(fieldlines(size(xi_0)))
call make_flock_of_fieldlines(fieldlines, xi_0, approx_iota, field, &
                              M_pol, N_tor, nfp, phi_tol)
! Now fieldlines(:)%I_ref, %delta_eta, %radial_drift etc. are available
```
