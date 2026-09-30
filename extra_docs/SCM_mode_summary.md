# Running MPAS-Atmosphere as a single-column model

This note summarizes what would have to change to run MPAS-Atmosphere in a
single-column (SCM) mode. It describes proposed changes only; nothing here has
been implemented.

MPAS has no SCM mode. The practical route, and the one WRF's single-column case
(`em_scm_xy`) uses, is a tiny doubly periodic planar mesh where every column
starts identical. Physics is column-by-column and the dynamics sees no
horizontal gradients, so the columns stay identical and the run behaves as one
column. The missing pieces are:
- initialization from a sounding;
- prescribed large-scale forcing;
- optionally, prescribed surface fluxes.

## 1. What already works (no code changes)

| Piece | Existing support |
|---|---|
| Planar, doubly periodic meshes | Used by the idealized init cases 4, 5, 6 and 10. The spherical curvature terms are inside `#ifdef CURVATURE`, which no build file defines ([mpas_atm_time_integration.F:6603](../src/core_atmosphere/dynamics/mpas_atm_time_integration.F#L6603)). A planar mesh given a real latitude is therefore harmless to the dynamics |
| Columns staying identical | Physics is column-wise with no random numbers (MYNN stochastic perturbations, `spp_pbl`, are off). Horizontal advection and diffusion are exactly zero for uniform fields, and the resolved vertical velocity stays 0 |
| Coriolis on an f-plane | `fEdge`/`fVertex` are read from the init file. Case 10 fixes them at 7.29e-5 |
| Constant geostrophic wind | `perturbation_coriolis` is hard-wired on ([mpas_atm_time_integration.F:6340](../src/core_atmosphere/dynamics/mpas_atm_time_integration.F#L6340)) and applies Coriolis to (u − `u_init`). `u_init`/`v_init` in the init file therefore act as a geostrophic wind profile that is constant in time. The code notes it is correct only for constant f, which holds on an f-plane |
| Relaxing to observed profiles | The existing nudging (u, θ, qv, ρ, aerosols, with a height taper and time filter), driven by reference files in the `lbc_in` format ([README.nudging.md](../src/core_atmosphere/dynamics/README.nudging.md)) |
| Solar geometry and date | Radiation uses `latCell`/`lonCell`, which can be set to the site in the mesh file. `config_perpetual_julday` and the other fixed-day options also apply |
| Surface | Fixed or time-varying SST (`sfc_update`), slab ocean, Noah with optional land nudging |
| Scale-aware physics | GF and MYNN take Δx from `config_len_disp`/meshDensity^¼. It can be set to the grid spacing the SCM is meant to represent, independent of the mesh |

## 2. What would need code changes

### 2.1 Initialization from a sounding (moderate)

**Why the existing cases fall short.**
- **Case 7 can't be reused.** Its static-field step assumes a unit sphere: it
  rescales every distance and area by `sphere_radius`
  ([mpas_init_atm_static.F:291](../src/core_init_atmosphere/mpas_init_atm_static.F#L291)).
- **Case 10 (LES) is close but incomplete.**
  [`init_atm_case_les`](../src/core_init_atmosphere/mpas_init_atm_cases.F#L6222)
  already handles planar geometry and hydrostatic balance, but it:
  - hard-codes an analytic sounding (`atm_get_sounding`);
  - adds random θ perturbations;
  - sets only `landmask`, `lu_index` and `xland` for the surface.

**What an SCM init would need to do:**
- read a sounding (height or pressure, θ, qv, u, v) from a file;
- set every field the physics needs uniformly: land use, soil type, soil
  temperature and moisture, SST, skin temperature, albedo and vegetation, snow,
  sea ice;
- set f from the site latitude, and `u_init`/`v_init` as the geostrophic wind;
- leave out any random perturbations, so the columns start exactly identical.

A no-Fortran alternative is a Python script that builds `init.nc` from a case-10
file plus the sounding. That's feasible but error-prone around the hydrostatic
balance and the derived fields.

### 2.2 Large-scale forcing (the largest piece)

A new tendency term added each step, similar to how nudging is applied:
- **Horizontal advection:** prescribed advective tendencies of θ and qv
  (optionally u, v and condensates).
- **Large-scale vertical velocity**, applied as a tendency −w_ls ∂φ/∂z on θ, qv
  and optionally u, v and hydrometeors. A horizontally uniform periodic domain
  has no mean ascent of its own, since the resolved w stays 0.
- **Geostrophic wind** that varies in time, replacing the constant
  `u_init`/`v_init`.

**Convection closure:** the θ forcing must also be added to `rthdynten`. GF and
new Tiedtke close on `rthdynten`/`rqvdynten` as the large-scale forcing, and in
an SCM the resolved dynamics contributes nothing to them.
- `rqvdynten` is computed from the change in qv over the whole step, so it would
  include qv forcing applied in the dynamics.
- `rthdynten` is computed from the advective flux divergence in
  `atm_compute_dyn_tend`, so it would not include θ forcing
  ([mpas_atm_time_integration.F:2807](../src/core_atmosphere/dynamics/mpas_atm_time_integration.F#L2807)).

### 2.3 Forcing input (moderate)

- New Registry fields and an input stream for time-varying forcing profiles:
  advective tendencies, w_ls, ug/vg, and optionally radiative heating.
- Interpolation in time between records, and from the forcing heights to model
  levels.
- Namelist switches for each forcing term.

The nudging and lateral-boundary input code already has a
read-two-records-and-interpolate pattern that could be copied.

### 2.4 Prescribed surface fluxes (optional; small to moderate)

Many standard cases (BOMEX, RICO and similar) prescribe the sensible and latent
heat fluxes. There is currently no way to do that with the physics PBL:
- **MM5 surface layers:** `sf_monin_obukhov` and `sf_monin_obukhov_rev` already
  contain WRF's `scm_force_flux` plumbing, which holds the fluxes at their input
  values. MPAS hard-wires it to 0
  ([mpas_atmphys_driver_sfclayer.F:33](../src/core_atmosphere/physics/mpas_atmphys_driver_sfclayer.F#L33)),
  and it isn't implemented in `sf_mynn`.
- **LES option:** `config_les_surface = 'specified'` with
  `config_surface_heat_flux`/`config_surface_moisture_flux` does prescribe
  fluxes, but only for the LES turbulence scheme
  ([mpas_atm_dissipation_models.F:1526](../src/core_atmosphere/dynamics/mpas_atm_dissipation_models.F#L1526)),
  not for the physics PBL.
- **With MYNN:** the fluxes would need overriding after the surface layer,
  keeping u* consistent. The LSM would need to be off where fluxes are
  prescribed.

### 2.5 Prescribed radiative heating (optional)

Some cases specify the radiative heating rather than computing it. With
radiation off, it can be supplied through the same forcing term, so it needs
nothing beyond 2.2 and 2.3.

## 3. Tooling (outside the MPAS source)

- **Mesh.** A small doubly periodic planar hex mesh, e.g. from MPAS-Tools
  `planar_hex`, with:
  - `latCell`/`lonCell`, the edge and vertex latitudes and longitudes, and f set
    to the site;
  - `meshDensity` = 1 and a sensible `nominalMinDc`.

  The smallest periodic mesh MPAS's halo exchange and connectivity tolerate has
  not been checked; WRF's SCM uses about 3×3.
- **Scripts** to write the forcing file, the nudging reference files in
  `lbc_in` format, and the sounding, and to extract column 1 from the output.

## 4. Things to watch

- **Exactly uniform start.** The columns stay identical only if the initial
  state is exactly uniform. Case 10's random perturbations must go.
- **Two meanings of Δx.** `config_len_disp` sets the Δx seen by scale-aware
  physics: GF's σ, MYNN's scale awareness and the surface-layer gustiness. It
  should be the grid spacing the SCM stands in for, not the mesh spacing.
- **Existing caveats still apply** in an SCM setting:
  - `config_pbl_interval` changes the PBL time step but not how often it runs;
  - GF's forcing comes from the previous step's tendencies.

  See [MYNN_summary.md](MYNN_summary.md) and
  [Grell_Freitas_summary.md](Grell_Freitas_summary.md).

## 5. Alternatives with no MPAS code changes

- **CCPP-SCM** includes CCPP versions of several of the same schemes.
- **WRF's single-column case** (`em_scm_xy`) uses the WRF versions of much of
  this physics.

Either avoids changing MPAS. Neither uses MPAS's dynamical core or its own
physics driver, which differs in details such as how `rthdynten` reaches GF.

## 6. Possible staging

1. Sounding-driven initialization plus a constant geostrophic wind and fixed
   SST, which gives radiative-convective SCM runs.
2. The forcing module and its input stream.
3. Prescribed surface fluxes.

## Status

This is an analysis based on reading the code. None of it has been implemented
or tested. The claims about planar-mesh behaviour, such as columns staying
identical and the minimum periodic mesh size, would need a short test run to
confirm.
