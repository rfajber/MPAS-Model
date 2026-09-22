# Held-Suarez forcing in MPAS-Atmosphere

This document describes the Held-Suarez idealized forcing implemented in
`mpas_atm_held_suarez.F`.

Held and Suarez (1994) proposed a benchmark for dry dynamical cores that
replaces the entire physics suite with two analytic terms: Newtonian relaxation
of the temperature field toward a zonally symmetric equilibrium, and Rayleigh
damping of the near-surface winds. Because the forcing is fully specified in a
few lines of algebra, different dynamical cores can be compared against each
other and against the published circulation statistics without any ambiguity
about the physics.


## The forcing

```
d(theta)/dt = -k_t(lat, sigma) * (theta - theta_eq(lat, sigma))
d(u)/dt     = -k_v(sigma) * u
```

with

```
k_t = k_a + (k_s - k_a) * max(0, (sigma - sigma_b)/(1 - sigma_b)) * cos^4(lat)
k_v = k_f *              max(0, (sigma - sigma_b)/(1 - sigma_b))
```

The published equilibrium state is written in temperature,

```
T_eq = max{ T_strat, [315 - dT_y sin^2(lat) - d(theta)_z log(p/p0) cos^2(lat)] (p/p0)^kappa }
```

but MPAS prognoses potential temperature, so the module works in theta
throughout. Dividing by the Exner function cancels the leading factor and leaves

```
theta_eq = max{ T_strat (p0/p)^kappa, 315 - dT_y sin^2(lat) - d(theta)_z log(p/p0) cos^2(lat) }
```

Nothing is approximated by this: it is the same field expressed in a different
variable, and the two agree to round-off. It also happens to be the cheaper form,
since the Exner factor only survives in the stratospheric floor.

MPAS uses `kappa = rgas/cp` with `cp = 7*rgas/2`, so `kappa = 2/7` exactly, which
is the value the benchmark specifies.


## Evaluation against the reference profile

Held and Suarez define the forcing against `sigma = p/p_s`. In a model with a
pressure-based vertical coordinate that has to be re-evaluated every timestep,
because the coordinate surfaces move with the flow.

MPAS uses a fixed height coordinate, so this module instead evaluates the
equilibrium state and both rate profiles **once, at the start of a run**, against
the static reference hydrostatic profile already carried in `pressure_base`. The
three resulting fields are constant in time, and the per-timestep forcing reduces
to a stored multiply-add.

Two consequences are worth being explicit about:

* The benchmark specifies a flat, topography-free sphere. With no terrain the
  reference surface pressure is `p0`, so `sigma = p_ref/p0` is exact for the
  configuration the benchmark actually describes. If you run this forcing on a
  mesh *with* topography it degrades gracefully but is no longer the published
  benchmark, since `sigma` is then measured against a uniform reference rather
  than the local surface pressure.
* The fields derive only from static data, so they are simply rebuilt at the
  start of every run rather than carried in the restart stream. Restarting is
  therefore exact with respect to the forcing, and there is no filter-style
  spin-up transient to worry about.


## Namelist options

All options live in the `held_suarez` record of `namelist.atmosphere`. The
defaults are the published benchmark values; they are exposed so that variants
can be explored, not because they normally need changing.

| Option | Default | Meaning |
| --- | --- | --- |
| `config_held_suarez` | `.false.` | Master switch |
| `config_hs_delta_T_y` | `60.` | Equator-to-pole temperature contrast, K |
| `config_hs_delta_theta_z` | `10.` | Vertical potential temperature contrast, K |
| `config_hs_T_strat` | `200.` | Stratospheric temperature floor, K |
| `config_hs_sigma_b` | `0.7` | Sigma level above which Rayleigh damping vanishes |
| `config_hs_tau_a` | `3456000.` | Free-atmosphere thermal relaxation timescale, s (40 days) |
| `config_hs_tau_s` | `345600.` | Equatorial surface thermal relaxation timescale, s (4 days) |
| `config_hs_tau_f` | `86400.` | Near-surface Rayleigh damping timescale, s (1 day) |

Timescales are in seconds for consistency with the rest of the model, with the
benchmark's day values noted above.


## Running the benchmark

The forcing **replaces** physics rather than supplementing it, so:

```
&physics
    config_physics_suite = 'none'
/
&held_suarez
    config_held_suarez = .true.
/
```

Everything else is a standard dry dynamical-core configuration. A few notes:

* **Initial conditions.** `init_atmosphere` cases 2 and 3 (Jablonowski-Williamson
  baroclinic wave with a perturbation) are a reasonable and conventional start.
  The benchmark statistics are insensitive to the initial state once the flow has
  spun up, since the forcing determines the equilibrium.
* **No topography.** The benchmark is defined on a flat sphere; see above.
* **The run is long.** The standard comparison is a time-mean, zonal-mean
  circulation computed over roughly 1000 days after discarding a spin-up of about
  200 days. The implementation is cheap; the integration is not, and that is where
  the real cost of this benchmark lies.
* **Diagnostics.** The quantities to compare against Held and Suarez (1994) are
  the time-mean zonal-mean zonal wind and temperature, and the eddy momentum and
  heat fluxes. These are computed in post-processing, not by this module.


## Implementation notes

The module lives in `src/core_atmosphere/dynamics/mpas_atm_held_suarez.F` and is
called from `atm_srk3` in `mpas_atm_time_integration.F`, immediately after the
nudging block and before the dynamics substep loop. It writes into the same
tendency arrays that carry physics tendencies into the dynamics
(`tend_ru_physics` and `tend_rtheta_physics`). Those arrays are zeroed even in
builds without physics, which is what makes the forcing safe to use with the
physics suite switched off.

The precomputed fields `hs_theta_eq`, `hs_k_t` and `hs_k_v` live in the
`held_suarez` var_struct, gated on the `held_suarez` package, so nothing is
allocated unless the forcing is switched on. `hs_k_v` is stored on edges, since
the Rayleigh drag acts on the edge-normal velocity; the cell-based rate is
averaged to edges once during setup.

The benchmark is dry. The coupling from a potential temperature tendency to the
prognostic `rho_zz*theta_m` retains the moisture factor for generality, but with
a zero water vapor mixing ratio it reduces to a single term, and no moisture or
mass forcing is applied.

### Known limitations

* The forcing has not been run. The equilibrium state and rate profiles have
  been checked against the published formulation, but no integration has been
  performed and the resulting circulation has not been compared against Held and
  Suarez (1994).
* There is no OpenACC path. The forcing is computed on the host while the
  dynamics tendency arrays live on the device, so its contribution would be
  silently discarded in a GPU build. The module aborts at startup when compiled
  with `MPAS_OPENACC` rather than producing a plausible but unforced simulation;
  adding the data movement would follow the pattern in `mpas_atm_nudging.F`.
* `sigma` is referenced to `p0` rather than to a local surface pressure, which is
  exact for the flat benchmark but not for a run with topography.
