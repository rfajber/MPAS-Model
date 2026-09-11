# Nudging in MPAS-Atmosphere

This document describes the nudging scheme implemented in `mpas_atm_nudging.F`,
and the steps required to build the driving data it needs from ERA5 or a
similar reanalysis.

Nudging relaxes the model state toward a time-varying driving state so that the
simulation tracks an observed evolution while the model's own physics and
small-scale dynamics remain free to act. It is typically used for dynamical
downscaling, for producing simulations that can be compared against observations
at specific times, and for long runs that would otherwise drift away from the
observed circulation.


## What the scheme does

For each nudged variable the scheme adds a relaxation term toward the driving
state, evaluated on the uncoupled (physically meaningful) fields:

```
d(u)/dt      = (1/tau_u)     * (u_drv      - u)
d(theta)/dt  = (1/tau_theta) * (theta_drv  - theta)
d(qv)/dt     = (1/tau_qv)    * (qv_drv     - qv)
d(rho_zz)/dt = (1/tau_rho)   * (rho_zz_drv - rho_zz)
```

Each variable has its own relaxation timescale, so the wind can be constrained
strongly while the moisture field is left nearly free, or vice versa. The
resulting increments are then mass-coupled for the dynamics.

Surface pressure is not a prognostic variable in MPAS, so it is nudged
indirectly through the dry density `rho_zz`: constraining the column mass
constrains the surface pressure.

**This is grid nudging, not spectral nudging.** The relaxation is applied at
every resolved scale. Because the driving state is smooth relative to the model
mesh, nudging toward it also damps small scales that the model generates itself.
If you need true scale selectivity — constraining only wavenumbers below a
cutoff, in the sense of von Storch et al. (2000) — that filter is not
implemented here. See *Extending the scheme* below.


## Where it acts

The relaxation term is computed once per model timestep, before the Runge-Kutta
loop in `atm_srk3`, and is added to the same tendency arrays that carry the
physics tendencies into the dynamics (`tend_ru_physics`, `tend_rtheta_physics`,
`tend_rho_physics`, and `scalars_tend`). It is applied over the whole mesh, to
every cell and edge owned by a task.

The driving state is read from the `lbc_in` stream at that stream's input
interval, and is linearly interpolated in time to the current model time between
the two bracketing input times.


## Namelist options

All options live in the `nudging` record of `namelist.atmosphere`.

| Option | Type | Default | Meaning |
| --- | --- | --- | --- |
| `config_apply_nudging` | logical | `.false.` | Master switch for the scheme |
| `config_nudging_tau_u` | real | `21600.` | Relaxation timescale for horizontal velocity, in seconds |
| `config_nudging_tau_theta` | real | `21600.` | Relaxation timescale for potential temperature, in seconds |
| `config_nudging_tau_qv` | real | `21600.` | Relaxation timescale for water vapor mixing ratio, in seconds |
| `config_nudging_tau_rho` | real | `21600.` | Relaxation timescale for dry density, and hence surface pressure, in seconds |
| `config_nudging_zbot` | real | `0.` | Height AGL below which no nudging is applied, in metres |
| `config_nudging_ztop` | real | `0.` | Height AGL above which the full nudging strength is applied, in metres |

A non-positive timescale disables nudging of that variable, so individual fields
can be switched off without rebuilding.

Between `config_nudging_zbot` and `config_nudging_ztop` the nudging strength
ramps linearly from zero to full. When `config_nudging_ztop <=
config_nudging_zbot` the taper is disabled and the full strength is applied at
every level, which is the default. Tapering the nudging out near the ground is
usually worthwhile: relaxing the low-level wind toward a coarse driving state
fights the model's own surface drag and boundary-layer scheme.

A minimal configuration:

```
&nudging
    config_apply_nudging     = .true.
    config_nudging_tau_u     = 21600.
    config_nudging_tau_theta = 21600.
    config_nudging_tau_qv    = -1.
    config_nudging_tau_rho   = -1.
    config_nudging_zbot      = 1000.
    config_nudging_ztop      = 2000.
/
```


## Input stream

Nudging reads its driving data from the **`lbc_in`** stream, the same stream used
to supply lateral boundary conditions to limited-area simulations. The stream is
defined in `streams.atmosphere`:

```xml
<immutable_stream name="lbc_in"
                  type="input"
                  filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="input_interval"
                  input_interval="3:00:00" />
```

Reusing `lbc_in` means the existing `init_atmosphere` machinery for generating
boundary data can be used unchanged to generate nudging data, and that the read,
the two-time-level bookkeeping and the time interpolation are all shared with the
limited-area code in `mpas_atm_boundaries.F`. The cost is that the stream keeps
the name `lbc_in` even in a global run, where nothing about it is lateral or a
boundary.

Two consequences worth noting:

* The driving file series must cover the whole simulation. The first read takes
  the latest time at or before the start time; every subsequent read takes the
  earliest time strictly after the current model time.
* `input_interval` sets how often new driving data is read, and therefore the
  interval over which the model linearly interpolates. It must match the
  interval at which you actually generated the files.

For a **global** nudged run set `config_apply_lbcs = false`. MPAS refuses to
start with `config_apply_lbcs = true` on a mesh that has no boundary cells.

For a **limited-area** run you may set both `config_apply_lbcs = true` and
`config_apply_nudging = true`, in which case the same driving files supply the
boundary relaxation and the interior nudging.


## Preparing driving data from ERA5

The driving files are ordinary MPAS LBC files: they hold the driving state
already interpolated onto your mesh and your model vertical levels. They are
produced by `init_atmosphere` with `config_init_case = 9`, exactly as boundary
data is produced for a regional run, and that case interpolates over the whole
mesh with no boundary masking, so it works on a global mesh.

The pipeline is:

```
ERA5 GRIB  ->  ungrib (WPS)  ->  intermediate files  ->  init_atmosphere
                                                            case 7 -> init.nc     (once)
                                                            case 9 -> lbc.*.nc    (series)
```

### 1. Obtain ERA5

Retrieve ERA5 in GRIB from the Copernicus Climate Data Store, at whatever
interval you intend to nudge on. Six-hourly is a common choice and matches the
typical relaxation timescale; three-hourly gives smoother interpolation at
roughly twice the storage cost.

You need both of the CDS datasets, over the whole period plus at least one time
past the end of your run:

* `reanalysis-era5-pressure-levels` — geopotential, temperature, U and V wind,
  and relative humidity (or specific humidity), on all 37 pressure levels.
* `reanalysis-era5-single-levels` — surface pressure, mean sea level pressure,
  2 m temperature, 2 m dewpoint, 10 m U and V, skin temperature, soil
  temperature and soil moisture for all four layers, land-sea mask, sea-ice
  cover, sea surface temperature, and surface geopotential.

For a global nudged run, request the full globe. For a regional run, request a
box comfortably larger than your mesh.

### 2. Run ungrib

Use the WPS `ungrib` program to convert the GRIB files into the WPS intermediate
format that `init_atmosphere` reads. Link the Vtable that matches what you
downloaded — `Vtable.ERA-interim.pl` is the usual choice for ERA5 on pressure
levels — and set `prefix` in `namelist.wps` to the value you will give
`config_met_prefix`, for example `ERA5`.

Note which humidity variable your Vtable produces. If ungrib writes `RH` you
must set `config_use_spechumd = false`; if it writes `SPECHUMD` you must set
`config_use_spechumd = true`. A mismatch here produces a plausible-looking but
badly wrong moisture field.

The result is a series of files named `ERA5:2010-10-23_00` and so on, one per
time.

### 3. Build the initial condition, once

Run `init_atmosphere` with `config_init_case = 7` in the usual way to produce
`init.nc` for your mesh. This is the standard real-data initialization step and
is unchanged by nudging; it establishes the static fields and, importantly, the
vertical grid `zgrid` that the driving data must be interpolated onto.

### 4. Build the driving file series

Run `init_atmosphere` a second time with `config_init_case = 9`:

```
&nhyd_model
    config_init_case   = 9
    config_start_time  = '2010-10-23_00:00:00'
    config_stop_time   = '2010-11-23_00:00:00'
/
&dimensions
    config_nvertlevels   = 55
    config_nsoillevels   = 4
    config_nfglevels     = 38
    config_nfgsoillevels = 4
/
&data_sources
    config_met_prefix   = 'ERA5'
    config_fg_interval  = 21600
    config_use_spechumd = false
/
```

`config_nfglevels` is the number of levels in the intermediate files: 37 ERA5
pressure levels plus one surface level, so 38. `config_fg_interval` is the
spacing of your intermediate files in seconds, 21600 for six-hourly.

The `preproc_stages` record does not need to be set for case 9. Unlike case 7,
which is driven by `config_static_interp`, `config_vertical_grid` and
`config_met_interp`, case 9 selects its own processing stages internally and
ignores those options.

In `streams.init_atmosphere`, the `input` stream must point at the `init.nc` you
built in step 3, not at the bare mesh file — case 9 reads `zgrid` and the static
fields from it. The `output_interval` of the `lbc` stream must equal
`config_fg_interval`; `init_atmosphere` checks this and aborts if they disagree.

```xml
<immutable_stream name="input"
                  type="input"
                  filename_template="init.nc"
                  input_interval="initial_only" />

<immutable_stream name="lbc"
                  type="output"
                  filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="output_interval"
                  packages="lbcs"
                  output_interval="6:00:00" />
```

This writes one `lbc.*.nc` file per driving time across the whole period.

### 5. Run the nudged simulation

Point the atmosphere model's `lbc_in` stream at those files, set
`input_interval` to match the interval you generated them at, and enable nudging
in `namelist.atmosphere` as shown above.


## Choosing relaxation timescales

The relaxation timescale sets how hard the model is pulled toward the driving
state. Some rules of thumb:

* A timescale of about 6 hours is a common starting point for the wind. It is
  strong enough to hold the large-scale circulation close to the reanalysis over
  a long run, and weak enough that the model still develops its own weather.
* The timescale should be comfortably longer than the model timestep and
  comfortably shorter than the run length. A timescale shorter than the driving
  data interval is rarely sensible, since it pulls the model toward an
  interpolated state more strongly than that state is actually known.
* Nudging the wind alone is the most conservative choice and is often enough to
  control the circulation, since the mass and temperature fields adjust to it.
* Moisture nudging interacts directly with the model's precipitation, and mass
  nudging can excite spurious gravity waves if applied too strongly. Enable
  these deliberately, start with longer timescales than you use for the wind,
  and check that the surface pressure tendency and precipitation fields still
  look reasonable.


## Implementation notes

The scheme lives in `src/core_atmosphere/dynamics/mpas_atm_nudging.F`. It is
called from `atm_srk3` in `mpas_atm_time_integration.F`, immediately after the
IAU block and before the dynamics substep loop. `mpas_atm_iau.F` is the closest
structural analogue and was used as the template.

Fields and the `lbc_in` stream are gated on the `nudging` package, activated
from `config_apply_nudging` in `atm_setup_packages`, so nothing is allocated
when nudging is switched off.

The conversion from a dry potential temperature tendency to a `theta_m` tendency
in this module deliberately differs from the corresponding conversion in
`mpas_atm_iau.F`. Differentiating `rho_zz*theta*(1 + rvord*qv)` directly gives a
final term `rvord*rho_zz*theta*d(qv)/dt`, where IAU uses
`rvord*theta*d(rho_zz*qv)/dt`. The two differ by `rvord*theta*qv*d(rho_zz)/dt`,
which vanishes only when the density field is left alone. Since nudging
increments the density, the difference matters here. See the comment in the
source.

### Extending the scheme

To make the nudging scale-selective, apply a low-pass filter to both the model
state and the driving state before differencing them, inside
`atm_add_tend_nudging` between the `mpas_atm_get_bdy_state` calls and the
relaxation loops. Nothing else in the plumbing needs to change. An iterated
Laplacian smoother over cells is the natural choice on an unstructured mesh: it
is local, reuses the existing halo exchanges, and works on variable-resolution
and regional meshes, at the cost of a smooth rather than sharp cutoff. A true
spherical-harmonic truncation would be closer to the classical formulation but
requires global reductions and only works on global meshes.

### Known limitations

* The scheme is not scale-selective, as described above.
* `scalars_tend` is only zeroed inside `physics_get_tend`, which is compiled out
  in builds without physics. Moisture nudging therefore assumes a
  physics-enabled build. This is the same exposure the IAU scheme already has.
* The OpenACC path moves the `lbc` fields to the device when nudging is active,
  but the scheme has not been tested on GPUs.
