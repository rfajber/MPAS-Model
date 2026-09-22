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

By default the relaxation is applied to the instantaneous departure at every
resolved scale — that is, plain grid nudging. Because the driving state is
smooth relative to the model mesh, nudging toward it that way also damps small
scales the model generates itself.

Setting `config_nudging_filter_tau` makes the scheme scale-selective *in time*
instead, which in practice also selects in space. See **Temporal low-pass
filtering** below. Selectivity in wavenumber directly, in the sense of von
Storch et al. (2000), is not implemented; see *Extending the scheme*.


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
| `config_nudging_filter_tau` | real | `0.` | Timescale of the temporal low-pass filter on the departure, in seconds; non-positive disables it |
| `config_nudging_filter_stages` | integer | `2` | Number of cascaded first-order sections in that filter |

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

To make that scale-selective, add a filter timescale:

```
    config_nudging_filter_tau    = 86400.
    config_nudging_filter_stages = 2
```


## Temporal low-pass filtering

Setting `config_nudging_filter_tau` to a positive value passes the departure
from the driving state through a low-pass filter in time before it is relaxed,
so that only the slowly evolving part of the departure is forced.

This is scale-selective in practice, because atmospheric scales are linked
advectively: an eddy of size `L` turns over in roughly `L/U`. At a typical
`U ~ 10 m/s`, a 1000 km wave evolves over about a day while a 100 km feature
evolves in about three hours. Filtering out the fast part of the departure
therefore leaves the model's own small-scale variability largely free while
still holding the large-scale flow to the driving data.

A filter timescale of one to two days is a reasonable starting point. As a
guide, a two-stage cascade with `config_nudging_filter_tau = 86400.` passes a
10-day signal at about 91% amplitude and attenuates a 3-hour signal to about
0.16% — a separation of roughly 600.

### The filter

Each stage is a first-order recursion, `f <- f + a·(input − f)`, with the output
of one stage feeding the next. The requested timescale is divided across the
cascade, so adding stages sharpens the rolloff and flattens the passband without
changing how far back the filter reaches.

This is an IIR filter, which matters for memory: its entire state is the stage
values themselves. A true boxcar running mean would instead need the whole
window of past samples retained, because the sample leaving the window changes
every step — the incremental update `mean += (new − old)/N` saves the summation
but not the storage. One or two extra stages buys a better-shaped filter far
more cheaply than a boxcar of any useful length, and without a boxcar's
accumulating round-off.

### Why the departure and not both states

The filter is applied to the departure, not to the model and driving states
separately. For a linear filter `LP[x_model] − LP[x_drv] == LP[x_model − x_drv]`,
so the two are mathematically identical, but filtering the departure needs one
set of filter states rather than two.

It also removes a trap. Any causal filter lags, so filtering the model state
in-model while pre-filtering the reanalysis offline with a centred window would
leave the two sides out of step and put a systematic timing error into the
forcing. Filtering the departure means whatever lag the filter has is common to
both terms and cancels exactly, so the result does not depend on the filter's
phase response at all. This is also why there is no reason to pre-filter the
ERA5 data offline.

### Memory and restarts

The filter stores one array per field per stage, over four fields (`u` on edges,
`theta`, `qv` and `rho_zz` on cells). Taking one cell field of
`nVertLevels × 8` bytes as the unit, and using `nEdges ≈ 3·nCells`, that is
6 units per stage, or about 2.6 KB per owned cell per stage at 55 levels in
double precision. For comparison the `lbc` pool already costs roughly 52 units,
about 23 KB per cell, so a two-stage filter adds around 23% to the memory the
nudging feature already uses. Fields are allocated only when the filter is
active.

The filter states are carried in the `restart` stream. On a cold start every
stage is seeded with the current departure, so the filter begins in equilibrium
rather than ramping up from zero; on a restart the states come from the restart
file. If you restart a run in which nudging was previously switched off, the
filter will re-seed and there will be a short transient in the forcing while it
spins up.


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
* The relaxation timescale and `config_nudging_filter_tau` are independent
  knobs, but the nudging cannot respond faster than the slower of the two: a
  filter timescale much longer than the relaxation timescale becomes the
  effective response time, no matter how small the relaxation timescale is.


## Implementation notes

The scheme lives in `src/core_atmosphere/dynamics/mpas_atm_nudging.F`. It is
called from `atm_srk3` in `mpas_atm_time_integration.F`, immediately after the
IAU block and before the dynamics substep loop. `mpas_atm_iau.F` is the closest
structural analogue and was used as the template.

There are two packages, both activated in `atm_setup_packages`. The `nudging`
package follows `config_apply_nudging` and gates the `lbc` fields and the
`lbc_in` stream; `nudging_filter` additionally requires
`config_nudging_filter_tau > 0` and gates the `nudging` var_struct holding the
filter states. Nothing is allocated when the corresponding feature is off, so
enabling nudging without the filter costs no filter memory.

The number of filter stages is a namelist-defined Registry dimension,
`nNudgeFilterStages`, and the filter states are declared with it as their
leading dimension so that each column's cascade is contiguous in memory.

The conversion from a dry potential temperature tendency to a `theta_m` tendency
in this module deliberately differs from the corresponding conversion in
`mpas_atm_iau.F`. Differentiating `rho_zz*theta*(1 + rvord*qv)` directly gives a
final term `rvord*rho_zz*theta*d(qv)/dt`, where IAU uses
`rvord*theta*d(rho_zz*qv)/dt`. The two differ by `rvord*theta*qv*d(rho_zz)/dt`,
which vanishes only when the density field is left alone. Since nudging
increments the density, the difference matters here. See the comment in the
source.

### Extending the scheme

Scale selectivity in space, rather than in time, would mean low-pass filtering
the departure spatially inside `atm_add_tend_nudging`, between the
`mpas_atm_get_bdy_state` calls and the relaxation loops. The same identity used
for the temporal filter applies, so only the departure needs filtering.

An iterated Laplacian smoother over cells is the natural choice on an
unstructured mesh: local, reuses the existing halo exchanges, and works on
variable-resolution and regional meshes, at the cost of a smooth rather than
sharp cutoff. Two practical obstacles are worth knowing before starting.
`u` is the edge-normal component, so smoothing it needs a vector Laplacian or a
reconstruct-and-project round trip rather than the scalar stencil that serves
the cell fields; and `mpas_dmpar_exch_group_add_field` takes a Registry field
name, so the smoothing work arrays must be Registry fields with a new exchange
group in `mpas_atm_halos.F`. A spherical-harmonic truncation would be closer to
the classical formulation but requires global reductions and only works on
global meshes.

### Known limitations

* Scale selectivity is available in time but not in wavenumber. The time-space
  correspondence is statistical, not exact: a terrain-locked mesoscale feature
  is stationary, so it has a long timescale and passes the filter and is still
  nudged, while a genuinely fast large-scale mode such as the diurnal cycle is
  filtered out and is not.
* `scalars_tend` is only zeroed inside `physics_get_tend`, which is compiled out
  in builds without physics. Moisture nudging therefore assumes a
  physics-enabled build. This is the same exposure the IAU scheme already has.
* The OpenACC path moves the `lbc` fields and the filter states to the device
  when nudging is active, but the scheme has not been tested on GPUs.
* The unfiltered scheme has been validated in a full simulation. The temporal
  filter has not: that code path, including the restart of the filter states,
  has never been exercised in a real run. Treat it as less proven than the rest
  when enabling `config_nudging_filter_tau`.
