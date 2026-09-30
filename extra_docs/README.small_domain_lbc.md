# A small limited-area domain as an SCM substitute

This document describes how to run MPAS-Atmosphere on a very small regional
domain driven by lateral boundary conditions (LBCs) from a reanalysis. The aim
is a column-like model of one site. It uses only existing code and tools: no
source changes are needed.

It is not a single-column model. The large-scale forcing is computed by the
model's own dynamics between boundaries that follow the reanalysis, instead of
being prescribed. For most SCM-style questions, such as how the physics responds
to observed large-scale conditions at a site, it plays the same role. See
**How this differs from an SCM** before relying on it. The alternative, a true
SCM mode that would need code changes, is summarized in
[SCM_mode_summary.md](SCM_mode_summary.md).


## Overview

| Step | Tool | Output |
| --- | --- | --- |
| 1. Cut a small region from a global mesh | MPAS-Limited-Area (external) | `site.static.nc` with `bdyMaskCell` |
| 2. Prepare reanalysis intermediate files | WPS `ungrib` | `ERA5:YYYY-MM-DD_HH`, one per driving time |
| 3. Initial conditions | `init_atmosphere`, case 7 | `site.init.nc` |
| 4. Boundary files | `init_atmosphere`, case 9 | `lbc.YYYY-MM-DD_hh.mm.ss.nc`, one per driving time |
| 5. (Optional) SST updates | `init_atmosphere`, case 8 | `site.sfc_update.nc` |
| 6. Run | `atmosphere_model`, `config_apply_lbcs = true` | history and diagnostics |
| 7. (Optional) Nudge the interior | the same run, `config_apply_nudging = true` | |


## Step 1: choose the resolution and cut the region

### How small can the domain be?

The boundary zone is fixed at **7 rows of cells**
([mpas_atm_boundaries.F:36](../src/core_atmosphere/dynamics/mpas_atm_boundaries.F#L36)):

| `bdyMaskCell` | Zone | What happens there |
| --- | --- | --- |
| 7, 6 (outermost 2 rows) | specified | the state is set from the LBC files |
| 5 to 1 (next 5 rows) | relaxation | the state is relaxed toward the LBC files, strongest at the outside |
| 0 | interior | free: only the dynamics and physics act |

Only cells with `bdyMaskCell = 0` evolve freely. A domain must therefore be at
least 15 cells across to have any interior. About **20–25 cells across** leaves a
small free core of a few to a few dozen cells, which is what you analyze. That is
a few hundred cells in total, so the run is very cheap and fits on one MPI task.

### Which mesh spacing?

The domain width is about 20–25 × Δx, so the resolution and the "column" size go
together:

| Global mesh (quasi-uniform) | Δx | Domain for ~22 cells across | Character |
| --- | --- | --- | --- |
| `x1.2621442` | 15 km | ~330 km | column of a regional model, parameterized convection still appropriate |
| `x1.655362` | 30 km | ~660 km | |
| `x1.163842` | 60 km | ~1300 km | GCM-like cells, but the "column" is now a sizeable region |

A variable-resolution mesh can also be used if it is refined over the site. Pick
Δx for the physics you want to test, then accept the domain size it implies.

### Cut the region

Use the MPAS-Limited-Area tool (`create_region`, from the MPAS-Dev
MPAS-Limited-Area repository). Give it the global **static** file rather than the
bare grid, so the land use, soil, terrain and other static fields come with the
region and do not need to be interpolated again.

A circular region is the simplest. The points file looks like this; check the
tool's own README for the exact syntax of your version:

```
Name: site
Type: circle
Point: 36.6, -97.5
radius: 170000
```

Here the point is (lat, lon) in degrees and the radius is in metres: 170 km is
about 11 cells of 15 km.

```
create_region site.pts x1.2621442.static.nc      # writes site.static.nc
```

Check the size of the free interior before going further:

```python
import netCDF4, numpy as np
d = netCDF4.Dataset('site.static.nc')
mask = d['bdyMaskCell'][:]
print('cells:', mask.size, ' interior (bdyMaskCell == 0):', (mask == 0).sum())
for z in range(1, 8):
    print('  zone', z, (mask == z).sum())
```

If the interior count is only a handful, increase the radius.


## Step 2: reanalysis intermediate files

Run WPS `ungrib` on the reanalysis for the whole period, at the interval you want
the boundaries updated (for example 1 or 6 hours for ERA5), with a prefix such as
`ERA5`. The procedure is the same as for any MPAS real-data run. See
[README.nudging.md](../src/core_atmosphere/dynamics/README.nudging.md),
"Preparing driving data from ERA5", for ERA5 specifics.

Request an area comfortably larger than the domain, so the horizontal
interpolation near the edges has data on all sides.


## Step 3: initial conditions (case 7)

In `namelist.init_atmosphere`:

```
&nhyd_model
    config_init_case  = 7
    config_start_time = '2020-06-01_00:00:00'
    config_stop_time  = '2020-06-01_00:00:00'
/
&data_sources
    config_met_prefix = 'ERA5'
/
&dimensions
    config_nvertlevels   = 55
    config_nfglevels     = 38         ! match the intermediate files
    config_nfgsoillevels = 4
/
&vertical_grid
    config_ztop = 30000.0
    config_blend_bdy_terrain = true   ! blend MPAS terrain toward the reanalysis in the boundary zone
/
&preproc_stages
    config_static_interp     = false  ! static fields came with the region in step 1
    config_native_gwd_static = false
    config_vertical_grid     = true
    config_met_interp        = true
    config_input_sst         = false
/
```

In `streams.init_atmosphere`, set the `input` stream to `site.static.nc` and the
`output` stream to `site.init.nc`.

`config_blend_bdy_terrain` makes the terrain in the boundary zone consistent with
the reanalysis terrain, which reduces spurious flow at the boundaries. On a
domain this small that is a large fraction of the cells, so it is worth enabling.

To take SST and sea ice from a different data set, see
[README.separate_sst_seaice.md](../src/core_init_atmosphere/README.separate_sst_seaice.md).


## Step 4: boundary files (case 9)

Change the namelist to cover the whole run and the driving interval:

```
&nhyd_model
    config_init_case  = 9
    config_start_time = '2020-06-01_00:00:00'
    config_stop_time  = '2020-06-11_00:00:00'
/
&data_sources
    config_met_prefix  = 'ERA5'
    config_fg_interval = 21600          ! seconds; must equal the lbc output_interval
/
```

In `streams.init_atmosphere`, point the `input` stream at **`site.init.nc`**
(case 9 needs `zgrid` and the static fields from it), and set the `lbc` stream:

```xml
<immutable_stream name="lbc"
                  type="output"
                  filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="output_interval"
                  packages="lbcs"
                  output_interval="6:00:00" />
```

`init_atmosphere` aborts if `output_interval` and `config_fg_interval` disagree.
The result is one `lbc.*.nc` file per driving time, from the start to the stop
time inclusive.


## Step 5 (optional): SST updates

For runs longer than a few days, build a surface-update file with case 8 over the
same period, and set `config_sst_update = true` in the model with the `surface`
input stream pointing at it. Remember to set that stream's `input_interval`. See
[README.separate_sst_seaice.md](../src/core_init_atmosphere/README.separate_sst_seaice.md),
"Using the file in the model".


## Step 6: run the model

In `namelist.atmosphere`:

```
&nhyd_model
    config_dt         = 90.0            ! about 6 x dx in km, as for any MPAS run
    config_start_time = '2020-06-01_00:00:00'
    config_run_duration = '10_00:00:00'
/
&limited_area
    config_apply_lbcs = true
/
&physics
    config_physics_suite = 'mesoscale_reference'   ! or any other choice
/
```

In `streams.atmosphere`, the `input` stream is `site.init.nc`, and the `lbc_in`
stream reads the boundary files at the interval they were written:

```xml
<immutable_stream name="lbc_in"
                  type="input"
                  filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="input_interval"
                  packages="limited_area;nudging"
                  input_interval="6:00:00" />
```

The default `input_interval` is 3 hours. It must match the interval of the files
you made.

The model checks at startup that `config_apply_lbcs = true` goes with a mesh that
has boundary cells, and that the `lbc_in` interval is set
([mpas_atm_boundaries.F:845](../src/core_atmosphere/dynamics/mpas_atm_boundaries.F#L845)).

A few hundred cells run comfortably on one MPI task, so no partition
(`graph.info`) file is needed.


## Step 7 (optional): nudge the interior

The branch's nudging reads the same `lbc_in` files, so the interior can also be
relaxed toward the reanalysis
([README.nudging.md](../src/core_atmosphere/dynamics/README.nudging.md)):

```
&nudging
    config_apply_nudging     = true
    config_nudging_tau_u     = 3600.    ! constrain the wind strongly
    config_nudging_tau_theta = -1.      ! leave temperature free
    config_nudging_tau_qv    = -1.      ! leave moisture free
    config_nudging_tau_rho   = -1.
/
```

Nudging the wind but not θ and qv is a common SCM-style choice. The large-scale
circulation follows the reanalysis, while the thermodynamic profiles respond to
the model's physics. A non-positive timescale disables nudging of that variable.
`config_nudging_zbot` can keep nudging out of the boundary layer.


## Analysing the output

Use the interior cells only. The boundary rows are set or relaxed from the
reanalysis and do not reflect the model's physics.

```python
import netCDF4, numpy as np
st = netCDF4.Dataset('site.static.nc')
mask = st['bdyMaskCell'][:]
lat, lon = np.degrees(st['latCell'][:]), np.degrees(st['lonCell'][:])
interior = np.where(mask == 0)[0]

# the interior cell closest to the site
site = (36.6, -97.5)
d2 = (lat[interior] - site[0])**2 + ((lon[interior] - site[1] + 180) % 360 - 180)**2
centre = interior[np.argmin(d2)]
print('interior cells:', interior.size, ' centre cell index:', centre)

h = netCDF4.Dataset('history.2020-06-01_00.00.00.nc')
theta_centre   = h['theta'][:, centre, :]                  # (time, level)
theta_interior = h['theta'][:, interior, :].mean(axis=1)   # interior mean
```

The interior mean is less noisy. The centre cell is the closest analogue of a
single column.


## How this differs from an SCM

- **The large-scale forcing is computed, not prescribed.**
  - Horizontal advection and vertical motion come from the resolved flow between
    the boundaries.
  - In an SCM both are imposed from observations.
  - Here the forcing is consistent with the model's own state, but it cannot be
    set directly or switched off term by term.
- **The columns are not identical.** The interior has horizontal structure, so
  results depend on which cells you average.
- **The boundaries control the interior.** Air crosses a 300 km domain in about
  8 hours at 10 m s⁻¹, so the state is tightly constrained by the reanalysis.
  That is usually desirable for comparison with site observations, but it limits
  free drift. Larger domains or weaker constraints give the physics more room.
- **Domain size and resolution are linked** (step 1). `config_len_disp` changes
  the Δx that scale-aware physics sees (GF, MYNN, surface-layer gustiness), but
  the dynamics also uses it for horizontal diffusion and divergence damping, so
  it is not a free knob.
- **Surface fluxes cannot be prescribed.** The surface is interactive: surface
  layer, land model and SST as usual.


## Things to be aware of

- **Fixed-day options.** The branch's fixed-day options (`config_perpetual_julday`
  and related) are inconsistent with LBCs, which carry real dates and seasons.
  Use real dates for the whole run
  ([README.fixed_day_realistic.md](../src/core_atmosphere/physics/README.fixed_day_realistic.md)).
- **Boundary update interval.** Between boundary files the LBCs are interpolated
  linearly in time. With 6-hourly data, the diurnal cycle at the boundaries is
  poorly sampled; hourly ERA5 is better for a small domain.
- **Terrain.** On a domain this small, a mismatch between MPAS and reanalysis
  terrain affects a large fraction of the cells. Keep
  `config_blend_bdy_terrain = true`, and prefer a site where the terrain is
  smooth.
- **Spin-up.** Allow a few hours, or a day for soil moisture and
  temperature-driven processes, before analysing.


## Known gotchas

| Symptom | Likely cause |
| --- | --- |
| Model stops: boundary cells found but `config_apply_lbcs = false` | the regional mesh needs `config_apply_lbcs = true` |
| Model stops about the `lbc_in` interval | `input_interval` not set, or it does not match the LBC files |
| `init_atmosphere` stops in case 9 about the interval | `config_fg_interval` differs from the `lbc` stream's `output_interval` |
| Case 9 fails reading `zgrid` or static fields | the `input` stream points at the static file instead of `site.init.nc` |
| No interior cells, or only one or two | region too small for the 7 boundary rows. Increase the radius |
| Missing `lbc.*` file partway through the run | the stop time in case 9 was earlier than the end of the run |
| Strong flow or noise near the edges | terrain mismatch; enable `config_blend_bdy_terrain` |


## Status

Not run. The steps use existing `init_atmosphere` cases 7, 8 and 9, the model's
limited-area support in `mpas_atm_boundaries.F`, and the branch's nudging, all as
documented elsewhere. The MPAS-Limited-Area tool is external: check its README
for the exact points-file syntax of your version.
