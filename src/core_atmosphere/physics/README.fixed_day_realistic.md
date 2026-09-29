# Running a fixed day over the real Earth

Instructions for running MPAS-Atmosphere with the insolation held at one day of
the year, but with the real topography, land surface, sea surface temperature
and sea ice. The result is a perpetual-March (or perpetual-July, or any other
day) climate on a realistic planet.

`config_perpetual_julday` handles the radiation and the ozone. The option itself
is described in `README.steady_forcing.md`. Over the real Earth, several other
things follow the calendar date. Some can be frozen with a namelist setting and
some cannot. This document covers those. For an idealized planet with no land,
use `README.aquaplanet.md` instead.

**The main rule: start the run on the calendar date that matches the fixed
day.** The initial SST, sea ice, snow, soil state and vegetation all come from
the start date. Several of them then stay at those values. If the start date
matches the insolation, all of them agree with each other from the first time
step.


## Step 1: choose the day and the matching start date

`config_perpetual_julday` counts days from 1 January 00 UTC, starting at zero.
This is the same count as `curr_julday` in `mpas_atmphys_manager.F`. So day `X`
is the calendar date `X` days after 1 January 00 UTC:

| `config_perpetual_julday` | Meaning | Start date, non-leap year | Start date, leap year |
| --- | --- | --- | --- |
| `80.` | vernal equinox | 22 March 00 UTC | 21 March 00 UTC |
| `172.` | June solstice | 22 June 00 UTC | 21 June 00 UTC |
| `266.` | September equinox (approx.) | 24 September 00 UTC | 23 September 00 UTC |
| `355.` | December solstice | 22 December 00 UTC | 21 December 00 UTC |

The declination and solar constant for each of these days are listed in
`README.steady_forcing.md`. Pick a year that has reanalysis or analysis data for
the initial conditions. Use a non-leap year so that the table above applies
directly.


## Step 2: initial conditions

Run `init_atmosphere` as usual for a real-data case (`config_init_case = 7`),
with the start date from step 1. Nothing special is needed: the topography,
land use, soil types, monthly vegetation and albedo climatologies, SST, sea ice
and snow are all real.

Do not zero out or edit any static fields. The monthly `greenfrac` and
`albedo12m` climatologies must stay intact. The Noah land surface uses their
annual minimum and maximum (`shdmin`, `shdmax`) to scale leaf area, emissivity
and roughness. If all twelve months are set to the same value, every vegetated
cell is treated as being at its seasonal peak.


## Step 3: namelist.atmosphere

```
&nhyd_model
    config_start_time   = '2011-03-22_00:00:00'
    config_run_duration = '365_00:00:00'
/
&physics
    config_physics_suite         = 'mesoscale_reference'

    config_perpetual_julday      = 80.
    config_perpetual_all_physics = true

    config_perpetual_greeness        = true
    config_perpetual_greeness_julday = 80.
    config_sst_update            = false

    config_fixed_co2             = true
    config_co2vmr                = 390.0e-6
/
```

Apart from `config_physics_suite`, these options have `in_defaults="false"` or
are not in the generated namelist. Add them by hand.

What each one does:

* **`config_perpetual_julday = 80.`** fixes the solar declination, the Earth-Sun
  distance and the ozone climatology at day 80. The diurnal cycle is kept.
* **`config_perpetual_all_physics = true`** fixes the day of the year for the
  non-orographic gravity-wave drag and Noah-MP. With the default Noah land
  surface, and with `config_ngw_scheme` off (its default), it has no effect. It
  is still worth setting, so that nothing quietly follows the calendar if you
  later turn the gravity-wave source on or switch to Noah-MP.
* **`config_perpetual_greeness = true`** with
  **`config_perpetual_greeness_julday = 80.`** holds the vegetation fraction and
  the background albedo at day 80. This is the important one over land; see
  below.
* **`config_sst_update = false`** is the default. It keeps SST and sea ice at
  their values in `init.nc`, which come from the start date. Do not turn it on
  unless you have prepared a surface update file with no seasonal cycle: a
  normal SST update file carries the real season.
* **`config_fixed_co2`** only matters with CAM radiation, where CO2 otherwise
  follows a scenario by model year. RRTMG, which every bundled suite uses, is
  already constant at 379 ppmv. Set a value if your protocol needs one; otherwise
  this line can be left out.


### Vegetation and albedo

`config_perpetual_julday` does **not** freeze the green vegetation fraction
(`vegfra`) or the background surface albedo (`sfc_albbck`). Both are
interpolated from the monthly `greenfrac` and `albedo12m` climatologies by
`physics_update_surface`, every `config_greeness_update` (24 hours by default).
Without a separate switch they follow the model date, so vegetation greens up
and browns off through the year even though the sun does not move.

`config_perpetual_greeness` is that switch. It is off by default. When it is on,
`physics_timetracker` passes `physics_update_surface` a fixed date built from
`config_perpetual_greeness_julday` in place of the model date:

* The day is counted from 1 January 00 UTC, the same way as
  `config_perpetual_julday`. Normally you give both the same value.
* Only the integer part is used, because the monthly interpolation works in
  whole days. Values must be in `[0., 365.)`. Switching the option on without
  setting a day, or with a day outside that range, stops the model at startup.
* The date is placed in a fixed non-leap year (2001), so the values are the same
  in every year of the run. A restart does not reset them to the restart date.
  Runs can therefore be split into segments of any length.
* The monthly climatologies are left alone, so `shdmin` and `shdmax` keep their
  real annual range.
* `config_sfc_albedo = false` still turns off the albedo update altogether. In
  that case the switch only affects `vegfra`.

The two days are separate options on purpose. A run can hold the insolation at
one day and the vegetation at another. The vegetation can also be frozen in an
otherwise seasonal run, which is the one use of `config_perpetual_greeness`
without `config_perpetual_julday`.

With the switch on, the start date no longer affects the vegetation, but it
still sets the initial SST, sea ice, snow and soil state. Starting on the
matching calendar date is still the right choice.


## Step 4: what still changes over time

With the settings above, these things still change during the run. None of them
come from the calendar. They are the model responding to fixed forcing:

* **Snow, soil moisture and soil temperature** are prognostic and move toward a
  perpetual-day balance. In perpetual March, Northern Hemisphere snow does not go
  through its normal spring melt. It moves toward whatever balance March
  insolation supports. Soil moisture can take a year or more to settle, so
  budget for a spin-up and do not use the first months for statistics. To hold
  them near their starting values instead, use land nudging (below).
* **The deep soil temperature `tmn`** stays at its initial value unless
  `config_deepsoiltemp_update` is on. If it is on, it relaxes toward running
  daily and annual means of the skin temperature. It always uses the real
  calendar to find the ends of days and years, whatever
  `config_perpetual_all_physics` says, because a frozen date breaks that
  averaging. In a perpetual run these means have no seasonal cycle, so `tmn`
  settles toward the perpetual-day balance. The annual mean is first updated
  at the end of the first calendar year.
* **SST and sea ice** stay fixed with `config_sst_update = false`. The ocean
  does not respond to the atmosphere at all. If you want it to, the slab ocean
  (`README.slab_ocean.md`) can be combined with real land, but that setup has
  not been worked through.

These still follow the calendar and cannot be switched off with a namelist
option:

* **CAM aerosols** are interpolated by month using the real date. Only CAM
  ozone is frozen. This does not apply to RRTMG.
* **Noah-MP urban irrigation** runs only in May–September, using the real month.
  This only matters with Noah-MP and urban physics enabled.
* **External data read during the run.** Lateral boundary conditions (regional
  runs), nudging reference fields (`README.nudging.md`) and any surface update
  stream all carry real dates and the real season. A fixed-day run should be
  global and not nudged toward analyses, unless the reference data have
  themselves been prepared with no seasonal cycle.


### Optional: nudging the land toward its starting state

```
&physics
    config_land_nudging     = true
    config_land_nudging_tau = 864000.
/
```

With `config_land_nudging` on, the model stores the Noah land state at the start
of the run and relaxes it back toward those values every time step, with the
timescale `config_land_nudging_tau` in seconds (864000 s is 10 days). The
nudged fields are the 16 values per cell that Noah carries from step to step:

| Field | What it is | Values per cell |
| --- | --- | --- |
| `skintemp` | skin temperature | 1 |
| `tslb` | soil temperature | 4 |
| `smois` | total soil moisture | 4 |
| `sh2o` | liquid soil water | 4 |
| `canwat` | water on the canopy | 1 |
| `snow` | snow water equivalent | 1 |
| `snowh` | snow depth | 1 |

This stops the land drifting away from the start date, for example the snow
building up without limit on a fixed winter day. It is off by default. Both
options have `in_defaults="false"`, so add them to `&physics` by hand.

How it works:

* **Land cells only** (`landmask = 1`). Over water the skin temperature is the
  SST's, and at sea-ice cells `tslb` is the ice temperature.
* **Every variable is pulled back by the same fraction each step,**
  `1 - exp(-dt/tau)`. Each field becomes a weighted average of its current and
  starting values. This is stable for any timescale, and it keeps the state
  physical without any clipping:
  - soil moisture stays at least as large as the liquid water, and both stay
    positive;
  - the snow density stays between its current and starting values.
* **The starting state is saved in the restart file.** Every segment of a
  restarted run is nudged toward the state at the start of the first segment.
  The restart file must come from a run that already had land nudging on;
  otherwise the stored state is missing and the model nudges toward zero.
* **Only the Noah land model is supported.** The model stops at startup with
  Noah-MP, which has its own snow layers, groundwater and vegetation state that
  would be left free. It also stops if `config_land_nudging_tau` is not
  positive.

Things to be aware of:

* **Water and energy are not conserved.** The nudging adds and removes heat and
  water with nothing to balance it, so the land budgets do not close.
* **Choose the timescale well above one day.** The skin temperature and the top
  soil layer have a strong diurnal cycle, and the stored state is a single
  snapshot at the start time: 00 UTC is night in Europe and Africa but daytime in
  the Americas. With a timescale of hours, the nudging damps the diurnal cycle
  and pulls each region toward its temperature at that local time. With a
  timescale of 10 days or more, the effect on the diurnal cycle is small while
  the slow drift of the deep soil and the snow is still held back.
* **Snow cover is not nudged directly.** `snowc` is recalculated from `snow` by
  Noah, and the snow age (`snotime`) is left alone because it is a counter, not
  a physical state.


## Step 5: check that the day is fixed

Before running for a year, do the five-day test in `README.fixed_day_test.md`
with this namelist. The top-of-atmosphere checks there do not depend on the land
surface, so they apply unchanged. If the start date is day 80, also start a
control run from a date well away from it, as that document says, so the two
runs can be told apart.

Add one check for the land surface. Put `vegfra` and `sfc_albbck` in
`stream_list.atmosphere.diagnostics` and confirm that they do not change between
output times:

```python
import glob, netCDF4, numpy as np

files = sorted(glob.glob("diag.*.nc"))
for name in ("vegfra", "sfc_albbck"):
    a = netCDF4.Dataset(files[0]).variables[name][0]
    b = netCDF4.Dataset(files[-1]).variables[name][0]
    print(name, "max change:", np.abs(b - a).max())
```

Both should be exactly zero. If they change after 24 hours of model time,
`config_perpetual_greeness` is not on, or is not in the `&physics` record.

The first diagnostics file may be written before the first surface update, in
which case it holds the values from `init.nc` for the start date. If the start
date and `config_perpetual_greeness_julday` differ, compare the second file with
the last file instead.


## Known gotchas

| Symptom | Likely cause |
| --- | --- |
| Vegetation and albedo follow the seasons | `config_perpetual_greeness` not set, or not in `&physics` |
| Model stops at startup mentioning `config_perpetual_greeness_julday` | the switch is on but no day was set, or the day is outside `[0., 365.)` |
| Land surface out of balance at the start, then a large drift | start date does not match `config_perpetual_julday` |
| Every vegetated cell has peak leaf area and roughness | all twelve months of `greenfrac` were set equal, collapsing `shdmin` and `shdmax` |
| SST and sea ice follow the seasons | `config_sst_update = true` with a normal surface update file |
| Circulation still has a seasonal cycle | nudging or lateral boundary data from real dates |
| Land nudged toward zero after a restart | the restart file came from a run without `config_land_nudging` |
| Model stops at startup mentioning `config_land_nudging` | the timescale is not positive, or the land model is not `sf_noah` |


## Status

None of this has been run. The list of things that follow the calendar comes from
reading the code: `physics_timetracker` and `physics_update_surface` in the
physics manager, the Noah and Noah-MP drivers, and the CAM radiation interface.
It has not been checked in a simulation.

`config_perpetual_greeness` has not been compiled or run: no Fortran compiler was
available when it was written. The conversion from the day to a calendar date was
checked separately against Python's `datetime` for every day of the year. The
`vegfra` check in step 5 is the way to confirm the rest.

`config_land_nudging` has not been compiled or run either. To check it, run two
short simulations with and without it and compare `tslb`, `smois` and `snow` at
land cells. With nudging on, the difference from the initial state should decay
on the chosen timescale instead of growing.
