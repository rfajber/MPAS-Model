# Steady forcing and prescribed surface conditions

This document describes how to run MPAS-Atmosphere with a forcing that does not
vary from year to year: fixed insolation and prescribed, constant sea surface
temperatures. This is the configuration used by aquaplanet experiments and by
idealized climate runs where a statistically steady state is wanted.

Only one of the pieces below needs new code. The rest are existing options or
preprocessing, and are documented here because the combination is not obvious.


## Fixed insolation

```
&physics
    config_perpetual_julday = 80.
/
```

This freezes the date the insolation is computed for. `radconst` evaluates the
solar declination and the Earth-Sun distance factor at
`config_perpetual_julday` instead of at the current model day, so both are
constant for the whole run and there is no seasonal cycle. A non-positive value,
the default, gives the usual seasonal cycle.

Nothing about the orbit is changed: the obliquity is still 23.5 degrees and the
eccentricity factor is still the Paltridge and Platt (1976) expression. Only the
date is held still. Any day of the year can be chosen, so the run can sit at an
equinox, at a solstice, or anywhere between:

| `config_perpetual_julday` | Declination | Solar constant |
| --- | --- | --- |
| `80.` | 0.00 deg | 1380.0 W m^-2 |
| `172.` | +23.50 deg | 1325.2 W m^-2 |
| `266.` | -1.38 deg | 1361.3 W m^-2 |
| `355.` | -23.49 deg | 1416.9 W m^-2 |

Day 80 is the vernal equinox in this formulation and is the natural choice for
an aquaplanet: it is where `sxlong`, and hence the declination, is exactly zero.

Note that it is **not** quite the same as the circular-orbit equinox this option
used to impose. The Earth-Sun distance is now that of 21 March, about 0.7% closer
than the mean, so the solar constant is 1380 rather than `solcon_0 = 1370`. If
you want exactly `solcon_0`, that corresponds to no real day and has to be done
by changing the eccentricity factor rather than the date.

**The diurnal cycle is retained.** Declination and the solar constant are set in
`radconst`, while the diurnal cycle is carried by the hour angle computed inside
the radiation schemes themselves, which this option does not touch. The sun
still rises and sets. This is deliberate: it matches the Aqua-Planet Experiment
specification of Neale and Hoskins (2000).

If you want to remove the diurnal cycle as well, that is a different and larger
change, since the hour angle is computed per column inside the WRF-derived
radiation code (`module_ra_rrtmg_sw.F`) rather than in the MPAS driver.

The option is declared with `in_defaults="false"`, matching the other
radiation options around it, so it does **not** appear in the generated
`namelist.atmosphere`. Add it to the `&physics` record by hand. A value of 367.
or more is rejected at startup, since the formulae take a day of the year.


## Prescribed constant SST

No code is required. `config_sst_update` defaults to `.false.`, which means the
`sst` field read from the initial conditions file is never updated and therefore
stays fixed for the whole run.

So prescribing an SST distribution is a matter of writing the field you want
into `init.nc` before the run, with Python or NCO. Keeping this outside the
model means any profile can be used without rebuilding — an analytic function of
latitude, the APE CONTROL profile, or an observed climatology.

For a full aquaplanet the land surface also has to be removed, which means
setting `landmask` to ocean everywhere, `ter` to zero, and making the remaining
land-surface fields consistent. That is more preprocessing, but still no model
code.


## Slab ocean

Instead of holding the SST fixed, the surface temperature can be made prognostic
by a slab ocean of constant depth driven by the net surface energy flux:

```
&physics
    config_slab_ocean  = true
    config_slab_depth  = 2.5
    config_slab_albedo = 0.07
/
```

This replaces the prescribed SST above rather than supplementing it, and is
incompatible with `config_sst_update`, which the model rejects at startup. The
slab conserves energy exactly, per cell, and provides the diagnostics needed to
verify that.

See `README.slab_ocean.md` for the formulation, the conservation argument, the
diagnostics and the limitations.


## Carbon dioxide

```
&physics
    config_fixed_co2 = true
    config_co2vmr    = 348.0e-6
/
```

`config_co2vmr` is a volume mixing ratio, so the value above is 348 ppmv; set
whatever your protocol calls for. The option is ignored unless
`config_fixed_co2` is true, and the model stops at startup if it is true and the
value is not positive.

This matters more than it first appears, because the two radiation schemes
disagree about CO2 and only one of them is steady:

* **CAM** (`config_radt_lw_scheme = 'cam_lw'`, `cam_sw`) interpolates an SRES A2
  scenario table in the model year, in `camrad`. CO2 therefore drifts upward as
  a run proceeds, spanning 289 ppmv in 1869 to 829 ppmv in 2100. For a long
  integration under nominally steady forcing this is a real trend, not a
  rounding detail. The table is also indexed without bounds checking, so a model
  year outside 1869-2101 reads past the end of the array.
* **RRTMG** (`rrtmg_lw`, `rrtmg_sw`) hardcodes 379 ppmv, the 2005 IPCC value, in
  a `data` statement inside `rrtmg_lwrad` and `rrtmg_swrad`. It is already
  steady, but it was not adjustable before this option.

Setting `config_fixed_co2 = true` replaces both: the scenario interpolation is
bypassed in CAM, and the hardcoded constant is overridden in RRTMG. One value
then applies whichever scheme is selected.

Like the other options here, both are declared `in_defaults="false"` and must be
added to `&physics` by hand. The default when `config_fixed_co2` is false is
unchanged behaviour in each scheme.


## Ozone

`config_perpetual_julday` also freezes the ozone. The monthly ozone climatology
is interpolated to the same locked day as the declination, so with the option
set there is no seasonal cycle in ozone either and `config_o3climatology` can be
left at its default of `.true.`.

That is worth preferring over switching the climatology off. With
`config_o3climatology = .false.` the model falls back to a single fixed vertical
profile, the same at every latitude; locking the date instead keeps the real
latitude-dependent distribution and merely stops it advancing through the year.
It is the radiation's own view of the date, so it follows `config_perpetual_julday`
directly and does not need `config_perpetual_all_physics`.

The lock is applied in two places, because the two radiation schemes reach the
climatology by different routes:

* **RRTMG** interpolates in `o3climatology_from_MPAS`, called once per radiation
  step from `physics_timetracker`. The manager computes a local `radt_julday`
  and passes that instead of `curr_julday`.
* **CAM** interpolates in `oznint`, called inside `radctl`. The locked day is
  passed down through `camrad` as an optional `julian_ozn` argument, so WRF
  callers of the shared file are unaffected and, when it is absent, behaviour is
  unchanged.

The CAM aerosol month interpolation still sees the real date; only the ozone is
locked.


## Freezing the date for the rest of the physics

Three other parameterizations carry an annual cycle of their own, and by default
they keep following the model date even when the radiation has been frozen:

* the **deep soil temperature** update, which relaxes `tmn` toward an annual
  cycle of the surface temperature;
* the **non-orographic gravity-wave drag**, whose source strength is read from
  the `ugwp_limb_tau` table as a function of day of the year;
* the **Noah-MP** land surface, which takes the day of the year directly.

```
&physics
    config_perpetual_julday      = 80.
    config_perpetual_all_physics = true
/
```

`config_perpetual_all_physics` extends the frozen date to all three. It is a
separate switch because freezing the radiation is the common case and freezing
the land surface with it is not always wanted: a run that wants a steady
radiative forcing over a seasonally evolving land surface is a legitimate thing
to ask for. The option is ignored without a positive `config_perpetual_julday`,
and rather than do nothing quietly, that combination is rejected at startup.

For an aquaplanet it makes little practical difference, since there is no land
and the deep soil temperature and Noah-MP never run. The gravity-wave source is
the one that still matters.

Internally the manager computes two days, `radt_julday` for the radiation and
`phys_julday` for these three, and leaves `curr_julday` itself untouched. Nothing
that genuinely needs the model date loses access to it.


## What still varies in time

One source of time variation is **not** affected by any of the options above:

* **Sea ice.** If `xice` is non-zero anywhere and `config_sst_update` is enabled,
  the ice fraction will evolve with the surface update stream. With
  `config_sst_update = .false.` it is held fixed along with the SST.


## A note on the solar constant

MPAS uses `solcon_0 = 1370 W m^-2` (`mpas_atmphys_constants.F`), whereas the
Aqua-Planet Experiment protocol specifies 1365. The 1365 value is present in
that file as a commented-out line. Changing it affects every configuration, not
just steady-forcing runs, so it is deliberately left alone here; adjust it only
if you need exact APE compliance, and be aware of the side effects.


## Status

Neither `config_perpetual_julday` nor `config_fixed_co2` has been exercised in
a simulation. The declination and solar constant in the table above were produced
by compiling `radconst` on its own and calling it, not by running the model.

The ozone lock is the least tested part. The RRTMG path is a single substitution
in `physics_timetracker` and is easy to read. The CAM path adds an argument to
`camrad` and `radctl` in `module_ra_cam.F`, a large WRF-derived file that has not
been compiled here; the argument lists were checked to correspond positionally,
but that is not the same as a build. If you are using RRTMG, which every bundled
physics suite does, the CAM path is not on your critical path.

`config_perpetual_all_physics` has not been run either. It changes which variable
three call sites read, and nothing about what those parameterizations then do.
