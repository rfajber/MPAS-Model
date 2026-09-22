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
    config_perpetual_equinox = true
/
```

This holds the insolation at equinoctial conditions on a circular orbit. The
solar declination is fixed at zero and the Earth-Sun distance correction is
dropped, so the solar constant reduces to `solcon_0`. The result has no seasonal
cycle.

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
`namelist.atmosphere`. Add it to the `&physics` record by hand.


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


## What still varies in time

Two further sources of time variation are **not** affected by
`config_perpetual_equinox`, and both need attention for a genuinely steady
forcing:

* **Ozone.** `config_o3climatology` defaults to `.true.`, which applies a
  monthly-varying ozone climatology. Set it to `.false.` to use a fixed vertical
  profile instead.
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

The `config_perpetual_equinox` option has not been exercised in a simulation.
The change it makes is small and confined to `radconst`, but it has not been run.
