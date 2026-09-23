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


## What still varies in time

Two further sources of time variation are **not** affected by
`config_perpetual_equinox` or `config_fixed_co2`, and both need attention for a
genuinely steady forcing:

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

Neither `config_perpetual_equinox` nor `config_fixed_co2` has been exercised in
a simulation. Both changes are small and confined -- `radconst` for the first,
the CO2 assignment in each radiation scheme for the second -- but neither has
been run.
