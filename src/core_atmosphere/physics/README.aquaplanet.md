# Running a slab aquaplanet

Step-by-step instructions for running MPAS-Atmosphere as an aquaplanet: a flat,
entirely ocean-covered planet with fixed equinoctial insolation and a prognostic
slab ocean.

The model-side configuration is straightforward and is described in
`README.slab_ocean.md` and `README.steady_forcing.md`. Most of the work is in
preparing the initial conditions, because MPAS is built for real-Earth
simulations and the land has to be removed deliberately.

**Read the ordering constraint in step 1 before starting.** Getting it wrong
produces a file that looks correct and runs, but has mountains hidden in its
vertical coordinate.


## Step 1: static fields, with the land removed

MPAS preprocessing runs in three stages, controlled independently by
`config_static_interp`, `config_vertical_grid` and `config_met_interp`. Run the
static stage on its own first:

```
&nhyd_model
    config_init_case = 7
/
&preproc_stages
    config_static_interp = true
    config_vertical_grid = false
    config_met_interp    = false
/
```

Then edit the resulting static file to make the planet an ocean. The key point:

> **`ter` must be zeroed before the vertical grid is generated.** MPAS uses a
> terrain-following height coordinate, so `zgrid` is built *from* `ter`. Zeroing
> the terrain afterwards leaves a vertical coordinate that still follows
> orography that no longer exists, and nothing in the model will complain.

```python
import netCDF4, numpy as np

nc = netCDF4.Dataset("static.nc", "a")
n  = len(nc.dimensions["nCells"])

iswater = int(nc.variables["iswater_lu"][...])   # 16 for USGS

nc.variables["ter"][:]      = 0.0                # flat, before the vertical grid
nc.variables["landmask"][:] = 0                  # 1 = land, 0 = ocean
nc.variables["ivgtyp"][:]   = iswater            # water land-use category
nc.variables["isltyp"][:]   = 14                 # water soil category
nc.variables["shdmin"][:]   = 0.0
nc.variables["shdmax"][:]   = 0.0
nc.close()
```

`ivgtyp` matters more than it looks: the surface layer scheme reads roughness
length, emissivity and thermal inertia from `LANDUSE.TBL` indexed by this
category, so it has to be the water entry or the ocean will behave like whatever
land type it inherited. Read `iswater_lu` from the file rather than hard-coding
16, since it depends on `mminlu`.


## Step 2: vertical grid and atmospheric state

Now run the remaining stages against the edited static file. `zgrid` is built on
the flat terrain, and the first-guess atmospheric state is interpolated onto it:

```
&preproc_stages
    config_static_interp = false
    config_vertical_grid = true
    config_met_interp    = true
/
```

Starting from real reanalysis for a given date is fine. The initial atmosphere is
not important — the forcing determines the equilibrium, and the land is gone, so
the state adjusts within the spin-up.


## Step 3: surface initial conditions

Set the starting ocean temperature and remove any sea ice:

```python
import netCDF4, numpy as np

nc  = netCDF4.Dataset("init.nc", "a")
lat = nc.variables["latCell"][:]                 # radians

# a smooth equator-to-pole profile; any function of latitude will do
sst = 271.0 + 29.0 * np.exp(-lat**2 / (2 * (26 * np.pi/180)**2))

nc.variables["sst"][:]      = sst
nc.variables["skintemp"][:] = sst
nc.variables["xice"][:]     = 0.0
nc.variables["snow"][:]     = 0.0
nc.variables["snowh"][:]    = 0.0
nc.variables["snowc"][:]    = 0.0
nc.close()
```

The profile above is the one used by the moist Held-Suarez test case
(Thatcher and Jablonowski 2016, their Eq. 6): 29 K equator-to-pole contrast on a
271 K floor. With a slab this is only the starting point — the ocean will move
away from it towards its own energy balance.

`xice = 0` matters. The slab does not handle sea ice and will happily take the
water below freezing, so any ice present at initialisation would be inconsistent
from the first step.


## Step 4: namelist.atmosphere

```
&physics
    config_physics_suite     = 'mesoscale_reference'
    config_lsm_scheme        = 'off'

    config_slab_ocean        = true
    config_slab_depth        = 2.5
    config_slab_albedo       = 0.07

    config_perpetual_equinox = true
    config_o3climatology     = false

    config_sst_update        = false
    config_sstdiurn_update   = false
/
```

None of these options appear in a generated `namelist.atmosphere`, so add them by
hand.

Notes on the choices:

* `config_lsm_scheme = 'off'` because there is no land. The **surface layer
  scheme must stay on** — it produces `hfx` and `lh`, without which the slab has
  no turbulent fluxes and will only respond to radiation.
* `config_o3climatology = false` switches to a fixed vertical ozone profile.
  Leaving it true reintroduces a seasonal cycle through the ozone field, which
  defeats the point of the fixed insolation.
* `config_sst_update` and `config_sstdiurn_update` must both be false. The model
  rejects either in combination with the slab at startup, because both would
  write `sst` behind the slab's back and break its energy budget.


## Step 5: streams.atmosphere

Generate a fresh `streams.atmosphere` rather than reusing an old one, so that
`slab_energy_input` and the surface flux terms appear in the `output` stream.
`output` is a mutable stream, so if you must reuse an existing file, add these by
hand:

```xml
<var name="sst"/>
<var name="slab_energy_input"/>
<var name="swdnb"/>  <var name="swupb"/>
<var name="lwdnb"/>  <var name="lwupb"/>
<var name="hfx"/>    <var name="lh"/>
```

Also make sure the `surface` stream is not supplying SST updates; with
`config_sst_update = false` it will not be read for that purpose, but an
`input_interval` left in place is a common source of confusion.


## Step 6: run, and check conservation first

Before looking at any climate, verify the surface is conserving energy. Between
any two output times, at every ocean cell:

```
slab_energy_input(t2) - slab_energy_input(t1)
    == rho_w * c_w * config_slab_depth * (sst(t2) - sst(t1))
```

with `rho_w = 1000` and `c_w = 4181`. This must hold to round-off.

```python
import netCDF4, numpy as np

nc = netCDF4.Dataset("history....nc")
C  = 1000.0 * 4181.0 * 2.5

e   = nc.variables["slab_energy_input"][:]
sst = nc.variables["sst"][:]

lhs = e[-1]   - e[0]
rhs = C * (sst[-1] - sst[0])
print("max |lhs - rhs| =", np.abs(lhs - rhs).max())
```

A non-trivial residual means something other than the slab is writing `sst`, or a
sign or unit error in the flux sum. Diagnose it before running anything long —
this check is cheap and it is the reason the diagnostic exists.


## Spin-up

With the default 2.5 m depth the surface relaxes radiatively on roughly three
weeks, so the ocean equilibrates in months rather than the decades a realistic
50 m mixed layer would need. Allow several months of simulation before treating
the climate as settled, and expect the global mean surface temperature to drift
substantially from its initial value: nothing constrains it except the top-of-
atmosphere balance.

Two things to watch in the first weeks:

* **Runaway cooling at the poles.** There is no sea ice, no freezing point and no
  ocean heat transport, so the polar ocean can become arbitrarily cold. This is
  expected behaviour for this configuration, not a bug, but it does mean the high
  latitudes are not physically meaningful.
* **Global mean drift.** Without a Q-flux the meridional gradient is set entirely
  by the local energy balance, so a warmer tropics and colder poles than observed
  is normal.


## Known gotchas

| Symptom | Likely cause |
| --- | --- |
| Model aborts at startup mentioning `config_sst_update` | the slab and a prescribed SST source are both enabled; disable the latter |
| Vertical levels follow terrain that should not exist | `ter` was zeroed after the vertical grid stage rather than before |
| Ocean has land-like roughness or albedo | `ivgtyp` was not set to the water category |
| Slab responds only to radiation, never to weather | the surface layer scheme is off, so `hfx` and `lh` are zero |
| Seasonal cycle persists despite fixed insolation | `config_o3climatology` is still true |
| `slab_energy_input` missing from output | reused an old `streams.atmosphere` |


## Status

This configuration has not been run. The individual pieces are documented in
`README.slab_ocean.md` and `README.steady_forcing.md`, and none of them have been
exercised in a simulation either. Treat these instructions as a starting point to
be corrected against the first real run, particularly the list of surface fields
in step 1, which was assembled from the code paths rather than from a working
aquaplanet setup.
