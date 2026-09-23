# Slab ocean in MPAS-Atmosphere

This document describes the slab ocean implemented in `physics_update_slab_ocean`
in `mpas_atmphys_update_surface.F`.

The slab makes the sea surface temperature prognostic instead of prescribed. A
layer of water of fixed depth sits under each ocean cell and its temperature
responds to the net energy flux through the surface, so the ocean surface is free
to adjust to the atmosphere rather than forcing it one way.

It is deliberately the simplest thing that conserves energy: no ocean heat
transport, no sea ice, no freezing, and no horizontal communication between
columns. Every ocean cell is an independent heat reservoir.


## What it does

At every cell with `landmask == 0`:

```
rho_w * c_w * h * d(sst)/dt = swdnb - swupb + lwdnb - lwupb - hfx - lh
```

using `rho_w = 1000 kg m^-3` and `c_w = 4181 J kg^-1 K^-1`, with `h` set by
`config_slab_depth`. The right-hand side is the net *downward* energy flux;
`hfx` and `lh` are defined positive upward and so enter negatively.

After stepping the temperature the routine also sets `skintemp` to the new `sst`,
since that is how the surface schemes see the ocean, and pins `sfc_albedo` to
`config_slab_albedo` so that `swupb` stays consistent with `swdnb`.


## Namelist options

All options live in the `physics` record of `namelist.atmosphere`. None appear in
the generated `namelist.atmosphere`, so they must be added by hand.

| Option | Default | Meaning |
| --- | --- | --- |
| `config_slab_ocean` | `.false.` | Master switch |
| `config_slab_depth` | `2.5` | Mixed layer depth, m |
| `config_slab_albedo` | `0.07` | Fixed surface albedo at slab points |

```
&physics
    config_slab_ocean  = true
    config_slab_depth  = 2.5
    config_slab_albedo = 0.07
/
```


## Energy conservation

This is the property the implementation is organised around, so it is worth
stating exactly what is guaranteed and why.

The slab is updated from `physics_timetracker`, which runs at the **top** of each
time step, before the physics drivers. Every flux it reads is therefore the value
the atmosphere and the surface layer actually exchanged over the *preceding*
step. Debiting the slab by exactly those fluxes over exactly that `dt` means
whatever the atmosphere gained from the surface last step, the slab loses this
step. Conservation is exact to round-off, and it is exact per cell, not merely in
the global mean.

Two consequences follow, both deliberate:

* **No flux is recomputed.** In particular `lwupb` was evaluated by the radiation
  scheme from the previous skin temperature, so the slab's emission lags its own
  temperature by one radiation call. Recomputing it from the new `sst` would make
  the physics more current, but the atmosphere would then have absorbed one value
  while the slab was debited another, opening a gap in the budget. Currency is
  the wrong thing to optimise here.
* **The slab must be the only writer of `sst`.** Anything else that writes the
  sea surface temperature silently breaks the budget, so those combinations are
  rejected at startup in `physics_namelist_check` rather than allowed to run:

  | Rejected with `config_slab_ocean = true` | Reason |
  | --- | --- |
  | `config_sst_update = true` | would overwrite the slab temperature from an input stream |
  | `config_sstdiurn_update = true` | the diurnal skin adjustment would add unbudgeted energy |
  | `config_slab_depth <= 0` | not a physical heat capacity |

### Verifying it

`slab_energy_input` accumulates the net downward energy the slab has taken up.
Between any two output times, at every ocean cell:

```
slab_energy_input(t2) - slab_energy_input(t1)
    == rho_w * c_w * config_slab_depth * (sst(t2) - sst(t1))
```

must hold to round-off. Both fields are written to the `output` stream. Note that
`output` is a mutable stream, so this holds for a freshly generated
`streams.atmosphere`; if you are reusing an existing one, add them by hand.


## Diagnostics

`slab_energy_input` (J m^-2) is the only new field. The six terms that drive it
are all written to the `output` stream, so a drift can be attributed rather than
merely observed:

| Term | Sign |
| --- | --- |
| `swdnb` | + |
| `swupb` | − |
| `lwdnb` | + |
| `lwupb` | − |
| `hfx` | − |
| `lh` | − |

`hfx` and `lh` were already in the `output` stream; the four radiative terms were
added for the slab. They are written by every run, not only slab runs. That is
deliberate rather than an oversight: the `output` stream is declared with
`runtime_format="separate_file"`, so its variable list is generated as bare names
in `stream_list.atmosphere.output`, a format with nowhere to record a package.
A `packages` attribute on a `<var>` in that stream is dropped by `streams_gen`
with a build-time warning and the field is listed anyway, so gating them there
would have been an illusion. Four extra 2-d fields is small against the 3-d
fields the stream already carries.

`slab_energy_input` is different: it belongs to the `slab_ocean` var_struct, so
the *field* is gated. Without the slab it is never activated and the stream
manager skips it, whatever the stream list says.

There is no domain-integrated energy budget, and the net flux itself is formed
inside the update and not retained. For a global budget, integrate
`slab_energy_input` against `areaCell` in post-processing.


## Thermal response

With the default 2.5 m depth the heat capacity is

```
C = rho_w * c_w * h = 1.05e7 J m^-2 K^-1
```

so a sustained net flux of 100 W m^-2 warms the slab by about 0.83 K per day.
Against the linearised longwave feedback near 288 K (`4 sigma T^3` is about
5.4 W m^-2 K^-1) the surface relaxes on a timescale of roughly 22 days.

The shallow default is deliberate. A realistic 50 m mixed layer has twenty times
the inertia and would take decades of simulation to equilibrate; at 2.5 m the
surface settles in months, which is what makes the configuration usable for
idealized experiments. It is correspondingly less realistic: the seasonal and
synoptic response of the surface is far faster than the real ocean's.


## Behaviour worth knowing

* **First time step.** On a cold start the flux fields are zero (Registry real
  fields default to `0.0`) and the physics drivers have not yet run, so the first
  update is a no-op. No spurious energy enters the budget.
* **Restart.** `sst` and `slab_energy_input` are both in the restart stream, so a
  restarted run resumes with the correct surface state and a continuous energy
  accumulation.
* **Albedo.** `config_slab_albedo` is reapplied every time step, so it overrides
  whatever the albedo climatology or land-surface scheme would otherwise set over
  water. This is what keeps `swupb` consistent with `swdnb` in the budget.


## Limitations

* **Sea ice is not handled.** The slab updates every point with `landmask == 0`
  regardless of `xice`, and nothing stops the temperature falling below the
  freezing point of sea water. This is intended, but it means results are not
  physically meaningful anywhere ice would form in reality.
* **No Q-flux.** With no ocean heat transport the meridional temperature gradient
  is set purely by the local energy balance, so expect a warmer tropics and colder
  poles than observed.
* **No horizontal coupling.** Columns do not exchange heat, so there is nothing to
  damp small-scale SST structure other than the atmosphere itself.
* **Not run.** The slab has not been exercised in a simulation.
