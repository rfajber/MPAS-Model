# Initial conditions with SST and sea ice from a separate file

This document describes how to build an MPAS initial-conditions file (`init.nc`)
where all the meteorology comes from one data set and the sea surface
temperature (`sst`) and sea-ice fraction (`xice`) come from another. For example,
the atmosphere and land from ERA5 and the ocean surface from OISST, or from an
idealized or perpetual-day SST field.

No code changes are needed. `init_atmosphere` already supports this through
`config_input_sst`.


## How it works

The real-data initialization (`config_init_case = 7`) reads two sets of WPS
intermediate files at the start time:

| Namelist option | File read | What is taken from it |
| --- | --- | --- |
| `config_met_prefix` | `<met_prefix>:YYYY-MM-DD_HH` | everything: 3-D atmosphere, surface pressure, soil, snow, `SKINTEMP`, and a first `SEAICE` |
| `config_sfc_prefix` | `<sfc_prefix>:YYYY-MM-DD_HH` | `SST` (or `SKINTEMP`) and `SEAICE` only, when `config_input_sst = true` |

The order of operations in `init_atmosphere` is:

1. **Meteorology** (`init_atm_case_gfs` in `mpas_init_atm_cases.F`).
   - Everything is interpolated from the met file.
   - `sst` is set to the met file's `SKINTEMP` everywhere. Case 7 never reads a
     field called `SST` from the met file.
   - `xice` is set from the met file's `SEAICE`.
2. **SST and sea ice from the second file**, when `config_input_sst = true`
   (`physics_initialize_real` in `mpas_atmphys_initialize_real.F`).
   - `interp_sfc_to_MPAS` (`mpas_init_atm_surface.F`) reads
     `<sfc_prefix>:YYYY-MM-DD_HH` for the start time and **replaces** `sst`
     (from a field named `SST` or `SKINTEMP`) and `xice` (from `SEAICE`). Each
     is replaced only if the file contains it.
   - Water cells with `sst < 271 K` are then set to `xice = 1`, and `xice` is
     clipped to [0, 1].
3. **Consistency with the new ocean surface** (`physics_init_sst`).
   - At open-water cells (`landmask = 0`, `xice` below the sea-ice threshold),
     `skintemp` is set to the new `sst`.
   - Any `xice` on land cells is reset to 0.
4. **Sea-ice cells** (`physics_init_seaice`). Cells whose new `xice` is at or
   above the threshold become sea-ice points:
   - `ivgtyp`, `isltyp`, `snoalb`, `vegfra` and `xland` are set to ice values;
   - the ice temperature and moisture profiles (`tslb`, `smois`, `sh2o`) are set;
   - `tmn` is set to 271.4 K.

   Water cells below the threshold get `xice = 0` and no snow.

So the land surface, soil, snow and atmosphere all come from the met file, while
`sst`, `xice`, and everything derived from them are consistent with the second
file. The land/water mask always comes from the MPAS static file (`landmask`),
not from either data set.


## Step 1: make the two sets of intermediate files

Both files are WPS intermediate format and must have the same date and hour as
`config_start_time`.

### Option A: both data sets are GRIB (WPS `ungrib`)

Run `ungrib.exe` twice, with a different Vtable and `prefix` each time:

```
# first pass: atmosphere and land
&ungrib
    out_format = 'WPS',
    prefix     = 'ERA5',
/
# link the atmospheric GRIB files and Vtable.ECMWF (or the Vtable for your source), run ungrib.exe

# second pass: ocean surface
&ungrib
    out_format = 'WPS',
    prefix     = 'SST',
/
# link the SST/ice GRIB files and Vtable.SST (or equivalent), run ungrib.exe
```

This gives `ERA5:2020-01-01_00` and `SST:2020-01-01_00`.

**Field names.** The second file must contain a field named exactly `SST` or
`SKINTEMP`, in kelvin, and `SEAICE`, as a fraction from 0 to 1. Standard WPS
Vtables use these names. Check with `rd_intermediate.exe SST:2020-01-01_00`.

The met file must still contain `SKINTEMP`. It supplies the skin temperature over
land and sea ice, and the soil initialization uses it.

### Option B: the SST and ice are in NetCDF (or you made them yourself)

Use `scripts/write_sst_intermediate.py`. It reads SST, and optionally sea ice
and a land mask, from a NetCDF file on a regular latitude-longitude grid, and
writes the intermediate file:

```
./scripts/write_sst_intermediate.py oisst.nc 2020-01-01_00:00:00 \
    --sst-var sst --ice-var ice                 # writes SST:2020-01-01_00
```

Options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--sst-var` | required | SST variable, in K or °C |
| `--ice-var` | none | sea-ice variable, as a fraction or %. Without it no `SEAICE` is written, and `xice` stays from the met file |
| `--landsea-var` | none | land mask (1 = land), written as `LANDSEA` |
| `--mask-ice-with-sst` | off | treat points where the SST is missing as land, and mask the sea ice with them (see below). Needs `--ice-var`; cannot be combined with `--landsea-var` |
| `--lat`, `--lon` | `lat`, `lon` | names of the 1-D coordinate variables |
| `--time-index` | 0 | which time to take if the variables have a time axis |
| `--prefix` | `SST` | output prefix; must match `config_sfc_prefix` |

The script needs `numpy` and `netCDF4`. It prepares the data so that
`init_atmosphere` interprets it correctly:
- **Units.** SST in °C is converted to kelvin, judged from the `units`
  attribute, or from the values if there is none. Sea ice in % is converted to a
  fraction and clipped to [0, 1]. **Any water cell that ends up below 271 K is
  made into sea ice**, so an SST left in °C would cover the ocean in ice.
- **Missing values.** Masked, NaN and fill values become `-1.e30`, the missing
  value, so the interpolation searches for the nearest valid point instead of
  using them. Never leave land or missing points as 0 or a °C fill value.
- **Grid orientation.** Latitude is flipped to run south to north, and longitude
  is shifted to [0, 360) and sorted. The grid must be regular and **global in
  longitude**, because the interpolation wraps around the globe; the script stops
  otherwise.
- **Check.** It reads the file back and prints each field's grid and value range.

About the land mask: `LANDSEA` is used only to exclude source land points when
interpolating `SEAICE`. It is not used for `SST`, which relies on the `-1.e30`
missing values instead.

**Masking the sea ice with the SST's land mask** (`--mask-ice-with-sst`). SST
and sea-ice data sets often have different land masks, or the sea-ice field has
none and holds 0 over land. Near coasts, `init_atmosphere` would then blend
land zeros into the ice fraction. With `--mask-ice-with-sst`, every point where
the SST is missing is taken as land, and the sea ice is masked with it in both
ways `init_atmosphere` recognizes:
- the ice is set to `-1.e30` there;
- the mask is written as `LANDSEA`.

Sea ice is then interpolated only from points with a valid SST, and the two
fields share one coastline.

Do not use it if your SST product leaves the SST missing under sea ice: those
points would be treated as land and their ice dropped. The script warns and
counts them when it sees ice where the SST is missing.

The script writes version 5 of the format, as read by
`mpas_init_atm_read_met.F`: big-endian Fortran sequential records, five per
field, on a lat-lon grid (`iproj = 0`). The reader also accepts Gaussian and
polar-stereographic grids, but the script does not write them.


## Step 2: namelist and streams

In `namelist.init_atmosphere`:

```
&nhyd_model
    config_init_case    = 7
    config_start_time   = '2020-01-01_00:00:00'
    config_stop_time    = '2020-01-01_00:00:00'
/
&data_sources
    config_met_prefix   = 'ERA5'      ! atmosphere, land, snow, SKINTEMP
    config_sfc_prefix   = 'SST'       ! SST and SEAICE
/
&preproc_stages
    config_static_interp = false      ! assuming static.nc already exists
    config_native_gwd_static = false
    config_vertical_grid = true
    config_met_interp    = true
    config_input_sst     = true       ! <-- read SST and SEAICE from config_sfc_prefix
    config_frac_seaice   = true       ! must match the model's config_frac_seaice
/
```

In `streams.init_atmosphere`, the `input` stream is your `static.nc`, and the
`output` stream is the new `init.nc`, as for any case-7 run.

Put both intermediate files in the run directory and run `init_atmosphere_model`.

Keep the rest of `&nhyd_model`, `&dimensions`, `&vertical_grid` and
`&interpolation_control` as you would normally for your met data. For example,
`config_nfglevels` and `config_nfgsoillevels` must match the met file, not the
SST file.


## Step 3: check the log and the output

In `log.init_atmosphere.0000.out` you should see:

```
--- read sea-surface temperature from auxillary file:
Processing file SST:2020-01-01_00
```

If instead you see `Error opening surface file SST:2020-01-01_00`, the file name
or its date/hour does not match `config_start_time`.

Then compare `init.nc` with the source data:

```python
import netCDF4, numpy as np
d = netCDF4.Dataset('init.nc')
sst, xice, tsk = d['sst'][0], d['xice'][0], d['skintemp'][0]
land = d['landmask'][:] == 1
water = ~land & (xice < 0.02)
print('SST over open water: ', sst[water].min(), sst[water].max())
print('skintemp == sst there:', np.allclose(tsk[water], sst[water]))
print('ice cells:', (xice >= 0.02).sum(), ' max xice on land:', xice[land].max())
```

With `config_frac_seaice = false` the ice threshold is 0.5, not 0.02.

`skintemp` should equal `sst` at open-water cells, `xice` should be 0 on land,
and the SST range should look like your SST data set, not like the met file's
skin temperature.


## Using the file in the model

- **Fixed SST.** With `config_sst_update = false` (the default), `sst` and `xice`
  stay at the values in `init.nc` for the whole run.
- **Time-varying SST.** To make them vary in time from the same SST source, also
  build a surface-update file with `config_init_case = 8`, using the same
  `config_sfc_prefix`, `config_fg_interval` and the `surface` output stream. Then
  set `config_sst_update = true` in the model, point its `surface` input stream
  at that file, and set that stream's `input_interval` (it is `none` by default).
- **Slab ocean** (`README.slab_ocean.md`). The slab starts from the `sst` in
  `init.nc`, and `config_sst_update` must stay `false`.
- **Fixed-day runs** (`README.fixed_day_realistic.md`). A separately prepared SST
  (for example a climatology for the chosen day) is exactly what this is for.
  Keep `config_sst_update = false`.


## Known gotchas

| Symptom | Likely cause |
| --- | --- |
| Log says `Error opening surface file ...` | file missing, wrong prefix, or its date/hour differs from `config_start_time` |
| No "read sea-surface temperature from auxillary file" line | `config_input_sst` not set, or not in `&preproc_stages` |
| `sst` in `init.nc` looks like the met file's skin temperature | the SST file has no field named `SST` or `SKINTEMP` (check with `rd_intermediate.exe`) |
| Sea ice unchanged from the met data | the SST file has no `SEAICE` field. `xice` then keeps the met-file value |
| Ocean covered in sea ice | SST in °C, or land/missing points set to 0 instead of `-1.e30` (any water cell below 271 K becomes ice) |
| Stripes or a seam at the date line | grid not global in longitude, or latitude stored north to south |
| Crash or nonsense reading the file | little-endian or wrong record layout. The file must be big-endian, version 5 |
| Unexpected sea ice even with `SEAICE` in the SST file | a leftover `SEAICE_FRACTIONAL:YYYY-MM-DD_HH` file in the run directory (see below) |

**`SEAICE_FRACTIONAL` files.** Case 7 also looks for a file named
`SEAICE_FRACTIONAL:YYYY-MM-DD_HH` in the run directory. If it exists, its
`SEAICE` replaces the met file's before step 2. With `config_input_sst = true`,
a `SEAICE` field in the SST file overrides it again. It matters only if your SST
file has no `SEAICE`, and it only supports polar-stereographic grids. Remove any
stray copies to avoid confusion.


## Status

This has not been run. It is based on reading the code in
`mpas_init_atm_cases.F` (`init_atm_case_gfs`), `mpas_atmphys_initialize_real.F`
(`physics_initialize_real`, `physics_init_sst`, `physics_init_seaice`),
`mpas_init_atm_surface.F` (`interp_sfc_to_MPAS`) and `mpas_init_atm_read_met.F`.

`scripts/write_sst_intermediate.py` was tested on synthetic NetCDF files:
- an OISST-like file (°C, %, `_FillValue` land, a time axis);
- an ERA5-like file (K, NaN land, latitude north to south, longitude −180 to 180).
- a file with SST masked on land and under ice, used for `--mask-ice-with-sst`.
  There the ice became missing exactly where the SST is, `LANDSEA` matched,
  and the ice-without-SST warning fired.

An independent decoder of its output reproduced the input values exactly after
the unit, orientation and fill-value changes. The output has not yet been read
by `init_atmosphere` itself or by WPS's `rd_intermediate.exe`, so check the log
messages in step 3 on first use.
