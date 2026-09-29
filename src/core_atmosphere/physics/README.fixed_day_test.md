# Testing a fixed-day run

Step-by-step instructions for a short test run with `config_perpetual_julday`,
which holds the insolation and the ozone climatology at one day of the year. The
aim is to show that the option actually does what `README.steady_forcing.md`
says it does. None of this has been confirmed in a simulation yet (see Status in
that file), so treat this as the first run, not a routine one.

The test compares two runs from the same initial conditions: a control run with
the seasonal cycle and a locked run held at day 80, the vernal equinox. Five
days is enough. Top-of-atmosphere insolation depends only on astronomy, not on
the atmosphere, so a short run gives exact answers to check against.


## Step 1: build and initial conditions

Build the atmosphere core from this branch as usual:

```
make gnu CORE=atmosphere
```

Any real-data `init.nc` will work, for example the 240-km tutorial case. **Pick
a start date well away from day 80.** If the run starts near the equinox, the
control and the locked run have almost the same insolation and the test cannot
tell them apart. A start in October or in January (day 1–20 or 280–330) gives a
declination between 10 and 23 degrees, which is easy to see.

Start at 00 UTC so that output times fall on day boundaries. That keeps the
daily means in step 4 simple.


## Step 2: two run directories

Make two copies of a working run directory, `control/` and `locked/`, with
identical `init.nc`, `namelist.atmosphere` and `streams.atmosphere`. Set a
five-day run in both:

```
&nhyd_model
    config_start_time   = '2010-10-23_00:00:00'
    config_run_duration = '5_00:00:00'
/
```

In `locked/` only, add to `&physics`:

```
&physics
    config_perpetual_julday = 80.
/
```

The option has `in_defaults="false"`, so it is not in the generated namelist and
has to be typed in by hand. Leave `config_o3climatology` at its default. The
ozone check in step 4 depends on it being on.

Leave out `config_perpetual_all_physics` for this test. It changes only the
non-orographic gravity-wave drag and Noah-MP. Neither shows up in a five-day
run, and the gravity-wave source is off by default (`config_ngw_scheme`).


## Step 3: output

The fields needed for the checks are not written by default. In both run
directories, add these lines to `stream_list.atmosphere.diagnostics`:

```
swdnt
acswdnt
o3clim
```

The `diagnostics` stream writes every 3 hours by default. That is fine: step 4
uses the files at 00 UTC for daily means and the rest for the diurnal cycle.

Also make sure `config_bucket_update` is `'none'` (the default), so that
`acswdnt` is not reset partway through the run.

Run both.


## Step 4: checks

Daily-mean insolation comes from `acswdnt`, which is `swdnt` multiplied by the
time step and summed (J m^-2). The difference between two 00 UTC files divided
by 86400 s is the daily-mean downward TOA shortwave in W m^-2. Latitude and
other mesh fields are not in the diagnostics files, so read them from `init.nc`.

```python
import glob, netCDF4, numpy as np

lat = np.degrees(netCDF4.Dataset("init.nc").variables["latCell"][:])
bands = np.arange(-90, 91, 10)
idx   = np.digitize(lat, bands) - 1
S     = 1380.0                                  # solcon at day 80

def daily_means(run):
    files = sorted(glob.glob(f"{run}/diag.*_00.00.00.nc"))
    acc   = np.array([netCDF4.Dataset(f).variables["acswdnt"][0] for f in files])
    q     = np.diff(acc, axis=0) / 86400.0      # (day, nCells), W m^-2
    zm    = np.array([[q[d, idx == b].mean() for b in range(len(bands) - 1)]
                      for d in range(q.shape[0])])
    return q, zm

qc, zc = daily_means("control")
ql, zl = daily_means("locked")

np.set_printoptions(precision=1, suppress=True, linewidth=150)
print("control, zonal mean by day:\n", zc)
print("locked,  zonal mean by day:\n", zl)
```

**Check 1: the locked run does not change from day to day.** Each row of `zl`
should match the others. The largest difference between days,
`np.abs(ql - ql[0]).max()`, should be at round-off level. In the control, the
rows should drift, by several W m^-2 per day at mid-latitudes.

**Check 2: the locked run is symmetric about the equator.** At day 80 the
declination is zero, so each band in `zl` should equal its mirror across the
equator, `zl[:, ::-1]`. The control run should not be symmetric, because an
October start favours the Southern Hemisphere.

**Check 3: the locked run matches the analytic value.** With zero declination,
the daily-mean TOA insolation is `S cos(lat) / pi`:

```python
expected = S * np.cos(np.radians(lat)) / np.pi
print("max relative error:",
      np.abs(ql[0] - expected).max() / expected.max())
```

Expect agreement to within about 1%, and about 439 W m^-2 at the equator. The
error is not zero because `swdnt` is computed once per radiation call
(`config_radtsw_interval`, 30 minutes by default) and held between calls. An
error of several percent, or a clear hemispheric asymmetry, means the
declination is not being locked. If the error is about 0.7% everywhere, the
solar constant is `solcon_0 = 1370` instead of the day-80 value of 1380, which
points to a build of the old `config_perpetual_equinox` code.

**Check 4: the diurnal cycle is still there.** Take the 3-hourly `swdnt` from
any single file. Roughly half the cells should be exactly zero (night), and the
largest value should be close to `S`. If `swdnt` is the same at every longitude,
something has flattened the hour angle, which this option is not meant to do.

**Check 5: ozone is frozen.** Compare `o3clim` between the first and last
diagnostics files:

```python
def o3(run, which):
    f = sorted(glob.glob(f"{run}/diag.*.nc"))[which]
    return netCDF4.Dataset(f).variables["o3clim"][0]

print("locked  change:", np.abs(o3("locked", -1)  - o3("locked", 1)).max())
print("control change:", np.abs(o3("control", -1) - o3("control", 1)).max())
```

Use index 1, not 0, because the first file is written before the first radiation
call. In the locked run the change should be exactly zero: the climatology is
interpolated to the same day every time. In the control it should be small but
not zero, since the monthly climatology is interpolated to the current day. If
the locked run's ozone changes, the RRTMG path in `physics_timetracker` is not
receiving `radt_julday`.

This check covers RRTMG only, which every bundled physics suite uses. The CAM
radiation path (`config_radt_sw_scheme = 'cam_sw'`) freezes ozone by a different
route and needs its own run if you use it.


## Step 5: startup checks

Two bad configurations should stop the model at startup, in
`physics_namelist_check`. Each takes seconds to test:

| Setting | Expected message |
| --- | --- |
| `config_perpetual_julday = 400.` | `config_perpetual_julday must be a day of the year, less than 367.` |
| `config_perpetual_all_physics = true`, with no `config_perpetual_julday` | `config_perpetual_all_physics requires a positive config_perpetual_julday` |

Look for the message in `log.atmosphere.0000.err`.


## Other days

To test a different day, change `config_perpetual_julday` and the value of `S`.
The table in `README.steady_forcing.md` gives the declination and solar constant
for days 80, 172, 266 and 355. Checks 2 and 3 as written only apply at an
equinox. At another day, the daily mean still has to be constant from day to
day (check 1). It also has to match the control's daily mean from a run started
on that day of the year, which is a stronger test than the analytic formula.


## If a check fails

| Symptom | Likely cause |
| --- | --- |
| Locked and control runs are identical | the option is not in `&physics`, is in a different record, or is non-positive |
| Locked insolation drifts from day to day | the build does not include this branch; rebuild after `make clean CORE=atmosphere` |
| Locked insolation about 0.7% low everywhere | the build has the old `config_perpetual_equinox` code, which used `solcon_0` |
| `acswdnt` jumps backwards between files | bucket reset is on; set `config_bucket_update = 'none'` |
| `o3clim` changes in the locked run | ozone is not using the locked day, or `config_radt_*_scheme` is CAM |
