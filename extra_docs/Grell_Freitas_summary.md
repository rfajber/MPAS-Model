# Grell–Freitas convection in MPAS-Atmosphere

## 1. Overview

MPAS runs a port of the WRF v3.5-era Grell–Freitas (GF) scheme, `module_cu_gf.mpas.F` ([module_cu_gf.mpas.F:1](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L1)). The code is a mass-flux scheme descended from Grell (1993) and Grell & Dévényi (2002, G3). Grell & Freitas (2014, ACP 14) added a "smooth transition to cloud-resolving scales" following Arakawa et al. (2011). The module header describes the scale-aware part as on and tested down to about 3 km, and the aerosol coupling as off. The two schemes, deep and shallow, run in the same call:

| Component | Routine | Precipitates? | Downdrafts? | Scale-aware? | Closure selector |
|---|---|---|---|---|---|
| deep | `cup_gf` | yes (c0 = 2e-3 m⁻¹) | yes | yes, via σ | `config_gfconv_closure_deep` |
| shallow | `cup_gf_sh` | no (c0 = 0) | no | no | `config_gfconv_closure_shallow` |

MPAS changes from WRF, by Laura Fowler, 2014–2016 ([module_cu_gf.mpas.F:27](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L27)):
- `dx` is replaced by per-cell `dxCell` and `areaCell`.
- `ichoice` is split into separate deep and shallow closure choices.
- The initial updraft radius is changed from 0.1/ε to 0.2/ε.
- The driver `gfdrv` is renamed `cu_grell_freitas`.
- The moisture-convergence calculation is fixed to use ρ interpolated to w levels.

**Where it is selected.** Setting `config_convection_scheme = 'cu_grell_freitas'` turns on the `cu_grell_freitas_in` package ([mpas_atmphys_packages.F:130](../src/core_atmosphere/physics/mpas_atmphys_packages.F#L130)). GF is also what the `convection_permitting` suite picks when the option is left at `suite` ([mpas_atmphys_control.F:169](../src/core_atmosphere/physics/mpas_atmphys_control.F#L169)). That suite pairs it with Thompson, MYNN PBL/surface layer, Noah and RRTMG. The `mesoscale_reference` suite uses `cu_ntiedtke` instead.

## 2. Deep convection (`cup_gf`)

### 2.1 Environment and "forced" profiles

Two sets of profiles go into the scheme ([module_cu_gf.mpas.F:285](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L285)):
- **Current state** (`t2d`, `q2d`): T and qv at the start of the physics step.
- **Forced state** (`tn`, `qo`): the current state advanced one `dt` with the non-convective tendencies:
  - tn = T + (θ̇_dyn + θ̇_rad + θ̇_PBL)·Π·dt
  - qo = qv + (q̇v_dyn + q̇v_PBL)·dt

The CAPE-type closures compare these two states. Cloud work functions are computed with moist static energy on "cup" (interface) levels (`cup_env`, `cup_env_clev`). Below `tcrit = 258 K` saturation is taken over ice.

### 2.2 Trigger and cloud model

1. **Source level k22** is the level of maximum h in the forced profile, below `zkbmax = 4000 m` AGL.
2. **Cloud base kbcon** is the first level where h(k22) exceeds h*. If the parcel has to be lifted more than `cap_max = 75 hPa` to reach it, k22 is moved up one level and the search repeats ([module_cu_gf.mpas.F:838](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L838), `cup_kbcon`). The test is repeated at 50 and 25 hPa caps (`ierr2`, `ierr3`). Over water, failing the 50 hPa test zeroes every closure except the vertical-mass-flux members (see 2.4).
3. **Entrainment.**
   - Base rate: ε₀ = 7e-5 m⁻¹ ([module_cu_gf.mpas.F:753](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L753)).
   - Humidity dependence: ε(k) = ε₀(1.3 − RH), so drier environments dilute more.
   - Mass-flux profile: detrainment is solved so the normalized updraft mass flux rises parabolically from 0.1 at the surface to 1 at kbcon, then decays back to 0.2 near cloud top. Profiles are 1-2-1 smoothed, and zu is floored at 0.05.
4. **Cloud top ktop** is the first level where the entraining updraft loses buoyancy. Clouds shallower than `depth_min = 1000 m` are rejected, as are clouds with ktop < kbcon+2.
5. **Updraft microphysics** (`cup_up_moisture`, `autoconv = 1`):
   - Condensate above saturation converts to rain at `c0 = 2e-3 m⁻¹` ([module_cu_gf.mpas.F:3876](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L3876)).
   - The Berry option (`autoconv = 2`) and aerosol-dependent evaporation (`aeroevap > 1`) are hard-wired off. The CCN value, fixed at 1500 ([module_cu_gf.mpas.F:264](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L264)), therefore has no effect.
   - Condensate detrains at levels between kbcon and ktop, and at ktop itself.
6. **Downdrafts.**
   - Downdrafts start at the h* minimum below min(0.6·cloud depth, `zcutdown = 3000 m`) AGL. Their mass detrains over the lowest `z_detr = 1250 m`.
   - The downdraft/updraft mass-flux ratio ε = 1 − ½(PE_shear + PE_cloudbase), clipped to 0.1–0.9 ([module_cu_gf.mpas.F:2090](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L2090)). This is the classic Fritsch–Chappell / Zhang–Fritsch precipitation efficiency as functions of cloud-layer wind shear and cloud-base height.
   - Evaporation of rain in the downdraft cools and moistens the subcloud layer.
7. **Environmental response.** Compensating subsidence and detrainment give per-unit-mass-flux tendencies `dellat`, `dellaq`, `dellaqc`, `dsubt` and `dsubq`. The "cloud-modified" CAPE xaa0 comes from applying those tendencies over a small `mbdt = 0.003·dt`. The sensitivity xk = (xaa0 − aa1)/mbdt feeds the CAPE closures.

### 2.3 Scale awareness (σ)

The updraft radius is taken as R = 0.2/ε₀ ≈ 2.86 km, and the updraft area fraction is f = πR²/Δx². When f exceeds 0.7 it is capped at 0.7, and R and ε are recomputed from the cap, so ε *grows* at very fine Δx ([module_cu_gf.mpas.F:753](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L753)). The closure mass flux is multiplied by σ = (1 − f)² ([module_cu_gf.mpas.F:3721](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L3721)).

| Δx | f | σ | ε (m⁻¹) |
|---|---|---|---|
| 60 km | 0.007 | 0.99 | 7.0e-5 |
| 30 km | 0.028 | 0.94 | 7.0e-5 |
| 15 km | 0.114 | 0.79 | 7.0e-5 |
| 10 km | 0.256 | 0.55 | 7.0e-5 |
| 8 km | 0.401 | 0.36 | 7.0e-5 |
| ≤ 6.05 km | 0.7 (cap) | 0.09 | 7.1e-5 (6 km), 1.4e-4 (3 km) |

**Δx in MPAS.** The driver sets Δx = `config_len_disp` / meshDensity^¼ ([mpas_atmphys_driver_convection.F:718](../src/core_atmosphere/physics/mpas_atmphys_driver_convection.F#L718)), not the local cell spacing. `config_len_disp` defaults to `nominalMinDc` from the grid file, so on variable-resolution meshes σ follows the mesh-density function. `areaCell` is passed through to `cup_gf` but not used.

**Resolved-ascent cutoff.** When σ < 0.091 (i.e. only at the f cap, Δx ≲ 6 km), the scheme checks for grid-scale ascent at saturation, ω < 0 and RH ≥ 0.99, at k22 or anywhere between k22 and kbcon. If it finds it, it flags `ierr = 1200` ([module_cu_gf.mpas.F:948](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L948)). The intent is to hand saturated resolved updrafts to the microphysics. See §7 for how the flag is later reset.

### 2.4 Closures (`cup_forcing_ens_3d`)

There are 16 closure members, but the scheme runs with `maxens = maxens2 = 1` and `irandom = 0`. In that configuration several members come out identical, so there are effectively four families ([module_cu_gf.mpas.F:2505](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L2505)).

| Members | Family | Mass flux m_b |
|---|---|---|
| 1, 2, 3, 13 | Grell (1993) quasi-equilibrium | −[(aa1 − aa0)/dt] / xk: CAPE generated by large-scale forcing is consumed within one step |
| 4, 5, 6 | low-level vertical mass flux (Brown 1979) | max(−ω/g) between level 2 and kbcon (members 4 and 5 are overwritten with 6) |
| 14 | low-level vertical mass flux | −ω/g at k22 |
| 7, 8, 9, 15 | moisture convergence (Krishnamurti et al. 1983) | mconv / (precip per unit m_b); mconv = ∫ ω ∂q/∂z /g, floored at 0 ([module_cu_gf.mpas.F:325](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L325)) |
| 10, 11, 12, 16 | CAPE removal (Arakawa–Schubert / Kain–Fritsch style) | −[aa0/1200 s] / xk, used only if xk < 0 |

- `ichoice_deep = 0` (the default) averages all 16 members. `ichoice_deep = n` (1–16) copies member n into every slot, so only that closure is used.
- The result is multiplied by σ. Mass fluxes above 100 kg m⁻² s⁻¹ are rejected (`ierr = 19`).
- **Land/ocean asymmetry** (`xland` > 1.5 counts as water). If the large-scale forcing *decreases* CAPE (aa1 < aa0), members 1–3, 10–13 are zeroed everywhere, and *all* members are zeroed over water ([module_cu_gf.mpas.F:2765](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L2765)).
- **Cap retests.** Over water, all members except the vertical-mass-flux ones (4–6, 14) are also scaled by `ens_adj`, which is 0 whenever the 50 hPa cap retest (`ierr2`) failed.
- Members with zero precipitation are discarded before averaging.

### 2.5 Limiters (`neg_check`)

The deep tendencies, including precipitation, are rescaled per column in two ways ([module_cu_gf.mpas.F:3496](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L3496)):
- so that heating stays below 600 K day⁻¹ and cooling stays below 300 K day⁻¹;
- so that no level ends up with qv < 1e-10 after one `dt`.

Precipitation is removed from the column immediately. GF has no rain or snow tendencies, unlike KF.

## 3. Shallow convection (`cup_gf_sh`)

It is always on: `ishallow = 1` is a compile-time parameter ([mpas_atmphys_vars.F:340](../src/core_atmosphere/physics/mpas_atmphys_vars.F#L340)).

- **Forcing.** The shallow scheme sees only the PBL tendencies:
  - tshall = T + θ̇_PBL·Π·dt
  - qshall = qv + q̇v_PBL·dt
  - dh/dt_PBL = c_p θ̇_PBL Π + L q̇v_PBL
- **Trigger.** k22 is the h maximum below 4 km; if that fails and `kpbl > 5`, k22 falls back to kpbl. The cap is 25 hPa, or down to p(kpbl) when kpbl > 5 ([module_cu_gf.mpas.F:4465](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L4465)). Minimum depth is 50 m.
- **Cloud model.** Entrainment is 1e-3 m⁻¹ × (1.3 − RH) ([module_cu_gf.mpas.F:4372](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L4372)). There is no rain (c0 = 0) and no downdraft. Detrained condensate is added to the *vapor* tendency ([module_cu_gf.mpas.F:4840](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L4840)), so shallow convection outputs only θ and qv tendencies; its qc/qi tendencies are zero.
- **Closures** (9 members, `xmbmax = 0.1 kg m⁻² s⁻¹`) ([module_cu_gf.mpas.F:4943](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L4943)):

| Members | Family | m_b |
|---|---|---|
| 1–3 | quasi-equilibrium on PBL forcing | −(aa1 − aa0)/(xk·dt), only if aa1 > aa0 and aa1 > 0 |
| 4–6 | CAPE removal over an advective time scale | −aa0/(xk·τ), τ = Δx/25 m s⁻¹ ([module_cu_gf.mpas.F:249](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L249)) |
| 7–9 | boundary-layer quasi-equilibrium (BLQE, Raymond 1995 via S. Freitas) | ∫ (dh/dt)_PBL dp/g below cloud base ÷ (h_c − h̄) at kbcon, only if k22 ≤ kpbl |

  `ichoice_shallow > 0` picks that member, `= 0` averages all 9, and `< 0` takes the maximum. The default is **8, which is BLQE**.
- **Limiters.** Heating/cooling is capped at ±100 K day⁻¹ by scaling m_b, and m_b is reduced so that qv stays above 1e-12.
- **No σ scaling.** `sig` is declared but never applied, so shallow convection stays at full strength at every Δx.

## 4. Coupling to MPAS

### 4.1 Call sequence and timing

Per model step, `physics_driver` ([mpas_atmphys_driver.F:244](../src/core_atmosphere/physics/mpas_atmphys_driver.F#L244)) calls radiation (when its alarm rings), then the surface layer, LSM/sea ice, PBL, GWDO and **convection** ([mpas_atmphys_driver.F:355](../src/core_atmosphere/physics/mpas_atmphys_driver.F#L355)). After that the dynamics runs, with the stored physics tendencies applied on every RK stage, and microphysics runs at the end of `atm_srk3`.

- **Call frequency.** GF runs when `l_conv` is true: every step if `config_conv_interval = 'none'`, otherwise on an alarm ([mpas_atmphys_manager.F:359](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L359)). Between calls, `rthcuten` etc. stay in `tend_physics` and are reapplied unchanged.
- **Time step used inside GF.** GF is always passed `dt = dt_dyn`, not `dt_cu` ([mpas_atmphys_driver_convection.F:484](../src/core_atmosphere/physics/mpas_atmphys_driver_convection.F#L484)). That value sets the forcing projection (tn, qo), the QE closure rate and the negative-qv checks.
- **Precipitation accounting.** `raincv` = pratec·dt_dyn is the step precipitation. `cuprec` = pratec, and `rainc` accumulates cuprec·dt_dyn every step, with the `config_bucket_rainc` bucket ([mpas_atmphys_driver_convection.F:1090](../src/core_atmosphere/physics/mpas_atmphys_driver_convection.F#L1090)).

### 4.2 Inputs, and which parameterization supplies them

| GF argument | MPAS source | Supplied by |
|---|---|---|
| u, v, w, T, qv, ρ, p, Π, dz, zgrid(1) | state/diag via `mpas_atmphys_interface` | dynamics |
| `rthften` ← `rthdynten` | advective θ_m tendency from `atm_compute_dyn_tend` ([mpas_atm_time_integration.F:6976](../src/core_atmosphere/dynamics/mpas_atm_time_integration.F#L6976)), converted to dry θ in `atm_srk3` ([mpas_atm_time_integration.F:2807](../src/core_atmosphere/dynamics/mpas_atm_time_integration.F#L2807)) | previous step's dynamics |
| `rqvften` ← `rqvdynten` | (qv_new − qv_old)/dt over the whole previous step, before microphysics; **zero unless `config_monotonic = true`** | previous step's dynamics and physics |
| `rthraten` | `rthratenlw + rthratensw` | RRTMG / CAM radiation (last radiation call) |
| `rthblten`, `rqvblten` | PBL tendencies from this step | YSU / MYNN |
| `kpbl` | PBL top index | YSU / MYNN |
| `xland` | land/water mask | static / sea ice |
| `hfx`, `qfx`, `gsw` | surface fluxes, SW at the surface | surface layer / LSM / radiation. **Passed but unused**: `use_excess = 0` disables the flux-based parcel excess, and `gsw` is never read. |
| Δx | `config_len_disp / meshDensity^¼` | mesh |

### 4.3 Outputs, and where they go

| Output | Destination | Consumer |
|---|---|---|
| `rthcuten`, `rqvcuten` | `tend_physics` → `tend_th`/`tend_scalars` × mass ([mpas_atmphys_todynamics.F:429](../src/core_atmosphere/physics/mpas_atmphys_todynamics.F#L429)) | dynamics (every RK stage) |
| `rqccuten` (T ≥ 258 K) or `rqicuten` (T < 258 K) | same | dynamics → the microphysics then processes the detrained qc/qi |
| `cuprec`, `raincv`, `rainc` | `diag_physics` | output; LSM (via total rain) |
| `qc_cu`, `qi_cu` (in-cloud condensate, gdc/gdc2) | `diag_physics` | output only; no radiation or cloud-fraction scheme reads them |
| `cubot`, `cutop`, `xmb_total`, `xmb_shallow`, `k22_shallow`, `kbcon_shallow`, `ktop_shallow` | `diag_physics` | output only |
| `rucuten`, `rvcuten`, `ktop_deep` | allocated and copied in and out | never computed by GF; momentum transport is applied only for Tiedtke ([mpas_atmphys_todynamics.F:449](../src/core_atmosphere/physics/mpas_atmphys_todynamics.F#L449)) |

Deep and shallow tendencies are summed as rthcuten = [outts + (subt + outt)·cuten]/Π, where `cuten = 1` only if the deep scheme produced precipitation ([module_cu_gf.mpas.F:436](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L436)). The shallow tendencies are always added.

### 4.4 Interactions with the other parameterizations

- **PBL (YSU/MYNN).**
  - PBL runs just before GF in the same step. Its θ/qv tendencies feed both the deep forcing and all of the shallow forcing, including BLQE, and `kpbl` sets the shallow source level and cap.
  - With `config_pbl_scheme = 'off'`, the shallow QE and BLQE closures have no forcing, so the default shallow closure (8) never fires.
  - MYNN also has its own EDMF shallow mass flux, which runs alongside GF shallow; the code does not coordinate the two.
- **Radiation.**
  - Radiative θ tendencies from the most recent radiation call are part of the deep forcing, so a radiation call that destabilizes the column raises the QE/Grell closure.
  - Radiation sees convection only through the qc/qi that GF detrains into the grid-scale fields and through the T/qv changes. RRTMG's `cld_fraction` uses the grid-scale qc/qi/RH, not `qc_cu`/`qi_cu`.
- **Microphysics (Thompson, WSM6, …).**
  - GF adds detrained condensate to qc or qi only, split at 258 K. No hydrometeor numbers are updated.
  - Under Thompson, qi is added without ni, which raises the mean ice size and speeds ice → snow conversion. Under `mp_thompson_aerosols`, nc is not updated either.
  - Rain produced in the updraft falls out immediately and is not passed to the microphysics.
  - Microphysics runs after the dynamics has applied the convective tendencies. The `rqvdynten` for the next GF call is computed just *before* microphysics, so it excludes microphysical condensation.
- **Dynamics.**
  - The deep closures use grid-scale w (ω = −ρ g w, with ρ interpolated to w levels) and the advective tendencies from the previous step. Resolved ascent at saturation can switch GF off at Δx ≲ 6 km (§2.3).
  - `config_monotonic = false` silently zeroes `rqvdynten`, which removes the moisture-advection part of the forcing.
- **Surface / LSM / slab ocean.** These act only indirectly, through the PBL tendencies, `xland` (which switches on the ocean-only closure suppression) and `kpbl`. `hfx`/`qfx` are passed but not used.

## 5. Namelist options

Options in `&physics` that affect GF:

| Option | Default | Allowed values | Effect on GF |
|---|---|---|---|
| `config_convection_scheme` | `suite` | `suite`, `cu_grell_freitas`, `cu_kain_fritsch`, `cu_tiedtke`, `cu_ntiedtke`, `off` | `cu_grell_freitas` turns on GF and the `cu_grell_freitas_in` package ([Registry.xml:2526](../src/core_atmosphere/Registry.xml#L2526)) |
| `config_physics_suite` | `mesoscale_reference` | `mesoscale_reference`, `convection_permitting`, `none` | `convection_permitting` picks GF when `config_convection_scheme = 'suite'` |
| `config_gfconv_closure_deep` | `0` | Registry lists only `0`; code accepts 1–16 | 0 = average of all 16 closure members; n = use member n only (table in §2.4) ([Registry.xml:2591](../src/core_atmosphere/Registry.xml#L2591)) |
| `config_gfconv_closure_shallow` | `8` | Registry lists only `8`; code accepts 1–9, 0 and negative values | 1–3 QE, 4–6 CAPE removal over Δx/25, 7–9 BLQE; 0 = average; < 0 = maximum ([Registry.xml:2596](../src/core_atmosphere/Registry.xml#L2596)) |
| `config_conv_interval` | `none` | `DD_HH:MM:SS` or `none` | interval between GF calls; tendencies are held in between. The scheme itself still uses `dt_dyn` ([Registry.xml:2491](../src/core_atmosphere/Registry.xml#L2491)) |
| `config_bucket_rainc` | `100.0` | positive real (mm) | `rainc` bucket size (with `config_bucket_update`) |
| `config_len_disp` (`&nhyd_model`) | `0.0` → `nominalMinDc` | positive real (m) | sets the Δx for σ and for the shallow τ = Δx/25 |
| `config_monotonic` (`&nhyd_model`) | `true` | logical | when false, `rqvdynten` = 0 and the moisture-advection forcing is lost |

Hard-wired in the source, not in the namelist ([module_cu_gf.mpas.F:70](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L70)):
- **Scheme switches:** `ishallow = 1`, `autoconv = 1` (constant c0), `aeroevap = 1`, `use_excess = use_excess_sh = 0`, `training = 0`.
- **Physical parameters:** `tcrit = 258 K`, `ccn = 1500`, `ccnclean = 250`, `beta = 0.02`.
- **Deep cloud model:** ε₀ = 7e-5 m⁻¹, `cap_max` = 75 hPa, `depth_min` = 1000 m, `zkbmax` = 4000 m, `zcutdown` = 3000 m, `z_detr` = 1250 m, `edtmin`/`edtmax` = 0.1/1.
- **Shallow cloud model:** ε = 1e-3 m⁻¹, cap = 25 hPa, `xmbmax` = 0.1.

## 6. Output fields (`cu_grell_freitas_in` package)

- **Tendencies (`tend_physics`):** `rthcuten`, `rqvcuten`, `rqccuten`, `rqicuten`, plus `rqvdynten`, `rucuten` and `rvcuten`. The last two stay zero for GF.
- **Precipitation:** `cuprec`, `raincv`, `rainc`, `i_rainc`.
- **GF diagnostics:**
  - cloud base and top indices: `cubot`, `cutop`;
  - cloud-base mass fluxes: `xmb_total` (deep), `xmb_shallow`;
  - shallow cloud levels: `k22_shallow`, `kbcon_shallow`, `ktop_shallow`;
  - `ktop_deep`, which GF never sets;
  - in-cloud condensate: `qc_cu`, `qi_cu`.

  ([Registry.xml:2867](../src/core_atmosphere/Registry.xml#L2867))

## 7. Caveats and open questions

- **`ierr = 1200` is undone.** `cup_forcing_ens_3d` resets any `ierr > 995` to 0 and sets aa0 = 0 ([module_cu_gf.mpas.F:2630](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L2630)). The updraft, downdraft and work-function calculations were skipped for that column while it was flagged. The closures and `cup_output_ens_3d` then run on partly initialized arrays, so it is not clear from the code that the "resolved ascent → GF off" test actually switches the scheme off. This path only matters at Δx ≲ 6 km.
- **Restricted `possible_values`.** Registry lists only 0 and 8 for the two closure options. Other values run, but they are untested in MPAS.
- **Timing mismatch.** With `config_conv_interval` longer than `dt`, GF computes its closure and limiters for one `dt_dyn` and the result is then held for the whole interval.
- **Stale and non-advective `rqvdynten`.** It comes from the previous step and includes all physics applied during that step, including convection itself, so it is not purely advective.
- **Shallow scheme at fine Δx.** Shallow convection is not scale-aware and is not suppressed when deep convection is active, so at convection-permitting Δx it stays fully active.
- **No momentum transport.** GF computes no convective momentum transport, even though `rucuten`/`rvcuten` are allocated for it.
- **Unused inputs.** `hfx`, `qfx`, `gsw`, `areaCell` and `ccn` are passed but have no effect in this configuration.
- **Minor code issues.** The `gdc`/`gdc2` zeroing loop uses `k = jts, kte` ([module_cu_gf.mpas.F:184](../src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F#L184)), and the driver passes `kms = kds`. Both are harmless in MPAS because jts = kts = kms = kds = 1.
