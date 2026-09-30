# MYNN PBL and surface layer in MPAS-Atmosphere

## 1. Overview

MPAS runs the MYNN-EDMF boundary-layer scheme (`bl_mynn`) and the MYNN surface layer (`sf_mynn`) from NCAR's shared **MMM-physics** repository. Neither core is stored in this tree. The build fetches it into `src/core_atmosphere/physics/physics_mmm/` at tag `20250616-MPASv8.3` (commit a4baf7f) ([Externals.cfg](../src/core_atmosphere/Externals.cfg)). Only the thin MPAS wrappers live in `physics_wrf/`:
- [module_bl_mynn.F](../src/core_atmosphere/physics/physics_wrf/module_bl_mynn.F), [bl_mynn_pre.F](../src/core_atmosphere/physics/physics_wrf/bl_mynn_pre.F) and [bl_mynn_post.F](../src/core_atmosphere/physics/physics_wrf/bl_mynn_post.F);
- [module_sf_mynn.F](../src/core_atmosphere/physics/physics_wrf/module_sf_mynn.F) and [sf_mynn_pre.F](../src/core_atmosphere/physics/physics_wrf/sf_mynn_pre.F).

References to the core files below are plain `file:line` text into that tag: `bl_mynn.F90` (1244 lines), `bl_mynn_subroutines.F90` (6565), `sf_mynn.F90` (2237) and `mynn_shared.F90`.

The code carries no version string. The latest change-log entry is "v4.6 / CCPP" (`bl_mynn_subroutines.F90:309`). MPAS's own changelog says both schemes were updated to the WRF 4.6 source (L. Fowler, 2024).

**Development history** (`bl_mynn_subroutines.F90:77-175`):
- Original code by Nakanishi.
- F90/WRF port by Pagowski.
- Later development by Olson, Kenyon, Angevine, Suselj, Puhales, Fowler and others.

**Key references:**
- Nakanishi & Niino (2004, 2006, 2009).
- Olson et al. (2019, NOAA Tech. Memo OAR GSD-61).
- The mass-flux routine calls itself "experimental" (`bl_mynn_subroutines.F90:4635`).

| Component | Where | What it does |
|---|---|---|
| Eddy diffusivity (ED) | `mym_turbulence`, `mym_predict` | Level 2.5 (or 2.6/3.0) closure with prognostic TKE |
| Mass flux (MF) | `dmp_mf` | 8 dry/moist updraft plumes for nonlocal transport and shallow cumulus |
| Subgrid clouds | `mym_condensation`, plus the MF cloud step | cloud fraction and condensate (`cldfrac_bl`, `qc_bl`, `qi_bl`) and their buoyancy effect |
| Surface layer | `sf_mynn_run` | Monin–Obukhov fluxes, u*, exchange coefficients, 2 m and 10 m diagnostics |

**Selection.** `config_pbl_scheme = 'bl_mynn'` turns on the `bl_mynn_in` package ([mpas_atmphys_packages.F:165](../src/core_atmosphere/physics/mpas_atmphys_packages.F#L165)). MYNN is also the `convection_permitting` suite default, together with GF, Thompson and `sf_mynn` ([mpas_atmphys_control.F:170](../src/core_atmosphere/physics/mpas_atmphys_control.F#L170)). Choosing `bl_mynn` **silently forces `config_sfclayer_scheme = 'sf_mynn'`**, whatever you set, including `off` ([mpas_atmphys_control.F:374](../src/core_atmosphere/physics/mpas_atmphys_control.F#L374)). YSU is not allowed with `sf_mynn`, but `sf_mynn` can run with the PBL off.

## 2. Turbulence closure

- **Prognostic variable:** q² = 2·TKE (`qke`). The bottom boundary is a surface-similarity production term; the top value is held fixed. After the solve, qke is clipped to [1e-4, 150] m² s⁻² (`bl_mynn_subroutines.F90:2321`).
- **Closure levels.** Set by `config_mynn_closure`, tested as a real number (`bl_mynn_subroutines.F90:1857, 2353, 2402`):

| closure | qke | qsq (qw variance) | tsq, cov | counter-gradient terms |
|---|---|---|---|---|
| ≤ 2.5 (includes the listed 2.0) | prognostic | diagnostic | diagnostic | none |
| 2.5 < c < 3 (e.g. 2.6) | prognostic | prognostic | diagnostic | none |
| ≥ 3.0 | prognostic | prognostic | prognostic | γθ, γq, Sm correction |

  "2.0" is not a true Level-2 scheme: qke is prognostic at every setting.
- **Constants** (`bl_mynn_subroutines.F90:360-382`):
  - Set: B1 = 24, B2 = 15, G1 = 0.235, Pr = 0.74, C2 = 0.729, C3 = 0.340, C5 = 0.2. C3 is the Canuto/Kitamura value, not NN09's 0.352.
  - Derived: A1 = 1.18, A2 = 0.665, C1 = 0.137.
  - The Canuto/Kitamura modification (`CKmod = 1`) is hard-wired. It scales A2 by 1/(1 + max(Ri, 0)), which removes the critical Richardson number.
  - TKE diffuses with 3·Km; the variances diffuse with 1·Km.
- **Stability functions.**
  - When q² is below its Level-2 equilibrium, the Level-2 Sm and Sh are used, scaled by √(q²/q²_eq) (Helfand–Labraga).
  - Sh is clipped to [0, 4] and Pr to ≤ 5.
  - A background floor applies Sm, Sh ≥ 0.03·min(10·a·w, 1) inside MF plumes and ≥ 0.05·cf in cloud (`bl_mynn_subroutines.F90:2007`).
  - Km = L·q·Sm and Kh = L·q·Sh.

## 3. Mixing length (`config_mynn_mixlength`)

All options blend a surface length `els` (κz with stability corrections), a turbulent length `elt` (∝ ∫q·z dz / ∫q dz, clipped to 10–400 m) and a buoyancy length `elb` (∝ q/N).

| Option | Description | Scale-aware? |
|---|---|---|
| 0 | Original NN09 (`bl_mynn_subroutines.F90:880`). Labelled "+ BouLac", but BouLac is never called. | no; ignores `config_mynn_scaleaware` |
| 1 | RAP/HRRR nonlocal form. Above PBLH, blended with 0.3 × the BouLac length. | L × Psig_bl |
| **2 (default)** | Ito et al. (2015)-style, mostly local. Uses `elb`/`elf` from τ·√TKE, with τ from w*, clipped to 30–150 s (stable) or 50–200 s (unstable); the MF a·w enters `elb`. Above PBLH, blends to a free-atmosphere length. | L = Psig_bl·L + (1 − Psig_bl)·L_LES, with L_LES = min(els/(1 + els/12 m), elb) |

## 4. Subgrid clouds (`config_mynn_cloudpdf`)

`mym_condensation` runs *before* the turbulence update each step, so it uses the previous step's variances, L and Sh.

| Option | σ(s) from | Cloud fraction and condensate |
|---|---|---|
| 0 | prognostic/diagnostic qsq, tsq, cov | Gaussian PDF (Sommeria–Deardorff / NN09) |
| 1 | L, Sh and gradients (Kuwano-Yoshida 2010) | Gaussian |
| **2 (default)** | √qsq, capped at 0.666·qsat, floored at 2.5–4 % of qsat, inflated for dz > 100 m | Chaboureau–Bechtold (2002): cf = 0.5 + 0.36·atan(1.8(Q1 + 0.2)). Q1 is boosted above the PBL where resolved condensate exists. Clouds are zeroed above a diagnosed tropopause. |
| −1, −2 | as 1, 2 | buoyancy effects kept, but output cloud zeroed |

- **Phase split:** liquid vs ice by the linear ramp (T − 240)/(269 − 240).
- **Buoyancy feedback:** the functions `vt` and `vq` (Bechtold & Siebesma 1998 non-Gaussian factor) enter N² and GH, so clouds change Sm, Sh, the mixing length and TKE buoyancy production.
- **MF cloud step.** Where plumes condense and the stratus fraction is below 0.5, the plume cloud *overwrites* `cldfrac_bl` and `qc_bl` (it can lower them). It puts all plume condensate into `qc_bl`, even when cold.
- **Unsupported values.** A value such as 3 matches no CASE and leaves the cloud outputs undefined (`bl_mynn_subroutines.F90:2647-2957`).

**These clouds are diagnostic only in MPAS.**
- `icloud_bl = 0` is a hard-wired parameter ([mpas_atmphys_vars.F:418](../src/core_atmosphere/physics/mpas_atmphys_vars.F#L418)).
- Nothing outside the PBL driver reads `cldfrac_bl`, `qc_bl` or `qi_bl`, so radiation never sees them. `icloud_bl` does nothing inside the scheme either.

## 5. Mass flux (EDMF) and scale awareness

- **Trigger** (`bl_mynn_subroutines.F90:4998`): all three must hold:
  1. surface buoyancy flux > 0.002 K m s⁻¹;
  2. maximum plume width > 300 m;
  3. superadiabatic lowest ~50 m.
- **Plumes:** 8 plumes with diameters from 300 m up to `maxwidth`. maxwidth = min(1.2Δx, 1000 m, 1.1·PBLH, a cloud-base limit, a flux-dependent width).
- **Areas:** N ∝ l^−1.9, with total area ≤ 0.1. The area is reduced for weak fluxes and tapered to zero as the lowest-level wind goes from 10 to 25 m s⁻¹.
- **Initial conditions:** w, θl and qt excesses come from surface-layer similarity (σw, σθ, σq).
- **Entrainment:** ε = 0.33/(w·l) (Tian & Kuang 2016 style), with a floor of 3e-4 m⁻¹, plus an increase above PBLH + 1.5 km.
- **Vertical velocity:** dw/dz = −2εw + bB/w, with w ≤ 3 m s⁻¹.
- **Plume microphysics:** saturation adjustment with no precipitation.
- **Plume top:** each plume stops where w ≤ 0. If *any* plume fails at the first level, the whole column's MF is turned off (`bl_mynn_subroutines.F90:5248`).
- **Flux limiter:** fluxes are rescaled when the MF heat flux near the surface would exceed 0.75 × the surface flux.
- **Options:**
  - `config_mynn_edmf_mom`: MF momentum transport.
  - `config_mynn_edmf_tke`: MF TKE transport.
  - `config_mynn_mixscalars`: MF and diffusive transport of hydrometeor and aerosol numbers.
  - `config_mynn_edmf_dd`: Stratocumulus-only downdrafts (`ddmf_jpl`). Experimental; see §10.
- **Environmental subsidence and detrainment are off.** They are off by a hard-wired `env_subs = .false.` (`bl_mynn_subroutines.F90:407`), so the `sub_thl`, `sub_qv`, `det_thl` and `det_qv` outputs are always zero.

**Scale awareness** (`config_mynn_scaleaware = 1`, `scale_aware`, `bl_mynn_subroutines.F90:6273`). Two Shin & Hong (2013)-type partition functions are used; both are 1 at coarse Δx.
- **ED part:** Psig_bl, with x = 2.5Δx / min(PBLH, 3000). It scales the mixing length only (options 1 and 2).
- **MF part:** Psig_shcu, with x = 2.5Δx / min(PBLH + 500, 3500). It multiplies all MF fluxes and areas.
- **Resolved-updraft shut-off.** The MF is further reduced where the resolved |w| below PBLH + 500 m exceeds 1 m s⁻¹, and switched off at 2 m s⁻¹.

| Δx | Psig_bl (PBLH 1 km) | Psig_shcu (1 km) | Psig_bl (2 km) | Psig_shcu (2 km) |
|---|---|---|---|---|
| 100 m | 0.66 | 0.29 | 0.41 | 0.19 |
| 250 m | 0.92 | 0.58 | 0.74 | 0.40 |
| 500 m | 0.99 | 0.81 | 0.92 | 0.65 |
| 1 km | 1.00 | 0.94 | 0.99 | 0.85 |
| 3 km | 1.00 | 0.99 | 1.00 | 0.98 |

At typical MPAS mesh spacings (≥ 3 km), scale awareness is effectively inactive. As elsewhere in MPAS physics, Δx = `config_len_disp` / meshDensity^¼ ([mpas_atmphys_driver_pbl.F:452](../src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F#L452)).

## 6. Tendencies, PBL height, other options

- **Solver** (`mynn_tendencies`): one density-weighted implicit tridiagonal solve per variable, with the MF term M(φ_u − φ̄) treated semi-implicitly. Surface fluxes enter as the lower boundary; momentum gets an implicit drag ρu*²/wspd.
- **Variables mixed:** u, v, θl, and qv and qc (or qt if `config_mynn_mixqt = 1`, followed by saturation adjustment). Also qi (diffusion only, no MF), and, with `mixclouds` and `mixscalars`, the numbers nc, ni, nwfa and nifa.
  - **Snow is never mixed.** It is hard-disabled (`cloudmix .AND. .false.`, `bl_mynn_subroutines.F90:3576`), and the wrapper passes zeros (`bl_mynn.F90:992`).
  - **Ozone is mixed, then its tendency is zeroed** (`bl_mynn.F90:1122`).
  - All moisture is handled internally as specific humidity; the MPAS pre/post wrappers convert to and from mixing ratio.
- **`config_mynn_dheat_opt = 1`** (default): heating from TKE dissipation, q³/(B1·L)/cp, capped at 0.002 K s⁻¹ and tapered off aloft.
- **`config_mynn_topdown = 1`:** an extra TKE source from cloud-top radiative cooling. It uses `rthraten` (LW + SW) near the PBL top when cloud is present.
- **`config_mynn_stfunc`:** only changes the φm used in the *surface TKE production* boundary term. 0 selects Kansas/Dyer–Hicks, 1 selects Cheng–Brutsaert/Grachev. It does not change the surface-layer fluxes, which come from `sf_mynn`.
- **`config_mynn_tkebudget = 1`:** outputs `qshear`, `qbuoy`, `qdiss`, `qwt` and `dqke`.
- **`config_mynn_edmf_output = 1`:** copies out the `edmf_*` diagnostics. The fields are allocated either way, because their package is `bl_mynn_in`.
- **PBL height** (`get_pblh`):
  1. The θv-increase height: the first level where θv exceeds the minimum below 200 m by 1.0 K (water) or 1.25 K (land). The comments say 1.5 K.
  2. The height where TKE falls to max(qke_sfc/40, 0.02), kept within ±350 m of the θv height.
  3. The two are blended with weight tanh((zi − 200)/400). kpbl is the level containing PBLH.

## 7. MYNN surface layer (`sf_mynn`)

Order of operations in `sf_mynn_run` (`sf_mynn.F90:283-935`):

1. **Gustiness.** A convective velocity w* (Beljaars 1995) uses the previous step's fluxes and PBLH. A Mahrt–Sun subgrid term adds 0.32·max(Δx/5 km − 1, 0)^0.33. Both are added to the wind speed, with wspd ≥ 0.1 m s⁻¹.
2. **Bulk Richardson number** Rib at the lowest level, clipped to ±4.
3. **Roughness lengths:**
   - **Water:** COARE 3.0, with z0 from a Charnock coefficient rising 0.011 → 0.018 between 10 and 18 m s⁻¹ plus a smooth-flow term; zt = zq = 5.5e-5·Re*^−0.6.
   - **Land:** z0 from the LSM/land-use table. zt = zq from Zilitinkevich (1995) with C_zil = 0.085, or from Andreas (2002) where snow depth ≥ 0.1 m.
   - **Not configurable in MPAS.** `isftcflx = 0` and `iz0tlnd = 0` are hard-wired parameters in the driver ([mpas_atmphys_driver_sfclayer.F:31](../src/core_atmosphere/physics/mpas_atmphys_driver_sfclayer.F#L31)), and COARE 3.0 vs 3.5 is fixed in the core. The driver comments ("Charnock and Carlson-Boland") are leftovers from the old sfclay scheme.
4. **z/L:**
   - First guess from the Li et al. (2010) regression, then fixed-point iteration (`zolrib`, ≤ 19 iterations), clipped to ±20.
   - MPAS passes `itimestep = initflag` (0 or 1), so the Li et al. first guess is used every step and the u*/T* first-guess branch is dead.
5. **Stability functions:**
   - Stable: Cheng & Brutsaert (2005).
   - Unstable: a COARE blend of Dyer–Hicks and Grachev et al. (2000) free-convection forms.
   - Evaluated from lookup tables for |z/L| ≤ 10.
6. **u\*:** under-relaxed, u* = ½u*_old + ½·κ·wspd/Ψm, so it lags sudden changes. Floor of 0.005 over land.
7. **Fluxes and coefficients.** hfx and qfx are computed; hfx ≥ −250 W m⁻² over land and qfx ≥ −0.02. Also computed: ch, chs, chs2, cqs2, cd, ck, and qsfc (recomputed over water).
8. **Diagnostics:** u10 and v10, with special handling when the lowest level is < 13 m, and t2, q2 and th2. **When an LSM is on, `atmphys_sfc_diagnostics` overwrites t2/q2/th2** (at every cell for Noah), so `sf_mynn`'s 2 m values survive only with `config_lsm_scheme = 'off'`.

**Fractional sea ice** (`config_frac_seaice`, threshold 0.02). `sf_mynn` is called twice: once treating ice cells as land, and once as open water (xland = 2, z0 reset to 1e-4, T = max(SST, 271.4 K)). The results are blended by ice fraction ([mpas_atmphys_driver_sfclayer.F:1016](../src/core_atmosphere/physics/mpas_atmphys_driver_sfclayer.F#L1016)). `znt` and `regime` are not blended.

**What the PBL actually uses from `sf_mynn`.** `bl_mynn` uses `ust`, `wspd` (gustiness included; it is the denominator of the surface drag), `hfx`, `qfx` and `tsk`. It recomputes 1/L from the fluxes every step. `ch`, `qsfc` and `znt` are passed in but never used, and `rmol` is used only at initialization.

## 8. Coupling to MPAS

### 8.1 Call sequence and timing

**Order within a step.** The physics driver runs radiation (on its alarm), then the surface layer, then the LSM/sea ice (which overwrite land hfx/qfx/qsfc), then **MYNN PBL**, then GWDO, then convection. The dynamics then applies the stored tendencies on every RK stage, and microphysics runs at the end of the step.

**When the PBL runs.** It is called every step whenever both the PBL and surface-layer schemes are not `off` ([mpas_atmphys_driver.F:320](../src/core_atmosphere/physics/mpas_atmphys_driver.F#L320)).

**`config_pbl_interval` is broken for MYNN (and YSU).** A non-`none` value sets `dt_pbl` and creates an alarm ([mpas_atmphys_manager.F:598](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L598)), but nothing ever checks that alarm. The PBL still runs **every step**, and its implicit solve and TKE update use the longer `dt_pbl`.

**Initialization.** `bl_mynn_init` runs at start-up ([mpas_atmphys_driver_pbl.F:746](../src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F#L746)). On the first step of a cold start (`initflag = 1`), `mym_initialize` iterates qke, tsq, qsq, cov and el to Level-2 equilibrium.

### 8.2 Inputs

| MYNN input | Source |
|---|---|
| u, v, w, θ, T, p, Π, ρ, dz, qv, qc, qi, qs, ni (and nc, nifa, nwfa with `mp_thompson_aerosols`) | model state via `mpas_atmphys_interface` |
| ust, wspd, znt, ch, qsfc, rmol | `sf_mynn` (blended over sea ice) |
| hfx, qfx, tsk | surface layer over water; LSM over land; blended by `driver_seaice` |
| `rthraten` = LW + SW | last radiation call (used only by top-down mixing and the downdrafts) |
| qke, tsq, qsq, cov, el_pbl, qke_adv, sh3d, sm3d, cloud fields | MYNN's own state from the previous call (`diag_physics`) |
| uoce, voce | set to **0** ([mpas_atmphys_driver_pbl.F:369](../src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F#L369)); ocean currents are not used |

### 8.3 Outputs and consumers

| Output | Consumer |
|---|---|
| `rublten`, `rvblten` | dynamics, interpolated to edges. GWDO adds its drag into these same arrays ([mpas_atmphys_driver_gwdo.F:680](../src/core_atmosphere/physics/mpas_atmphys_driver_gwdo.F#L680)) |
| `rthblten`, `rqvblten`, `rqcblten`, `rqiblten` | dynamics; GF (deep forcing and all of the shallow forcing); new Tiedtke |
| `rqsblten` (always 0), `rniblten`, and `rncblten`/`rnifablten`/`rnwfablten` (aerosol Thompson only) | dynamics ([mpas_atmphys_todynamics.F:406](../src/core_atmosphere/physics/mpas_atmphys_todynamics.F#L406)) |
| `hpbl`, `kpbl` | the next step's surface layer (w*); GWDO (`bl_ysu_gwdo`, `bl_ugwp_gwdo`); GF shallow source level and cap |
| `kzh`, `kzm` (exchange coefficients) | output only |
| `qke`, `tke_pbl`, `el_pbl`, `cldfrac_bl`, `qc_bl`, `qi_bl`, `edmf_*`, budget terms | output and restart state only |

### 8.4 Interactions with other parameterizations

- **Microphysics.**
  - The `bl_mynn_in` package itself enables the `qc`, `qi`, `qs` and `ni` scalars (Registry.xml, moist and number arrays). With WSM6 or Kessler, `ni` (and for Kessler also `qi` and `qs`) exist as passive fields that MYNN mixes but no microphysics uses.
  - With Thompson, MYNN mixes qc, qi and ni. With `mp_thompson_aerosols` it also mixes nc, nwfa and nifa.
  - Snow is never mixed, and MYNN's subgrid clouds are not passed to the microphysics.
- **Radiation.** Radiation influences MYNN through `rthraten` (top-down option only). MYNN has no path back to radiation (`icloud_bl = 0`): radiation sees only the resolved qc and qi, which MYNN has mixed.
- **Convection (GF).**
  - GF runs *after* MYNN and uses `rthblten`/`rqvblten` as forcing. The default shallow closure (BLQE) is driven entirely by them, and `kpbl` sets GF's shallow source level.
  - MYNN's EDMF shallow-cumulus plumes and GF's shallow scheme both run, with no coordination, so shallow cumulus can be represented twice.
- **Surface layer / LSM.** The surface layer uses the previous step's MYNN PBLH for its w* gustiness. The LSM supplies land fluxes. MYNN recomputes stability from those fluxes rather than taking the surface layer's z/L.
- **GWDO.** It reads `kpbl`/`hpbl` and writes into MYNN's `rublten`/`rvblten`/`rthblten`.

## 9. Namelist options

All of these are in `&physics`. Registry.xml lines 2599–2682 give defaults and allowed values. Integer switches are converted to logicals as `== 1` ([module_bl_mynn.F:309-335](../src/core_atmosphere/physics/physics_wrf/module_bl_mynn.F#L309)), so any other value means off.

| Option | Default | Allowed | Effect |
|---|---|---|---|
| `config_pbl_scheme` | `suite` | `suite`, `bl_ysu`, `bl_mynn`, `off` | `bl_mynn` selects MYNN and forces `sf_mynn` |
| `config_sfclayer_scheme` | `suite` | `suite`, `sf_monin_obukhov`, (`sf_monin_obukhov_rev`), `sf_mynn`, `off` | overridden to `sf_mynn` when the PBL is `bl_mynn` |
| `config_physics_suite` | `mesoscale_reference` | `mesoscale_reference`, `convection_permitting`, `none` | `convection_permitting` gives `bl_mynn` + `sf_mynn` |
| `config_mynn_closure` | 2.5 | 2.0, 2.5, 3.0 (real; 2.6 also works) | closure level (§2). 2.0 behaves as 2.5 |
| `config_mynn_mixlength` | 2 | 0, 1, 2 | mixing-length formulation (§3) |
| `config_mynn_cloudpdf` | 2 | 0, 1, 2 (−1 and −2 also work) | subgrid-cloud PDF (§4) |
| `config_mynn_scaleaware` | 1 | 0, 1 | Psig_bl on the mixing length and Psig_shcu on the MF (§5) |
| `config_mynn_edmf` | 1 | 0, 1 | mass-flux plumes |
| `config_mynn_edmf_mom` | 0 | 0, 1 | MF momentum transport |
| `config_mynn_edmf_tke` | 0 | 0, 1 | MF TKE transport |
| `config_mynn_edmf_dd` | 0 | 0, 1 | Stratocumulus downdrafts (experimental, buggy; §10) |
| `config_mynn_edmf_output` | 0 | 0, 1 | fill the `edmf_*` diagnostics |
| `config_mynn_mixscalars` | 1 | 0, 1 | mix number concentrations and aerosols (with `mixclouds`) |
| `config_mynn_mixclouds` | 1 | 0, 1 | mix qc, qi, nc, ni. Snow is never mixed |
| `config_mynn_mixqt` | 0 | 0, 1 | mix total water plus saturation adjustment, instead of each species separately |
| `config_mynn_dheat_opt` | 1 | 0, 1 | TKE dissipative heating |
| `config_mynn_topdown` | 0 | 0, 1 | cloud-top radiatively driven TKE source |
| `config_mynn_stfunc` | 1 | 0, 1 | φm in the surface TKE boundary term: Dyer–Hicks or Cheng–Brutsaert |
| `config_mynn_tkeadvect` | false | logical | uses `qke_adv` in place of `qke`. **MPAS never advects `qke_adv`** (it is a diagnostic field, not a scalar), so this only reloads the previous value |
| `config_mynn_tkebudget` | 0 | 0, 1 | TKE-budget output |
| `config_pbl_interval` | `none` | `DD_HH:MM:SS` | **broken**: changes `dt_pbl` but not the call frequency (§8.1) |
| `config_len_disp` (`&nhyd_model`) | 0 → `nominalMinDc` | m | Δx for the scale-aware functions, plume widths and surface-layer gustiness |
| `config_frac_seaice` | true | logical | two-call ice/water surface layer (§7) |
| `config_do_DAcycling` | false | logical | intended to keep TKE across DA cycles, but see §10 |

**Hard-wired, not in the namelist:**
- `icloud_bl = 0` and `spp_pbl = 0` (no stochastic perturbations) ([mpas_atmphys_vars.F:417](../src/core_atmosphere/physics/mpas_atmphys_vars.F#L417)).
- Surface-layer `isftcflx = 0`, `iz0tlnd = 0` and `isfflx = 1`.
- Inside the MMM code: `env_subs = .false.`, `CKmod = 1`, 8 plumes, `Atot = 0.1`, `pgfac = 0`.

## 10. Caveats and possible bugs

These come from reading the code; none was confirmed in a run.

1. **Dew and frost fluxes are thrown away** (`bl_mynn_subroutines.F90:3471-3474`). The intended limit, `qvflux = max(qvflux, min(0.9·qv1 − 1e-8, 0)/dtz)`, reduces to `max(qvflux, 0)` whenever qv1 > ~1e-8. Every negative surface moisture flux is therefore set to zero in the default path (`mixqt = 0`). The LSM still books the deposition, so water is not conserved at night and in frost conditions.
2. **Wrong columns in `sf_mynn` with more than one OpenMP thread.** The wrapper passes the full `(ims:ime, kms:kme, jms:jme)` arrays to `sf_mynn_pre`, whose dummies are `(its:ite, 1:kte)` ([sf_mynn_pre.F:74](../src/core_atmosphere/physics/physics_wrf/sf_mynn_pre.F#L74), [module_sf_mynn.F:247](../src/core_atmosphere/physics/physics_wrf/module_sf_mynn.F#L247)). For any thread whose block does not start at cell 1, the lowest-level fields come from the wrong cells. Runs with one thread per MPI task are unaffected.
3. **`config_pbl_interval`** changes the time step but not the call frequency (§8.1).
4. **`config_mynn_tkeadvect` does no advection** in MPAS (§9).
5. **DA cycling TKE is zeroed.** In the cold-start block the test `.not.restart .or. .not.cycling` is always true, so qke is zeroed even when cycling asked to keep it (`bl_mynn.F90:699-721`).
6. **EDMF downdrafts (`edmf_dd = 1`):**
   - `ddmf_jpl` receives the whole 2-D `rthraten` instead of the column (`bl_mynn.F90:932`), so it reads the wrong elements.
   - It reads out of bounds (`DOWNA(K-1)` at K = kts), prints unconditionally, and uses svp1 = 0.6112 where 0.61 is meant.
   - Its momentum source has the opposite sign from the scalar terms.
   - Treat it as unusable.
7. **Density-weighted EDMF means are not divided by ρ.** `edmf_qt`, `edmf_thl` and `edmf_w` are Σρ·a·φ/Σa (`bl_mynn_subroutines.F90:5466`). The resulting factor of ~1.0–1.2 feeds the MF cloud σq (via QTp − q̄t), `qc_bl`, and the a·w terms in the mixing length.
8. **Heating counted twice for negative condensate.** `moisture_check` is passed θl and adds latent heating for any negative-condensate fix; the θ tendency then adds it again (`bl_mynn_subroutines.F90:3978-3996`).
9. **Units.** The surface θl flux uses hfx/(ρcp) without 1/Π, and dissipative heating is also added without 1/Π. The plume θv omits 1/Π on the condensate term, and the environment θv has no condensate loading while the plume θv does.
10. **Snow conversion typo.** MPAS `bl_mynn_pre.F` converts snow with qs/(1 + qs) instead of qs/(1 + qv) ([bl_mynn_pre.F:135](../src/core_atmosphere/physics/physics_wrf/bl_mynn_pre.F#L135)). This is harmless in practice because snow is not mixed.
11. **Dead or ignored settings:**
    - `mixlength = 0` has no BouLac and no scale awareness.
    - `icloud_bl` does nothing, and the `sub_*`/`det_*` outputs are always zero.
    - Ozone and snow are not mixed.
    - `ch`, `qsfc` and `znt` are passed to the PBL but unused.
    - Surface-layer `wstar` and `qstar` are discarded by the wrapper.
    - The Registry lists `rmol` as unitless; it is m⁻¹.
12. **Fragile MF switch.** One weak plume failing at level 2 turns off the MF for the whole column. With `cloudpdf = 2`, if no tropopause is found, all subgrid clouds above level 2 are removed.
13. **Doubled shallow convection.** MYNN EDMF shallow cumulus and GF shallow convection are not coordinated (§8.4).
