# RRTMG inputs in MPAS-Atmosphere

This note traces what MPAS passes into RRTMG LW and SW, and where each field comes from. It does not describe the radiative transfer itself. Emphasis is on clouds and aerosols.

## 1. Call sequence

On a radiation step (`config_radtlw_interval`/`config_radtsw_interval`, default 30 min), `physics_driver` does the following, in order ([mpas_atmphys_driver.F:207-244](../src/core_atmosphere/physics/mpas_atmphys_driver.F#L207)):

1. **`MPAS_to_physics`** ([mpas_atmphys_interface.F](../src/core_atmosphere/physics/mpas_atmphys_interface.F)) converts the model state at the start of the step into WRF-style column arrays.
2. **`driver_cloudiness`** diagnoses `cldfrac` and makes the "radiative" water copies `qvrad`, `qcrad`, `qirad` and `qsrad` (see [Cloud_fraction_summary.md](Cloud_fraction_summary.md)).
3. **`driver_radiation_sw` → `radiation_sw_from_MPAS` → `rrtmg_swrad`** ([mpas_atmphys_driver_radiation_sw.F:933](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_sw.F#L933)).
4. **`driver_radiation_lw` → `radiation_lw_from_MPAS` → `rrtmg_lwrad`** ([mpas_atmphys_driver_radiation_lw.F:879](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_lw.F#L879)).

Before the step, the physics manager refreshes the ozone climatology to the current (or perpetual) day whenever RRTMG will be called ([mpas_atmphys_manager.F:377](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L377)).

**What happens between radiation calls.** RRTMG is not called at all. `rthratenlw`/`rthratensw` are reapplied unchanged on every step, and so are the surface fluxes `gsw`/`glw` given to the LSM. The SW heating is **not rescaled for the changing solar zenith angle** between calls.

## 2. Atmospheric state

Everything comes from `MPAS_to_physics`, at time level 1, before this step's dynamics and microphysics.

| RRTMG argument | MPAS source | Notes |
|---|---|---|
| `t3d` | T = θ·Π, with θ = θ_m/(1 + (R_v/R_d)·qv) | layer temperature |
| `t8w` | `t2_p`: T at w-levels | linear in z. Extrapolated at the surface and model top, as in WRF `phy_prep` |
| `p3d`, `p8w` | `pres_hyd_p`, `pres2_hyd_p` | **Hydrostatic** pressure, integrated downward from an extrapolated model-top pressure using ρ_d(1 + qv) ([mpas_atmphys_interface.F:534](../src/core_atmosphere/physics/mpas_atmphys_interface.F#L534)). There is no condensate loading, and it is not the model's nonhydrostatic pressure |
| `pi3d` | Exner function | converts the heating rate to a θ tendency |
| `dz8w` | layer thickness from `zgrid` | |
| `qv3d` | `qvrad_p` = qv (≥ 0) | converted to H₂O vmr |
| `xlat`, `xlong` | `latCell`, `lonCell` | also used for the latitude-varying decorrelation length |

**Above the model top.** RRTMG adds buffer layers of about 4 hPa each ([module_ra_rrtmg_lw.F:11837](../src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_lw.F#L11837)):
- **Temperature:** from a standard profile, shifted to match the model top.
- **H₂O and trace gases:** held at their top-level values.
- **Ozone:** from the ozone source (§5).
- **Clouds:** none.

## 3. Clouds and hydrometeors

| RRTMG argument | Source | Set by |
|---|---|---|
| `cldfra3d` | `cldfrac` | `driver_cloudiness` (`config_radt_cld_scheme`) |
| `qc3d` | `qcrad_p` = qc, plus the Thompson cloud-fraction additions if `cld_fraction_thompson` | microphysics, PBL and convection (via the model qc), and `cal_cldfra3` |
| `qi3d` | `qirad_p` = qi, plus the same Thompson additions | same |
| `qs3d` | `qsrad_p` = qs | microphysics |
| `re_cloud`, `re_ice`, `re_snow` | `re_cloud`, `re_ice`, `re_snow` in `diag_physics` | microphysics, only when `has_reqc/i/s = 1` |
| `icloud` | 1 | hard-wired in [mpas_atmphys_vars.F](../src/core_atmosphere/physics/mpas_atmphys_vars.F) |
| `cldovrlp`, `idcor` | `config_radt_cld_overlap`, `config_radt_cld_dcorrlen` | namelist |

**Not passed to RRTMG:**
- **Rain (qr) and graupel (qg).** Their argument slots exist in the LW wrapper but MPAS does not use them, so rain and graupel are radiatively invisible.
- **Subgrid cloud from other schemes.** MYNN's `cldfrac_bl`/`qc_bl`/`qi_bl` (`icloud_bl = 0`) and GF's `qc_cu`/`qi_cu` never reach RRTMG. Convection and the PBL affect radiation only through the resolved qc and qi they produce or mix.

**How the cloud fields are used.**
- **In-cloud paths.** RRTMG builds grid-mean liquid, ice and snow paths from qc, qi and qs, divides them by max(0.01, cf), and gives them to McICA with the overlap choice. Condensate in a cell with cf = 0 has no effect, and `config_radt_cld_overlap = 'none'` removes all clouds.
- **Optics options.** With defaults, `inflg = 2`, `iceflg = 3` (Fu ice optics) and `liqflg = 1` (Hu & Stamnes). Snow mass is lumped into the ice path.
- **Effective radii.** Controlled by `config_microp_re`, default **false**, so the microphysics radii are *not* used by default ([mpas_atmphys_manager.F:815](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L815)):

  | Case | Liquid radius | Ice radius | Snow |
  |---|---|---|---|
  | default | `relcalc`: 8 µm over warm land, rising to 14 µm by −20 °C; 14 µm over ocean, sea ice or snow cover | `reicalc`: Kristjansson–Mitchell temperature table (×1.0315, ≤ 140 µm) | lumped into ice |
  | `config_microp_re = true`, Thompson or WSM6, both LW and SW RRTMG | `re_cloud` from the microphysics (computed at the end of the previous step) | `re_ice` | separate path with `re_snow` (`iceflg` 4/5) |

  In the second case, floor values are replaced by defaults: 10.5 µm (ocean) or 7.5 µm (land), and the ice table. Changing `re` changes optical depth at fixed water path, so this switch has a first-order effect on cloud radiative forcing.
- **Where the Thompson aerosol-aware scheme enters the clouds.** With `mp_thompson_aerosols`, the prognostic nc and nwfa change the droplet sizes and hence `re_cloud`. Those radii reach radiation only if `config_microp_re = true`. This is the only path by which aerosol affects radiation through the clouds; there is no separate aerosol indirect-effect code in the radiation.

## 4. Aerosols (direct effect)

**LW:** no aerosols at all. `rrtmg_lwrad` has no aerosol arguments.

**SW:** `aer_opt` is set in `radiation_sw_from_MPAS` ([mpas_atmphys_driver_radiation_sw.F:531](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_sw.F#L531)):

| Microphysics | `aer_opt` | What RRTMG SW receives |
|---|---|---|
| anything except `mp_thompson_aerosols` | 0 | tauaer = 0, ssa = 1, g = 0: **no aerosol** |
| `mp_thompson_aerosols` | 3 | Layer and band aerosol optical properties diagnosed from Thompson's prognostic aerosol numbers (see below) |

**The `aer_opt = 3` chain:**
1. **AOD from aerosol numbers.** `gt_aod` ([module_mp_thompson_aerosols.F](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_aerosols.F)) computes the 550 nm AOD per layer:

   τ(k) = [b_wfa(RH, T)·nwfa + b_ifa(RH, T)·nifa]·ρ_a·Δz

   `nwfa` is the "water-friendly" and `nifa` the "ice-friendly" aerosol number, both in kg⁻¹. The unit extinctions come from an 8-RH × 4-T lookup table, with RH clipped to 10–98 %.
2. **Column AOD.** The column sum is saved to `taod5502d`, and the layer values to `taod5503d`, for output.
3. **Spectral properties.** `calc_aerosol_rrtmg_sw` ([module_ra_rrtmg_sw_aerosols.F](../src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_sw_aerosols.F)) spreads the 550 nm layer AOD across the 14 SW bands. Using the Ruiz-Arias WRF parameterization, it takes the Ångström exponent, single-scattering albedo and asymmetry parameter from RH and an aerosol type.
   - Type: `xland = 1` (land and sea ice) → "rural"; `xland = 2` → "maritime". There is no urban type.
   - Option settings are hard-wired ([mpas_atmphys_vars.F:667](../src/core_atmosphere/physics/mpas_atmphys_vars.F#L667)): `taer_aod550_opt = 2` (use the gridded AOD) and `angexp`/`ssa`/`asy` option 3 (from RH and type). The constant fallbacks (AOD 0.12, Ångström exponent 1.3, SSA 0.85, g = 0.9) are therefore unused.
4. **Into RRTMG.** The arrays are passed as `tauaer3d`, `ssaaer3d` and `asyaer3d`, and RRTMG uses them with `iaer = 10`.

**Where nwfa and nifa come from.** They are prognostic Thompson scalars:
- initialized from the climatology in the initial conditions;
- advected by the dynamics;
- mixed by the MYNN PBL;
- changed by microphysical activation and scavenging, and by a surface emission term in the microphysics.

The SW aerosol effect therefore responds to the model's own aerosol evolution. On this branch nwfa and nifa can also be nudged (commit eaf50068).

**The `aerosols` state array** (a CAM climatology) is fetched by both drivers but used **only by CAM radiation**, not by RRTMG.

## 5. Gases

| Gas | Source |
|---|---|
| H₂O | `qvrad` (model qv) |
| O₃ | `config_o3climatology = true` (default): the monthly CAM ozone climatology (`ozmixm`, read at init), interpolated to the current day (`o3clim`, [mpas_atmphys_manager.F:377](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L377)) and then to model and buffer-layer pressures (`vinterp_ozn`). With `false`, RRTMG's built-in annual-mean profile (`inirad`) is used |
| CO₂ | 379 ppmv (IPCC 2005) by default ([module_ra_rrtmg_lw.F:11680](../src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_lw.F#L11680)), or `config_co2vmr` when `config_fixed_co2 = true` (this branch) |
| CH₄, N₂O, O₂, CFC-11, CFC-12, CFC-22, CCl₄ | constant 2005 values: CH₄ 1774 ppb, N₂O 319 ppb, O₂ 0.209488, CFCs in ppt. Not configurable |

All gases except H₂O and O₃ are vertically uniform.

## 6. Surface and solar inputs

| RRTMG argument | Source | Set by |
|---|---|---|
| `tsk` | `skintemp` | LSM over land; SST update, slab ocean or sea-ice scheme over water/ice |
| `albedo` (SW) | `sfc_albedo` | Land-use table at init, updated by Noah/Noah-MP, the sea-ice driver and `update_surface` (0.08 open water, 0.80 ice, blended by ice fraction; `config_slab_albedo` on this branch). **One broadband value is used for all four of RRTMG's direct/diffuse × UV-vis/near-IR albedos** ([module_ra_rrtmg_sw.F:10394](../src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_sw.F#L10394)), and ocean albedo has no zenith-angle dependence |
| `emiss` (LW) | `sfc_emiss` | same routines: 0.98 over water, table values over land |
| `xland`, `xice`, `snow` | `sfc_input` | land/water mask, ice fraction, snow water. Used for aerosol type and default droplet radii |
| `coszr`, `declin`, `solcon`, `gmt`, `xtime`, `julday` | `radconst` and the clock | Solar geometry. On this branch `config_perpetual_julday` can lock the declination and solar constant to a fixed day (commit 395aebd8). `julday` is used by McICA only for the latitude-varying decorrelation length |

## 7. Outputs and their consumers

| Output | Consumer |
|---|---|
| `rthratenlw`, `rthratensw` | dynamics θ tendency every step ([mpas_atmphys_todynamics.F:472](../src/core_atmosphere/physics/mpas_atmphys_todynamics.F#L472)); GF deep forcing; MYNN top-down mixing |
| `glw`, `gsw`, `swddir`, `swddni`, `swddif` | LSM, sea ice and slab ocean surface energy budget; surface layer (via LSM fluxes) |
| TOA and surface fluxes (all-sky and clear-sky), `olrtoa`, `lwcf`, `swcf`, `taod5502d`/`taod5503d` | diagnostics only |
| `rre_cloud`, `rre_ice`, `rre_snow` | diagnostics only: the radii RRTMG actually used |

## 8. Namelist options affecting RRTMG inputs

| Option | Default | Effect on inputs |
|---|---|---|
| `config_radt_lw_scheme`, `config_radt_sw_scheme` | `suite` → `rrtmg_lw`/`rrtmg_sw` | selects RRTMG (or CAM, or off) |
| `config_radtlw_interval`, `config_radtsw_interval` | `00:30:00` | call frequency; tendencies are held in between |
| `config_radt_cld_scheme` | `suite` → `cld_fraction` | how `cldfrac` and the radiative qc/qi are built |
| `config_radt_cld_overlap` | `maximum_random` | McICA overlap. `none` means clear sky |
| `config_radt_cld_dcorrlen` | `constant` | decorrelation length for the exponential overlaps |
| `config_microp_re` | `false` | use microphysics effective radii (Thompson or WSM6, both RRTMG) |
| `config_microp_scheme` | `suite` | `mp_thompson_aerosols` is the only way to get SW aerosols (and aerosol-aware radii) |
| `config_o3climatology` | `true` | climatological ozone vs. RRTMG's default profile |
| `config_fixed_co2`, `config_co2vmr` | `false`, 379e-6 | constant CO₂ override (this branch) |
| `config_perpetual_julday` (this branch) | off | fixed-day solar geometry and ozone date |

Hard-wired:
- the trace-gas values other than CO₂;
- `icloud = 1`;
- the RRTMG optics flags;
- the aerosol option settings and type mapping;
- the buffer-layer spacing (4 hPa).

## 9. Caveats

1. **Default radii ignore the microphysics.** With `config_microp_re = false`, cloud radiative properties depend only on water path, cloud fraction and a temperature/land formula for the radii. This holds even with Thompson, and even with the aerosol-aware Thompson scheme.
2. **No SW aerosol without `mp_thompson_aerosols`,** and no LW aerosol at all.
3. **Rain and graupel are never seen by RRTMG.** Snow is seen, lumped into the ice path by default.
4. **Subgrid clouds from MYNN and GF are not seen** (§3).
5. **Radiation lags the state.** It uses the state at the start of a radiation step and holds the heating for up to one radiation interval, with no zenith-angle rescaling of the SW heating in between.
6. **Pressure differs from the dynamics.** The radiation pressure is a hydrostatic reconstruction, not the model's nonhydrostatic pressure.
7. **Broadband albedo.** A single value is used for all spectral and direct/diffuse components.
8. **Cloud fraction and radiative water copies** can be modified by `cld_fraction_thompson` independently of the model's water (see [Cloud_fraction_summary.md](Cloud_fraction_summary.md)).
