# Thompson microphysics in MPAS-Atmosphere

## 1. Overview

MPAS runs a port of WRF's Thompson scheme (`module_mp_thompson.F`, now referenced to WRF v4.1.4). The module header cites Thompson et al. (2008, MWR 136) for the base scheme and Thompson & Eidhammer (2014, JAS 71) for the aerosol-aware version. Both modes run the same code (`mp_gt_driver` → `mp_thompson`), and one module logical, `is_aerosol_aware`, switches between them ([module_mp_thompson.F:86](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L86)).

| `config_microp_scheme` | Mode | How it is selected | Package |
|---|---|---|---|
| `mp_thompson` | non-aerosol (WRF mp_physics=8) | `convection_permitting` suite, or set explicitly | `mp_thompson_in` |
| `mp_thompson_aerosols` | aerosol-aware (WRF mp_physics=28) | set explicitly only (no suite uses it) | `mp_thompson_aers_in` |

The flag is `.true.` only after `init_thompson_aerosols_forMPAS` runs ([mpas_atmphys_init_microphysics.F:133](../src/core_atmosphere/physics/mpas_atmphys_init_microphysics.F#L133)), and that routine is called only for `mp_thompson_aerosols`. In both modes the scheme is sub-stepped `n_microp = max(nint(dt_dyn/90),1)` times per call.

## 2. Common framework

Number concentrations are passed per kg and multiplied by ρ inside the scheme to get m⁻³.

| Species | Mass | Number, `mp_thompson` | Number, `_aerosols` | Size distribution |
|---|---|---|---|---|
| cloud water | qc | none; nc = Nt_c | nc, capped at 1999e6 m⁻³ | gamma, integer ν_c = MIN(15, NINT(1000e6/nc)+2); mean D 1–100 µm |
| rain | qr | nr | nr | exponential spheres; mvd 37.5 µm–2.5 mm |
| cloud ice | qi | ni | ni | exponential spheres, ρ_i = 890; mean D 5–300 µm |
| snow | qs | — | — | Field et al. (2005) two-gamma sum, m = 0.069 D² |
| graupel | qg | — | — | exponential, ρ_g = 500, N0 diagnosed |
| aerosol | — | fixed (nwfa 11.1e6, nifa 5e3 m⁻³) | nwfa, nifa | — |

- **Nt_c.** In MPAS this is a per-cell field, not a constant: 300e6 m⁻³ where landmask = 1 and 100e6 m⁻³ where landmask = 0. That gives ν_c = 5 over land and 12 over ocean ([mpas_atmphys_init_microphysics.F:81](../src/core_atmosphere/physics/mpas_atmphys_init_microphysics.F#L81)). `mu_c` is also computed and passed in, but the scheme never reads it.
- **Snow moments.** M_n = a(n,T_c)·M₂^b(n,T_c), using Field's cubic fits ([module_mp_thompson.F:1824](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L1824)).
- **Graupel intercept.** N0 = 10^zans1. zans1 is an empirical function of r_g and of the mvd of supercooled rain above the 270.65 K level. N0 is clamped to 1e4–3e6 m⁻⁴ and cannot increase downward. The code gives no reference for the formula ([module_mp_thompson.F:1909](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L1909)).
- **Fall speeds.** v = aD^b e^(−fD) √(ρ₀/ρ), with rain (4854, 1, 195; Ferrier 1994), snow (40, 0.55, 100), graupel (442, 0.89), ice (1847.5, 1) and cloud (0.317e8 D²).
- **Thresholds.** A species is treated as absent at ≤ R1 = 1e-12. New ice crystals have mass 1e-12 kg (≈12.9 µm).
- **Lookup tables.** Collection efficiencies are computed at init: Ef_rw from Beard & Grover (1974) and Pruppacher & Klett fits, Ef_sw from Wang & Ji (2000). Four `.DBL` files from `build_tables` are read: rain–graupel collection, rain–snow collection, Bigg freezing (1–45 K supercooling × IN from 1 to 1e6 m⁻³), and ice→snow ([module_mp_thompson.F:847](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L847)). A missing file is fatal. None of these tables depends on the mode.

## 3. The standard (non-aerosol) scheme

**Fixed inputs.** On every call, each column gets nc = Nt_c, nwfa = 11.1e6 m⁻³ and nifa = 5e3 m⁻³, and none of them is written back ([module_mp_thompson.F:1168](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L1168)). nc is forced back to Nt_c at lines 1679, 2871 and 3124, and again in `calc_effectRad`.

1. **Warm rain.**
   - Autoconversion follows Berry & Reinhardt (1974), using droplet-gamma characteristic diameters, for r_c > 0.01 g m⁻³ ([module_mp_thompson.F:1971](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L1971)). With nc fixed, it depends only on r_c and on land vs. ocean.
   - Accretion uses the tabulated Ef_rw.
   - Self-collection uses Ef_rr = 1 − exp(2300(mvd_r − 1.95 mm)). This turns into breakup for drops above 1.95 mm (Seifert 1994; Verlinde & Cotton 1993).
2. **Ice initiation (T < 0 °C).**
   - *Drop freezing* reads the Bigg tables at a fixed 1 IN L⁻¹, with no temperature shift ([module_mp_thompson.F:2310](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L2310)). Frozen raindrops lighter than a 250 µm graupel particle become ice and heavier ones become graupel.
   - *Deposition nucleation* triggers when s_i ≥ 0.25, or at water supersaturation below 253.15 K. It tops ice up to the Cooper (1986) curve, min(250e3, 5·e^(0.304(273.15−T))) m⁻³, with 1e-12 kg per crystal ([module_mp_thompson.F:2368](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L2368)).
   - *Hallett–Mossop* makes 3.5e8 splinters per kg of graupel rime, between −3 and −8 °C, peaking at −5 °C.
   - *Homogeneous freezing.* Aqueous aerosols do not freeze. Below HGFR = 235.16 K, cloud or rain amounts too small for the tables freeze outright. After sedimentation, all cloud water below HGFR freezes and all cloud ice above 0 °C melts ([module_mp_thompson.F:3546](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L3546)).
3. **Ice-phase growth and conversion.**
   - *Ice deposition* uses the Srivastava & Coen (1992) prefactor with capacitance 0.5. A fraction P(μ_i+2, λD0s) stays as ice and the rest goes to snow (Harrington et al. 1995) ([module_mp_thompson.F:2395](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L2395)).
   - *Ice → snow.* Ice larger than 200 µm converts to snow, and 99% of the mass converts if the mean diameter exceeds 1 mm.
   - *Snow deposition* uses a capacitance that rises from 0.15 at −1.5 °C to 0.5 at −30 °C.
   - *Graupel* only sublimates; it never grows by deposition.
   - *Riming.* Snow uses Ef_sw from the table. Graupel's Ef_gw is computed inline from the Stokes number, with no reference given. When riming exceeds twice deposition, 15–95% of snow rime goes to graupel and snow fall speed is boosted by up to 1.5×.
   - *Collection.* Snow collects ice with Ef = 0.05. Rain collecting ice (Ef = 0.95) produces graupel. Rain–snow collection (0.95) and rain–graupel collection (0.75) come from explicit-integral tables.
4. **Melting.** Snow and graupel melt with ventilation, and snow melting is enhanced by collected liquid. Above 0 °C, rain collecting graupel adds 5·tnr_gacr drops by breakup; Mansell's size-dependent variant is commented out. When dt > 120 s, riming above 0 °C goes to rain instead.
5. **Limiters.**
   - The ice terms that consume vapor are scaled to 0.999 of the vapor excess over ice saturation.
   - Each species' sinks are scaled to the mass available. Number rates are not scaled by either limiter.
   - The mean-size bounds are then enforced through the number tendencies.
6. **Condensation** uses a 3-iteration Newton saturation adjustment ([module_mp_thompson.F:3027](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L3027)).
7. **Rain evaporation.**
   - Srivastava & Coen (1992) with ventilation, capped at (q_vs − q_v)/dt. Number is lost in proportion to mass.
   - Under melting graupel, the rate is multiplied by 0.01 + 0.98 T_c/20 (marked "TEST", 2013).
8. **Sedimentation.**
   - Explicit upwind with sub-steps. Rain and ice move mass and number, with the number fall speed weighted by D^(b/2). Snow and graupel move mass only.
   - Cloud water sediments only in the lowest ~500 m where w < 0.1 m s⁻¹. What falls out of the lowest level is not counted as precipitation ([module_mp_thompson.F:3276](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L3276)).

## 4. What the aerosol-aware scheme changes

**Species.** nwfa is water-friendly aerosol, treated as sulfate-like CCN with κ = 0.4. nifa is ice-friendly dust.

**Initialization and sources.**
- *Cold start.* nwfa and nifa come from the GOCART monthly climatology, which `init_atmosphere` interpolates. If a field is exactly zero everywhere, the model uses an exponential height profile instead (50e6 + 300e6·e^(−z/h) for nwfa; 0.5e6/1.5e6 for nifa) ([mpas_atmphys_init_microphysics.F:171](../src/core_atmosphere/physics/mpas_atmphys_init_microphysics.F#L171)).
- *Restart.* Aerosol fields come from the restart file.
- *Surface source.* nwfa2d = nwfa(1)·0.000196·airmass·0.5e-10 is set once, at cold start, and nifa2d = 0. It is added to the lowest level after each call ([module_mp_thompson.F:1208](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L1208)). The code comment calls this "very far from ideal".

**Droplet activation.** Activation happens on net condensation.
- *Table lookup.* `activ_ncloud(T, w, nwfa)` interpolates the CCN_ACTIVATE_DATA activated fraction bilinearly in ln N and ln w. The table comes from a parcel model (Feingold & Heymsfield, modified by Eidhammer). Temperature uses the nearest 10 K entry, and radius (0.04) and κ (0.4) are fixed ([module_mp_thompson.F:4381](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L4381)).
- *Vertical velocity.* w is the grid-scale value, with a floor of 0.011 m s⁻¹, so droplets activate even in subsiding air.
- *Update.* nc is raised to the activated number if that is larger, not incremented, and nwfa loses the same amount ([module_mp_thompson.F:3049](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L3049)).

**Ice from nifa.** `iceDeMott` implements DeMott et al. (2010): N_IN = 5.94e-5 (−T_c)^3.33 · n_d^(−0.0264T_c+0.0033), where n_d is dust in cm⁻³, floored at 0.5 ([module_mp_thompson.F:4662](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L4662)). The Phillips (2008) branch is commented out, so this applies at any saturation. It is used in two places:
- It picks the Bigg freezing-table row, shifting drop freezing by ±3 K.
- It replaces Cooper as the deposition-nucleation target. Each crystal nucleated this way consumes one nifa.

**Aqueous homogeneous freezing.** This follows Koop et al. (2001), with J computed from the water-activity difference. It needs T < 238 K, s_i ≥ 0.4 and total ice-phase number ≤ 999e3 m⁻³. Aerosol radius is 0.025 µm and each new crystal is 1e-13 kg. It consumes nwfa ([module_mp_thompson.F:2383](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L2383)).

**Droplet number.** Autoconversion, accretion, riming and freezing all deplete nc. On partial evaporation, the `tnc_wev` table removes droplets smaller than the size that evaporates completely in one step.

**Aerosol budget.**
- *Sinks.* nwfa is removed by activation and Koop freezing, nifa by DeMott nucleation. Rain, snow and graupel scavenge both. Scavenging uses Eff_aero (Wang et al. 2010, after Slinn 1983), which includes Brownian diffusion, interception and impaction, with diameters of 0.04 µm (nwfa) and 0.8 µm (nifa) ([module_mp_thompson.F:2678](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L2678)).
- *Sources.* Each evaporating cloud drop or raindrop returns one nwfa. Nothing is returned by collection, freezing or sublimation, and nifa is never returned at all.
- *Transport.* Aerosols do not sediment. They are advected as scalars and mixed by MYNN, but not by YSU or convection. On this branch they can also be nudged (`config_nudging_tau_nwfa/nifa`).

| | `mp_thompson` | `mp_thompson_aerosols` |
|---|---|---|
| Droplet number | Nt_c = 100/300 cm⁻³ (ocean/land) | prognostic, activated from nwfa |
| ν_c | 12 / 5 | 3–15 |
| Drop-freezing IN | 1 L⁻¹, no shift | DeMott(nifa), ±3 K |
| Deposition nucleation | Cooper, capped at 250 L⁻¹ | DeMott (2010), consumes nifa |
| Aqueous homogeneous freezing | none | Koop (2001), consumes nwfa |
| Droplets lost to evaporation | only on complete evaporation | `tnc_wev` table; returned to nwfa |
| re_cloud | from Nt_c | from nc |
| SW aerosol | none (τ = 0) | `gt_aod` → RRTMG SW |

## 5. Coupling to MPAS

- **State update.** Microphysics updates the state directly.
  - Masses and θ are written back up to `mp_top_level` (default 45 km), and condensate is zeroed above it. Number fields are written on all levels.
  - `r*mpten` diagnostic tendencies are (after − before)/dt_dyn.
- **Precipitation.** RAINNC accumulates rain + snow + graupel + ice, SNOWNC snow + ice, and GRAUPELNC graupel ([module_mp_thompson.F:1190](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L1190)).
- **Effective radii.** RRTMG uses the scheme's radii only if `config_microp_re = .true.` (the default is false) and both LW and SW are RRTMG ([mpas_atmphys_manager.F:814](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L814)). Clips: cloud 2.49–50 µm, ice 4.99–125 µm, snow 9.99–999 µm.
- **AOD.** Aerosol mode with `rrtmg_sw` only.
  - `gt_aod` turns nwfa, nifa and RH into a 550 nm AOD.
  - `calc_aerosol_rrtmg_sw` splits it into bands using Shettle–Fenn rural (land) or maritime tables ([mpas_atmphys_driver_radiation_sw.F:530](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_sw.F#L530)).
  - Longwave radiation gets no aerosol.
- **Cloud fraction and reflectivity.** `cld_fraction_thompson` and `calc_refl10cm` take no nc or aerosol input, so they behave the same in both modes.
- **Initial and boundary data.** `init_atmosphere` supplies only qv, qc, qr, nifa and nwfa. There is no nc, ni, nr, qi, qs or qg in either the initial conditions or the LBCs.

## 6. Caveats and open questions

- **Integer bin index.** `nic1` is declared INTEGER ([module_mp_thompson.F:256](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L256)), so 7.93 truncates to 7. For nc = 100 and 300 cm⁻³ the lookups then use the ≈189 and ≈682 cm⁻³ bins, which biases the cloud-freezing and evaporation tables in both modes.
- **Final-clamp units.** The final update clamps per-kg nwfa/nifa with m⁻³ bounds, and caps nc1d at Nt_c_max, which is in m⁻³ ([module_mp_thompson.F:3581](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L3581)).
- **Rain–snow table.** `qr_acr_qs` inverts the bm_s = 2 moment test relative to the run-time code ([module_mp_thompson.F:3791](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson.F#L3791)). The size of the effect on the table was not assessed.
- **DeMott floor.** The 0.5 cm⁻³ dust floor is about 100× the nifa floor, so in clean air the floor, not nifa, sets the IN concentration.
- **Resolution-dependent emission.** nwfa2d scales with cell area, so the emission rate changes with mesh resolution.
- **Uninitialized fields.** Cold-start nc and `lbc_nc`/`lbc_ni`/`lbc_nr` are never supplied. What values they end up with was not traced.
- **Scaled mass, unscaled number.** The vapor limiter scales nucleated mass but not number, so aerosol consumption uses the unscaled number.
- **Melting number loss.** The instant melt above 0 °C drops the number of ice that sedimented in during the step.
- **Evaporation factor.** The factor under melting graupel has no lower bound and could turn negative just below 0 °C.
- **Stale metadata.** The Registry unit for nt_c (nb kg⁻¹) is wrong; the code uses m⁻³. `mu_c`, `L_nwfa` and `L_nifa` are dead. The comments mention a Koop J reduction that is not visible in the code, and they cite Bigg as both 1953 and 1954.
- **Table age.** Whether a run's `.DBL` files predate the 2017 `freezeH2O` change cannot be told from the source.