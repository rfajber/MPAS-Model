# Radiative cloud fraction in MPAS-Atmosphere

## 1. Overview

MPAS diagnoses a 3-D cloud fraction, `cldfrac`, whose only physical use is in the radiation schemes. `driver_cloudiness` computes it immediately before the SW and LW radiation calls, and only on steps when radiation runs ([mpas_atmphys_driver.F:215](../src/core_atmosphere/physics/mpas_atmphys_driver.F#L215)). The code is in [mpas_atmphys_driver_cloudiness.F](../src/core_atmosphere/physics/mpas_atmphys_driver_cloudiness.F) and [module_mp_thompson_cldfra3.F](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F). No microphysics, PBL or convection scheme reads `cldfrac`.

| `config_radt_cld_scheme` | Method | Needs condensate? | Changes the condensate radiation sees? | Selected by |
|---|---|---|---|---|
| `cld_incidence` | 0/1: cloudy if qc + qi > 1e-6 | yes | no | explicit, or automatic fallback when radiation is on and the scheme is `off` |
| `cld_fraction` | Xu & Randall (1996) RH–condensate formula | yes | no | `mesoscale_reference` and `convection_permitting` suites |
| `cld_fraction_thompson` | Thompson: Sundqvist-type RH scheme with a grid-spacing-dependent critical RH | no | **yes**: adds "radiative" qc/qi to RH-only cloud layers | explicit only |
| `off` | none | — | — | `none` suite. Replaced by `cld_incidence` with a log message if LW or SW radiation is on ([mpas_atmphys_control.F:350](../src/core_atmosphere/physics/mpas_atmphys_control.F#L350)) |

Registry's `possible_values` omits `cld_fraction_thompson` ([Registry.xml:2564](../src/core_atmosphere/Registry.xml#L2564)), but the control check accepts it ([mpas_atmphys_control.F:340](../src/core_atmosphere/physics/mpas_atmphys_control.F#L340)).

**Inputs.** All schemes use the state at the start of the step, via `MPAS_to_physics`. The driver first copies qv, qc, qi and qs into separate "radiative" arrays `qvrad_p`, `qcrad_p`, `qirad_p` and `qsrad_p` ([mpas_atmphys_driver_cloudiness.F:134](../src/core_atmosphere/physics/mpas_atmphys_driver_cloudiness.F#L134)). Only the Thompson option modifies them; the model's prognostic water is never changed.

## 2. `cld_incidence`

`cldfrac = 1` where qc + qi > 1e-6 kg kg⁻¹, otherwise 0 ([mpas_atmphys_driver_cloudiness.F:237](../src/core_atmosphere/physics/mpas_atmphys_driver_cloudiness.F#L237)).
- **Warm-only microphysics.** If the microphysics has no cloud ice (`f_qi` false, e.g. Kessler), only qc is tested.
- **Snow is ignored** in the test, but RRTMG still counts snow mass in the ice path of cells that are cloudy.
- **Result:** every condensate-bearing layer is overcast, and grid-mean condensate is spread uniformly over the cell.

## 3. `cld_fraction` (Xu & Randall 1996)

Code: [mpas_atmphys_driver_cloudiness.F:284](../src/core_atmosphere/physics/mpas_atmphys_driver_cloudiness.F#L284).

**Saturation.**
- Saturation mixing ratios over water and over ice follow Murray (1966).
- They are blended by the ice fraction of the condensate, w = (qi + qs)/(qc + qi + qs), to give qvs = (1 − w)·qvsw + w·qvsi. RH = qv/qvs.
- With no condensate, w = 0, so saturation is taken over water.

**Cloud fraction.** With q_cld = qc + qi + qs (snow included):
- If q_cld < 1e-12: cf = 0. There is no cloud without condensate, whatever the RH.
- If RH ≥ 1: cf = 1.
- Otherwise:

  cf = RH^0.25 · [1 − exp(−100·q_cld / ((qvs − qv)^0.49))],

  with the exponent floored at −6.9, and cf < 0.01 reset to 0.
- The constants α₀ = 100, γ = 0.49 and p = 0.25 are the Xu–Randall values, hard-wired.

**What this means.** Cloud fraction rises with both condensate amount and RH. Small amounts in sub-saturated air give tiny fractions that are zeroed. With qvs = 10 g kg⁻¹, q_cld = 1e-5 gives cf ≈ 0.03 at RH 0.9 (kept) and ≈ 0.09 at RH 0.99. q_cld = 1e-6 gives < 0.01 below RH 0.99, so it is dropped.

## 4. `cld_fraction_thompson` (`cal_cldfra3`)

Code: [module_mp_thompson_cldfra3.F:44](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L44). It was copied from WRF 3.8.1. The header says it was tested only with Thompson microphysics, "should not be used with other cloud microphysics schemes", and is "not intended for combining with any cumulus or shallow cumulus parameterization".

### 4.1 Cloud fraction from RH

- **Overcast from condensate:** cf = 1 wherever qc > 1e-6, qi ≥ 1e-7 or qs > 1e-5 ([module_mp_thompson_cldfra3.F:101](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L101)).
- **Saturation reference:** otherwise, saturation is over water above −12 °C, over ice below −20 °C, and blended linearly in between.
- **Critical RH** depends on grid spacing Δx in km, with the higher value over ocean ([module_mp_thompson_cldfra3.F:90](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L90)):

  RH₀(land) = 0.781 + √(1/(35 + 0.5Δx³)),  RH₀(ocean) = 0.831 + √(1/(70 + 0.5Δx³)).

  Δx = `config_len_disp`/meshDensity^¼, converted to km. A 2024 bug fix added the m→km conversion ([mpas_atmphys_driver_cloudiness.F:124](../src/core_atmosphere/physics/mpas_atmphys_driver_cloudiness.F#L124)).

  | Δx (km) | 1 | 3 | 5 | 10 | 15 | 30 | 60 | 120 |
  |---|---|---|---|---|---|---|---|---|
  | RH₀ land | 0.95 | 0.93 | 0.88 | 0.82 | 0.81 | 0.79 | 0.78 | 0.78 |
  | RH₀ ocean | 0.95 | 0.94 | 0.92 | 0.87 | 0.86 | 0.84 | 0.83 | 0.83 |

- **T ≥ −12 °C (Sundqvist):** cf = 1 − √((1 − RH)/(1 − RH₀)), with RH ≤ 0.999.
- **−70 °C < T < −12 °C:** only where RH > RH₀(ocean), which is used even over land. Here cf = 1 − √((RH_max − RH)/(RH_max − RH₀)), with RH_max = qvsw/qvsi, so cf reaches 1 at water saturation.
- **T ≤ −70 °C:** no RH cloud.
- **Cap:** RH-derived fractions are limited to 0.9 ([module_mp_thompson_cldfra3.F:134](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L134)).

### 4.2 Column clean-up (`find_cloudLayers`)

1. **Tropopause.** Searching down from the model top, it is the first level between 4 and 19 km where dθ/dz < 10 K / 1.5 km, measured over two levels. Fractional clouds (0 < cf < 0.999) above it are removed.

   If no level qualifies, k_tropo = 3, so *all* fractional clouds above level 3 are removed ([module_mp_thompson_cldfra3.F:250](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L250)).
2. **Mixed layer.** Fractional clouds are removed in the well-mixed layer near the surface. That layer runs up to where dθ > 0.05 K km⁻¹·dz first appears, minus 2 levels ([module_mp_thompson_cldfra3.F:273](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L273)).

### 4.3 Adding radiative condensate

For RH-only clouds, the scheme adds condensate to `qcrad`/`qirad`. Only the radiation codes see it.
- **Condensate estimate:** max_wc = |qvs(top − 1) − qvs(base)| for each contiguous cloud layer of 2 or more levels, distributed in proportion to cumulative dz and scaled by (1 − entr), with entr = 0.5.
- **Ice layers** (between the −12 °C level and the tropopause): qi += 0.1·cf·iwc for 0.01 < cf < 0.99, or 0.01·iwc for overcast cells with little ice. T must be ≥ 203 K ([module_mp_thompson_cldfra3.F:419](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L419)).
- **Liquid layers** (from the −12 °C level down to the mixed-layer top, 253–298 K): qc += cf²·lwc for fractional cells, or 0.1·lwc for overcast cells with little qc ([module_mp_thompson_cldfra3.F:463](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L463)).
- **Single-level clouds** get a flat 1e-5·cf.
- **No additions where water already exists:** layers whose existing path exceeds 1 kg m⁻² are left alone.
- **Path cap:** the column path in fractional-cloud cells is normalized to at most 1 kg m⁻² ([module_mp_thompson_cldfra3.F:500](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L500)).

## 5. How radiation uses it

### 5.1 RRTMG LW and SW (McICA)

**In-cloud water paths.** Grid-mean paths are converted to in-cloud paths by dividing by max(0.01, cf) ([module_ra_rrtmg_lw.F:12141](../src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_lw.F#L12141)):
- liquid from `qcrad`;
- ice from `qirad` (plus `qsrad` if the ice-optics flag is 3);
- snow as a separate path only when microphysics effective radii are in use (`config_microp_re = true`).

The McICA sub-column generator then places cloud using cf and the overlap assumption. Two consequences:
- **Condensate in a cell with cf = 0 is radiatively invisible.** Examples: qc + qi ≤ 1e-6 under `cld_incidence`, cf < 0.01 under `cld_fraction`, and snow-only cells under `cld_incidence`.
- **For a fixed water path, a smaller cf makes a denser but smaller cloud.** The cloud fraction therefore changes the radiative effect even when the water is the same.

**Overlap** (`config_radt_cld_overlap` → `icld`, [mpas_atmphys_driver_radiation_lw.F:865](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_lw.F#L865)):

| Value | icld | Meaning |
|---|---|---|
| `none` | 0 | **clear sky only**: RRTMG ignores all clouds. It does not mean "no overlap assumption" |
| `random` | 1 | random |
| `maximum_random` | 2 | maximum–random (default) |
| `maximum` | 3 | maximum |
| `exponential` | 4 | exponential |
| `exponential_random` | 5 | exponential–random |

For the two exponential options, `config_radt_cld_dcorrlen` sets the decorrelation length:
- `constant`: 2.5 km.
- `latitude_varying`: the NASA GMAO latitude and day-of-year function. On this branch it follows `config_perpetual_julday` ([mpas_atmphys_driver_radiation_lw.F:843](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_lw.F#L843)).

The decorrelation length is ignored for the other overlap options.

**Effective radii.** By default (`config_microp_re = false`) RRTMG does **not** use microphysics effective radii. Liquid radius comes from a land/ocean/temperature formula (`relcalc`: 8 µm over land at T ≥ 0 °C, ramping to 14 µm by −20 °C; 14 µm over ocean, sea ice or snow), and ice from a temperature table (`reicalc`, Kristjansson–Mitchell). Snow mass is added to the ice path.

With `config_microp_re = true`, Thompson or WSM6, and RRTMG for both LW and SW ([mpas_atmphys_manager.F:815](../src/core_atmosphere/physics/mpas_atmphys_manager.F#L815)), radiation uses the microphysics re_cloud, re_ice and re_snow. Where cf > 0 but re is at its floor, which is typical of RH-only Thompson cloud whose condensate came from §4.3, defaults are used: liquid 10.5 µm (ocean) or 7.5 µm (land), and ice from the temperature table ([module_ra_rrtmg_lw.F:12072](../src/core_atmosphere/physics/physics_wrf/module_ra_rrtmg_lw.F#L12072)).

### 5.2 CAM LW and SW

CAM receives the same `cldfrac` but the **original** `qc_p`/`qi_p`, not the radiative copies ([mpas_atmphys_driver_radiation_lw.F:925](../src/core_atmosphere/physics/mpas_atmphys_driver_radiation_lw.F#L925)). With `cld_fraction_thompson`, CAM therefore sees the extra cloud fraction without the extra condensate. The overlap namelist options do not apply to CAM.

### 5.3 Diagnostics

`cldfrac_low_UPP`, `cldfrac_mid_UPP`, `cldfrac_high_UPP` and `cldfrac_tot_UPP` are the column maxima of `cldfrac` in the bands p ≥ 642 hPa, 350–642 hPa and 150–350 hPa, and over the whole column ([mpas_cloud_diagnostics.F:94](../src/core_atmosphere/diagnostics/mpas_cloud_diagnostics.F#L94)). That is maximum overlap, independent of the radiation overlap setting.

`cldfrac` itself is refreshed only on radiation steps, so between calls the output field is stale.

## 6. Interactions with other parameterizations

- **Microphysics.** The microphysics supplies qc, qi and qs, which drive every option.
  - With Thompson or WSM6 and `config_microp_re = true`, radiation also uses the microphysics effective radii. The default is off.
  - `cld_fraction_thompson` is documented as Thompson-only, but MPAS does not enforce this.
  - Under Kessler, `cld_incidence` tests qc only, and `cld_fraction` sees qi = qs = 0.
- **Convection.** No scheme contributes a convective cloud fraction. Parameterized convection affects cloud fraction only through the condensate it detrains into qc/qi:
  - GF detrains condensate into qc or qi; its in-cloud `qc_cu`/`qi_cu` are output only.
  - KF also adds qr and qs.
  - With `cld_fraction`, detrained condensate in dry air gives small or zero fractions.
  - The Thompson scheme header warns it is not designed to combine with cumulus schemes.
- **MYNN PBL.** MYNN computes its own subgrid cloud (`cldfrac_bl`, `qc_bl`, `qi_bl`), but `icloud_bl = 0` is hard-wired, so it never reaches `cldfrac` or radiation. Subgrid PBL clouds (stratocumulus, shallow cumulus) affect radiation only through `cld_fraction_thompson`'s RH-based cloud, or through the resolved qc that MYNN mixes.
- **Radiation timing.** Radiation calls follow `config_radtlw_interval`/`config_radtsw_interval`, so cloud fraction is only as fresh as the last radiation step.

## 7. Namelist options

| Option | Default | Allowed values | Effect |
|---|---|---|---|
| `config_radt_cld_scheme` | `suite` → `cld_fraction` | `cld_incidence`, `cld_fraction`, `cld_fraction_thompson` (missing from the Registry list), `off` | cloud-fraction method. `off` becomes `cld_incidence` when radiation is on |
| `config_physics_suite` | `mesoscale_reference` | `mesoscale_reference`, `convection_permitting`, `none` | both real suites use `cld_fraction`; `none` gives `off` |
| `config_radt_cld_overlap` | `maximum_random` | `none`, `random`, `maximum_random`, `maximum`, `exponential`, `exponential_random` | RRTMG McICA overlap. **`none` removes clouds from radiation** |
| `config_radt_cld_dcorrlen` | `constant` | `constant` (2.5 km), `latitude_varying` | decorrelation length, used only with the exponential overlaps |
| `config_radtlw_interval`, `config_radtsw_interval` | `00:30:00` | `DD_HH:MM:SS` or `none` | how often the cloud fraction is recomputed |
| `config_len_disp` (`&nhyd_model`) | 0 → `nominalMinDc` | m | Δx for the Thompson critical RH |
| `config_microp_re` | false | logical | use microphysics effective radii (Thompson or WSM6) in RRTMG. It needs RRTMG for both LW and SW |

**Hard-wired, not in the namelist:**
- `cld_incidence` threshold 1e-6.
- Xu–Randall α₀ = 100, γ = 0.49, p = 0.25, and the 0.01 cut-off.
- Thompson condensate thresholds, the 0.9 cap, entr = 0.5 and the RH₀ coefficients.

## 8. Caveats

1. **`config_radt_cld_overlap = 'none'` means clear-sky radiation,** not a neutral overlap choice.
2. **Condensate can be invisible to radiation.** Anything below the thresholds in §5.1 is dropped. Under `cld_fraction`, clouds cannot form in saturated air without condensate. Under `cld_fraction_thompson`, the reverse happens: RH-only clouds are given made-up condensate.
3. **Thompson radiative condensate is inconsistent with the model water.** The extra qc/qi exist only in RRTMG, and CAM does not get them at all (§5.2).
4. **Fragile tropopause test** in `cld_fraction_thompson`. If no tropopause is found, all fractional clouds above level 3 are removed. The same fallback exists in MYNN's cloud PDF.
5. **Ocean threshold used over land.** In the −12 to −70 °C branch, RH₀(ocean) is applied everywhere, even over land ([module_mp_thompson_cldfra3.F:129](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L129)).
6. **The 2016 MPAS "bounds fix" changed behaviour.** `k = min(k_m12C, k_m12C+2)` always equals k_m12C, so the liquid-layer search now starts two levels lower than in WRF ([module_mp_thompson_cldfra3.F:326](../src/core_atmosphere/physics/physics_wrf/module_mp_thompson_cldfra3.F#L326)).
7. **Stale cloud base for single-level clouds.** In `find_cloudLayers`, a level with no cloud falls into the "single-level cloud" branch using the previous `k_cldb`. The 1e-5·cf top-up can therefore be applied at the last cloud base found rather than at the current level. This only adds water where it is below 1e-6.
8. **Cloud radiative effects from subgrid schemes are missing.** MYNN's `cldfrac_bl` and GF's in-cloud water are never used, so a convection-permitting-suite run represents shallow and PBL cloud radiatively only through the resolved qc/qi.
9. **Stale output.** `cldfrac` is not updated between radiation calls.
