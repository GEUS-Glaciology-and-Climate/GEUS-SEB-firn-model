## Changes

Three parts:
- **Shortwave penetration** is now on by default, with realistic parameters.
- **Closure checks:** the model writes mass and energy closure errors, and the test-run report checks them.
- **Bug fixes**, several of them found by the new closure checks.

### 1. Shortwave penetration on by default, with realistic parameters
SW penetration (#16) sent **all** the net shortwave radiation into the snow and ice with Beer–Lambert extinction, with `ext_ice = 5` m⁻¹. That gave too much internal melt (×3.5 total melt at DY2 in #16), so it was off by default. Now:
- **`sw_surface_fraction = 0.36`** (new constant): the fraction of the net shortwave absorbed at the surface (mostly near infrared, absorbed in the top millimetres). The other 0.64 penetrates (Greuell and Konzelmann, 1994; Woods and Hewitt, 2023, *The Cryosphere* 17, 1967–1987, Table 1).
- **`ext_ice = 1.5` m⁻¹** (was 5): the value used by Woods and Hewitt (2023), following Taylor and Feltham (2005). They note it is at the lower end of the 0.9–8 m⁻¹ observed by Cooper et al. (2021) for 350–750 nm on the Greenland ice sheet. The code's own reference range was 0.6–1.5 m⁻¹ (Greuell and Maykut, 1977).
- **`ext_snow = 20` m⁻¹** (unchanged), in the 4–40 m⁻¹ range of Greuell and Maykut (1977).
- **`sw_penetration = 1`** by default.

The same surface fraction is used for snow and ice. For snow it is probably too low: Brandt and Warren (1993) find the solar heating concentrated in the top few millimetres. With 20 m⁻¹, though, the penetrating part is absorbed within the top ~10–20 cm.

**Effect** (see the tables below):
- **Total melt** (surface + internal) changes little: −3 % and −2 % at KAN_L, ±0.003 m w.e. at DY2.
- **Where it melts** changes: at KAN_L about half of the melt is now internal (1.4–2.2 m w.e.), since 64 % of the net shortwave penetrates into the bare ice and ice is only weakly extinguishing (whether this split is realistic still needs to be checked against observations); at DY2 it is 0.02–0.05 m w.e.
- **Runoff and SMB** change by less than 0.05 m w.e.

### 2. Closure checks
At every time step the model now writes three outputs (`output_variables_info.csv`):
- **`mass_error_mweq`:** the change of the column mass (snow bucket included) minus what entered or left it (snowfall, deposition and sublimation, rain, runoff, ice added at the column bottom);
- **`energy_error`:** the same for the column enthalpy (sensible heat of snow and ice relative to ice at 0 °C, plus latent heat of the liquid water), with melt energy, absorbed shortwave, rain, runoff, new and sublimated snow, and ice added at the bottom;
- **`ground_flux_error`:** the conductive flux used by the surface energy budget minus the flux received by the first layer in the heat diffusion.

`tools/test_runs/compare_runs.py` also shows the total melt (surface + internal) next to the surface and internal melt, in the table, and `*_cumulative_melt_smb.png` shows internal melt and internal + surface melt of both runs on the same panel (the gap between them is the surface melt). Internal melt is meltwater like surface melt and reaches the SMB through runoff. It also adds a closure table to `summary.md`, flags any run outside the tolerances (mass 1e-9 m w.e., energy and flux 1e-6 W m⁻²), and prints a warning. On the test runs all three are at round-off level (< 4e-14 m w.e., < 1e-9 W m⁻², < 4e-12 W m⁻²). For comparison, the bugs fixed in #19 and #20 gave 16–21 W m⁻² and 0.05 m w.e.

### 3. Bug fixes
1. **`variables_preparation`: `pres` used before it is defined.** When the forcing has relative humidity instead of specific humidity (AWS forcing), `SpecHumSat(RH, T, pres, c)` was called before `pres` was read: `NameError`. CARRA forcing has specific humidity, so it was not affected.
2. **Old-snow roughness length not capped.** Lefebre et al. (2003) go from `z0_old_snow` at 600 kg m⁻³ to `z0_ice` at 920 kg m⁻³. The density used (`rho[0, k]`) includes the liquid water, so wet snow can exceed 920 kg m⁻³ and get a roughness length larger than the ice one. It is now capped at `z0_ice`. This is the only change that affects the test runs (KAN_L_2012, < 0.0005 m w.e.).
3. **`ind_ice` for SW penetration.** It was `np.argmax(depth_m > snowthick)`, which returns 0 (all ice) when `snowthick` reaches the bottom of the column. That case became possible with the diagnosed `snowthick` of #20 at sites without ice layers. It is now the first ice layer of the column (`first_ice_layer`, the same criterion as `snowthick`). SW penetration is off by default, so there is no effect on the test runs.
4. **The "Problem closing energy budget" warning was never printed.** The test was `findbalance == c.iter_max_EB`, but `range(1, iter_max_EB)` stops at `iter_max_EB − 1`. It now prints when the iteration stopped without converging. In the 34-year runs at DY2, KAN_L, KAN_U and SWC this happens a few dozen times (out of 102 000 steps), and the unclosed imbalance counted as melt is about 10⁻⁴ m w.e. in total.
5. **Initial layers denser than ice.** The initial state used to be stored as snow (`snowc`, `rhofirn` = bulk density), except when the whole profile was exactly 900 kg m⁻³. Layers denser than the model's ice (e.g. up to 939 kg m⁻³ in the KAN_U profile) were kept as snow denser than ice: negative pore space, and a jump of the layer volume at the first densification step, which caps `rhofirn` at `rho_ice`. Now, in `IniVar`, every initial layer at or above `rho_ice` is ice (`snic`) with density `rho_ice`. The all-900 profile used at KAN_L (ablation profile) is converted as before. At KAN_U, one layer (2.55 m w.e.) is now ice.
6. **First time step after the initialization.** The heat diffusion coefficients, surface conductance and `snowthick` of the first step now come from the initial state, as at the end of every step (`update_tempdiff_params_opt`, `calc_snowthick_m`).
   - Before, a start from the initial-state csv files had `grndc = T`, `grndd = 0`: no heat exchange between the surface and the subsurface during the first step, and `snowthick = snowthick_ini`.
   - A start from a spin-up pickle used the saved coefficients without the matching surface conductance (#19), and `snowthick` was not set at all (`np.empty`).
   - This changes runoff at DY2_2012 by −0.001 m w.e. through the first steps; the other test cases change by less than 0.0005 m w.e.
7. **Merged layers gained heat from their liquid water.** Merging two layers (`merge_small_layers`, `merge_layer`, the merge in `sublimation`, and adding the snow bucket to the first layer) weighted the temperatures by the total mass, liquid water included. The water is at 0 °C and has no sensible heat relative to ice at 0 °C, so a cold layer merged with a wet one was warmed without any heat source: up to 185 kJ m⁻² (17 W m⁻²) in one step at KAN_L_2019. Found by the energy closure. The temperature of merged layers now conserves the enthalpy of the frozen mass (`mix_temperature`).
8. **One heat capacity throughout.** `refreeze`, `melting`, `sw_heating` and the merges used the heat capacity at the current temperature (`cpiceF(T)`) over temperature changes of up to ~10 K, while it depends on temperature. They now all use the enthalpy (the integral of `cpiceF`) and its exact inverse (`enthalpy_per_kg`, `temperature_from_enthalpy`). This is how the energy closes to round-off; before, single refreezing events were off by up to 68 kJ m⁻² (6 W m⁻²).
9. **Negative precipitation and negative melt.**
   - The CARRA forcing has tiny negative precipitation (~1e-8 m w.e.) from the interpolation. The subsurface model ignored it but the SMB output counted it; it is now clipped at 0 when the forcing is read.
   - Likewise, when the surface energy budget iteration does not converge, the remaining negative imbalance gave negative melt (a few dozen steps in 34 years, about 1e-4 m w.e.). It is now 0, and the imbalance is reported by the warning (fix 4).
10. **`snowthick_ini` removed** from the constants and `main_SEB_firn.py`: the initial snow thickness now comes from the initial column (fix 6).

With `sw_penetration = 0`, the bug fixes alone change melt, runoff, refreezing and SMB by at most 0.001 m w.e. on the test runs (attribution run), and the closure is also at round-off level.

## Effect on the test runs

Totals over each test case in m w.e. (fix/seb-small-fixes vs main). Constant Tdeep = 263 K, no spin-up (initial profiles from input/initial state).

| Test case | Variable | main | fix/seb-small-fixes | Difference | Difference (%) |
|---|---|---|---|---|---|
| DY2_2012 | Total melt | 0.704 | 0.707 | +0.003 | +0.4 |
| DY2_2012 | Surface melt | 0.704 | 0.661 | -0.042 | -6.0 |
| DY2_2012 | Internal melt | 0.000 | 0.046 | +0.046 | n/a |
| DY2_2012 | Runoff | 0.023 | 0.023 | -0.000 | -1.1 |
| DY2_2012 | Refreezing | 0.702 | 0.710 | +0.008 | +1.1 |
| DY2_2012 | Sublimation | 0.000 | 0.001 | +0.001 | +149.1 |
| DY2_2012 | SMB | 0.564 | 0.565 | +0.001 | +0.2 |
| DY2_2019 | Total melt | 0.265 | 0.262 | -0.003 | -1.0 |
| DY2_2019 | Surface melt | 0.265 | 0.247 | -0.018 | -7.0 |
| DY2_2019 | Internal melt | 0.000 | 0.016 | +0.016 | n/a |
| DY2_2019 | Runoff | 0.000 | 0.000 | +0.000 | n/a |
| DY2_2019 | Refreezing | 0.248 | 0.248 | +0.000 | +0.1 |
| DY2_2019 | Sublimation | 0.006 | 0.006 | +0.001 | +11.7 |
| DY2_2019 | SMB | 0.278 | 0.279 | +0.001 | +0.2 |
| KAN_L_2012 | Total melt | 4.420 | 4.297 | -0.123 | -2.8 |
| KAN_L_2012 | Surface melt | 4.420 | 2.905 | -1.515 | -34.3 |
| KAN_L_2012 | Internal melt | 0.000 | 1.391 | +1.391 | n/a |
| KAN_L_2012 | Runoff | 4.345 | 4.298 | -0.047 | -1.1 |
| KAN_L_2012 | Refreezing | 0.190 | 0.176 | -0.013 | -7.0 |
| KAN_L_2012 | Sublimation | 0.041 | 0.042 | +0.001 | +1.4 |
| KAN_L_2012 | SMB | -3.775 | -3.727 | +0.047 | +1.3 |
| KAN_L_2019 | Total melt | 4.688 | 4.609 | -0.079 | -1.7 |
| KAN_L_2019 | Surface melt | 4.688 | 2.423 | -2.265 | -48.3 |
| KAN_L_2019 | Internal melt | 0.000 | 2.186 | +2.186 | n/a |
| KAN_L_2019 | Runoff | 4.506 | 4.515 | +0.009 | +0.2 |
| KAN_L_2019 | Refreezing | 0.132 | 0.121 | -0.011 | -8.1 |
| KAN_L_2019 | Sublimation | -0.036 | -0.032 | +0.004 | +11.9 |
| KAN_L_2019 | SMB | -4.293 | -4.298 | -0.005 | -0.1 |


### Closure checks

Mass and energy closure of the subsurface column and agreement of the conductive flux between the surface energy budget and the subsurface (OK: within tolerance, FAIL: not, n/a: not written by this version).

| Test case | Check | Tolerance | main | fix/seb-small-fixes |
|---|---|---|---|---|
| DY2_2012 | Mass error, max per step | 1e-09 m w.e. | n/a | OK 3.7e-14 |
| DY2_2012 | Energy error, mean of abs | 1e-06 W m-2 | n/a | OK 2.5e-10 |
| DY2_2012 | Energy error, max per step | 1e-06 W m-2 | n/a | OK 5.5e-10 |
| DY2_2012 | SEB vs subsurface conductive flux, max | 1e-06 W m-2 | n/a | OK 4.9e-13 |
| DY2_2019 | Mass error, max per step | 1e-09 m w.e. | n/a | OK 2.8e-14 |
| DY2_2019 | Energy error, mean of abs | 1e-06 W m-2 | n/a | OK 2.5e-10 |
| DY2_2019 | Energy error, max per step | 1e-06 W m-2 | n/a | OK 6.2e-10 |
| DY2_2019 | SEB vs subsurface conductive flux, max | 1e-06 W m-2 | n/a | OK 4.9e-13 |
| KAN_L_2012 | Mass error, max per step | 1e-09 m w.e. | n/a | OK 3e-14 |
| KAN_L_2012 | Energy error, mean of abs | 1e-06 W m-2 | n/a | OK 2.7e-10 |
| KAN_L_2012 | Energy error, max per step | 1e-06 W m-2 | n/a | OK 4.6e-10 |
| KAN_L_2012 | SEB vs subsurface conductive flux, max | 1e-06 W m-2 | n/a | OK 3.5e-12 |
| KAN_L_2019 | Mass error, max per step | 1e-09 m w.e. | n/a | OK 2.9e-14 |
| KAN_L_2019 | Energy error, mean of abs | 1e-06 W m-2 | n/a | OK 2.6e-10 |
| KAN_L_2019 | Energy error, max per step | 1e-06 W m-2 | n/a | OK 4.7e-10 |
| KAN_L_2019 | SEB vs subsurface conductive flux, max | 1e-06 W m-2 | n/a | OK 2.8e-12 |

### DY2_2012

<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2012_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2012_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2012_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2012_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2012_lwc_profiles.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2012_lwc_profiles.png?raw=true" alt="">

### DY2_2019

<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2019_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2019_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2019_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2019_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2019_lwc_profiles.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/DY2_2019_lwc_profiles.png?raw=true" alt="">

### KAN_L_2012

<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2012_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2012_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2012_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2012_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2012_lwc_profiles.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2012_lwc_profiles.png?raw=true" alt="">

### KAN_L_2019

<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2019_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2019_cumulative_melt_smb.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2019_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2019_cumulative_runoff_refreezing.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/fix/seb-small-fixes/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2019_lwc_profiles.png?raw=true" alt="">
<img src="https://github.com/GEUS-Glaciology-and-Climate/GEUS-SEB-firn-model/blob/main/docs/PR_results/2026-11-10_seb-small-fixes/KAN_L_2019_lwc_profiles.png?raw=true" alt="">

🤖 Generated with [Claude Code](https://claude.com/claude-code)
