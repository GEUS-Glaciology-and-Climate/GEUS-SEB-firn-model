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
