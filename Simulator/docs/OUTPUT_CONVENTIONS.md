# Output conventions

What the engine writes and how to read it: the meaning of each CSV column, and the naming and caching rules that govern the files themselves. 

## Output columns

| column                                            | meaning                                                                                                                        |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `Time`                                            | The recording instant `t`.                                                                                                     |
| `CGFCloningNormalization`                         | The growth-rate estimate `g^{(Z)}` of Eq. 26, i.e. `log Z_s / t`. Subtract the `s=0` run to get the CGF.                       |
| `MeanPopulationPotential`                         | The instantaneous average of the total potential, Eq. 24. Its time average is the growth-rate estimate `g^{(\Phi)}` of Eq. 26. |
| `AcceptedCurrentRate`                             | Mean walker current divided by `t`, the biased current density `j_k` of Eq. 25. Tends to `chi'(s)`.                            |
| `ReferenceCurrentRate`                            | The same for the untilted reference trajectory. No reference is evolved at `k = 0`, so this is `0` there.                      |
| `PotentialSpread`                                 | Largest deviation of any walker potential from the population mean, Eq. 27.                                                    |
| `EffectiveSampleSize`                             | Mean pre-cloning ESS (Eq. 22) over the recording window.                                                                       |
| `ResamplingSteps` / `ResamplingChecks`            | Number of resamplings performed / examinations made since `t=0`.                                                               |
| `AcceptedWalkerEvents` / `RejectedWalkerEvents` * | Number of accepted / rejected walker events, i.e. configuration jumps, since `t=0`.                                            |
| `DensityProfile`                                  | Average density profile across all walkers at the recording instant `t`, JSON vector.                                          |

* `RejectedWalkerEvents` are only non-zero for `uniformized` dynamics, where impossible configuration jumps are included deliberately.

All equations in this table refer to equations in `theory_algorithm_implementation.pdf`.

## Reading a result

Every column is a plain number except `DensityProfile`, which holds a JSON vector and has to be decoded:

```python
import json
import pandas as pd

rows = pd.read_csv('data/<name>.csv')
potential = rows['MeanPopulationPotential']       # one number per row
profile = rows['DensityProfile'].map(json.loads)  # one list of L densities per row
```

The counters -- `ResamplingSteps`, `ResamplingChecks`, `AcceptedWalkerEvents`, `RejectedWalkerEvents` -- are cumulative from the start of the run, so use `.diff()` for per-window counts. 

## Files, caching and restarts

A result is a CSV of rows and a JSON of metadata. The name is `<model>_<O|C>_<Mechanism>_<rates>_L..._M..._T..._E..._s..._k..._<boundary-parameters>_seed...`: the model family (`NNN`, `KLS` or `WASEP`), `O` for an open chain or `C` for a closed ring, the event mechanism (`direct`, `heap` or `uniformized`), then the rate description — `_eps..._del...` for KLS, the eight weights for NNN, nothing for WASEP. In `<boundary-parameters>`, closed runs carry `_rho...` (initial filling, which is conserved) whereas open runs carry the four reservoir rates followed by `_init...`, the initial configuration. A ring ignores the initial condition -- its particle number is already fixed by the filling -- so the tag is named only for open chains, exactly as in the cache signature.

The JSON layout is `flat-wasep-v1`, meaning every parameter of the run sits at the top level under its own descriptive key — `field_E`, `bias_s`, `alpha_left_injection`, `kls_epsilon` and so on. Three nested objects sit beside them: `cache_parameters`, what a later run must match to reuse and extend this one; `diagnostics`, the settings the engine filled in for whatever was left unset, along with a few statistics about the run performed; `checkpoint`, a base64 snapshot of the complete engine state. Re-running the same parameters with a larger `tmax` restores that state and simulates only the missing time, giving a result bitwise identical to the run done in one go; the shorter run it supersedes is deleted once the longer one has been written. Every parameter that affects the model takes part in the match; only `tmax` may differ. 