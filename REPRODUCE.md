# Reproduce — Forecasting ICANS after CAR-T (Amidi et al., JITC 2022)

Everything runs from the **committed de-identified** data — no download needed.
MATLAB (R2019b+; verified on R2026a). From the repo root:

```matlab
run_all        % Step 1: NTonset_PLR  -> LR_prob_allp_Lasso_LOO.mat (lasso LOO)
               % Step 2: Forecast_Matricies -> HMM forecast metrics (bias/MAD/wMAPE)
```

To rebuild the augmented-state file from the de-identified source first:
```matlab
build_data_new    % regenerates data_new.mat (Y/Ya/Z) from AllBWH_5days_deID.xlsx
```

| Paper item | Script | Input (committed) | Output |
|---|---|---|---|
| ICANS-grade trajectories `Y` (n=199) | `build_data_new.m` | `AllBWH_5days_deID.xlsx` | `data_new.mat` |
| Lasso-LR onset probabilities (LOO) | `NTonset_PLR.m` (`fcnGetData.m`) | `AllBWH_5days_deID.xlsx`, `data_new.mat`, `time_LR.mat` | `LR_prob_allp_Lasso_LOO.mat` |
| Forecast bias / MAD / wMAPE by horizon (Fig 4) | `Forecast_Matricies.m` | `LR_prob…mat`, `data_new.mat`, `AllBWH_5days_deID.xlsx` | metrics `Baias`,`MAD`,`wMAPE` |
| Time-to-end-of-ICANS forecast | `time_to_end.m` | same | figure |

The HMM transition/emission matrices come from `fcnEstimateTransitionAndEmissionMx.m`
+ `fcnGetMatrices_BaseCase.m`; forward prediction from `fcnForwardPredict.m`.

**De-identification is lossless:** the de-identified `AllBWH_5days_deID.xlsx`
reproduces the raw `Y` grade matrix exactly (max abs diff = 0, verified in MATLAB).
See `DATA_SOURCE.md`.
