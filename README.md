# ICANS-forecasting-after-CAR-T-cell-therapy

MATLAB code and de-identified data for:

> Amidi Y, Eckhardt CA, Quadri SA, Malik P, Firme MS, Jones DK, Jain A, Danish HH,
> Rubin DB, Jacobson CA, Cash SS, Lee JW, Dietrich J, **Westover MB**.
> *Forecasting immune effector cell-associated neurotoxicity syndrome after chimeric
> antigen receptor t-cell therapy.* **J Immunother Cancer 2022;10(11):e005459.**
> [doi:10.1136/jitc-2022-005459](https://doi.org/10.1136/jitc-2022-005459) · PMID 36450377

A combined **hidden Markov model + lasso-penalized logistic regression** that forecasts
the onset and day-by-day trajectory of ICANS in patients receiving CAR-T-cell therapy
(199-patient BWH cohort, leave-one-patient-out cross-validation).

## Reproduce (one command)

Everything runs from the **committed de-identified** data — no download needed.
MATLAB R2019b+ (verified on R2026a), Statistics & Machine Learning Toolbox.

```matlab
run_all      % NTonset_PLR (lasso LOO) -> Forecast_Matricies (HMM forecast metrics)
```

See **[REPRODUCE.md](REPRODUCE.md)** for the figure/table → script → input map, and
**[DATA_SOURCE.md](DATA_SOURCE.md)** for data provenance and de-identification.

## Data

All analysis data are small and committed (`AllBWH_5days_deID.xlsx`, `data_new.mat`,
`LR_prob_allp_Lasso_LOO.mat`, `time_LR.mat`). The patient identifier is a surrogate
`SID`; there are **no MRNs, names, or dates**. The raw MRN-keyed hospital export is
**not** part of this repository (git-ignored) and is never published; de-identification
is lossless (reproduces the raw grade matrix exactly). See [DATA_SOURCE.md](DATA_SOURCE.md).

## Published dataset

Published on BDSP: **https://bdsp.io/content/jyuf7m1cislzg90ew5k6/1.0.0/**
DOIs: [10.60508/4f5c-0z25](https://doi.org/10.60508/4f5c-0z25) (v1.0.0) · [10.60508/5ppk-g283](https://doi.org/10.60508/5ppk-g283) (core).
The de-identified data are also at `s3://bdsp-opendata-credentialed/icans-forecasting-car-t/`.

## License

See [LICENSE.txt](LICENSE.txt).
