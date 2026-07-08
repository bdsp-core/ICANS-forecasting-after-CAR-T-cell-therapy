% Full ICANS-forecasting reproduction from the committed de-identified source.
close all; clear; clc;
fprintf('--- Step 1: NTonset_PLR (lasso LOO -> LR_prob) ---\n'); NTonset_PLR;
fprintf('--- Step 2: Forecast_Matricies (HMM forecast metrics) ---\n');
run('Forecast_Matricies.m');
fprintf('REPRODUCED forecast error by horizon (mean wMAPE, h=1..7):\n');
disp(mean(wMAPE,2)');
fprintf('mean |bias| by horizon:\n'); disp(mean(abs(Baias),2)');
fprintf('ALL_STEPS_OK\n');
