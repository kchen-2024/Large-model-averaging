# Implementation of Large Model Averaging (LaMA)

This repository contains the MATLAB implementation of **Large Model Averaging (LaMA)**, along with several classic model selection and model averaging methods for comparison. 

## Data Sources

This analysis uses the `mtcars` and `UScrime` datasets from **R**'s built-in data library (the latter from the `MASS` package).

## Code Structure

- **`plot_empirical_overfit.m`**: Plot the comparison trend between in-sample loss and out-of-sample loss.
- **`plot_risk_double.m`**, **`plot_risk_delete_1.m`**, & **`plot_risk_weight_var.m`**: Plot the limiting behavior of out-of-sample risk under equal, random, and variance-penalized weights, respectively.
- **`plot_sectional_bias_var.m`** & **`plot_sectional_MA_single.m`**: Plot the cross-sections of the three-dimensional graph for risk, bias‑variance decomposition, and model averaging vs. single model.
- **`Mallows_Cp_order.m`**: Prioritize regressors via Mallows' Cₚ‑based forward selection.
- **`test_real_data.m`**: Perform predictions on the sorted dataset.
- **`compare_seven_methods.m`**: Generate samples for the numerical examples and perform model averaging using seven methods:
  1. AIC model selection (AIC)
  2. Smoothed AIC averaging (SAIC)
  3. BIC model selection (BIC)
  4. Smoothed BIC averaging (SBIC)
  5. Mallows model averaging (MMA)
  6. Jackknife model averaging (JMA)
  7. Large model averaging (LaMA)

## Workflow

- **Simulation**: Run `compare_seven_methods.m`.
- **Empirical Analysis**: Modify the data paths as needed, then run `Mallows_Cp_order.m` and `test_real_data.m`.
