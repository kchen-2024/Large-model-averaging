# Implementation of Large Model Averaging (LaMA)

This repository contains the MATLAB implementation of **Large Model Averaging (LaMA)**, along with several classic model selection and model averaging methods for comparison. 

## Data Sources

This analysis uses the `mtcars` and `UScrime` datasets from **R**'s built-in data library (the latter from the `MASS` package).

## Code Structure

 - "**plot_emprical_overfit.m**": Plot the comparison trend between in-sample loss and out-of-sample loss.
 - "**plot_risk_delete_1.m**" &  "**plot_risk_weight_var.m**": Plot the limiting behavior of out-of-sample risk under simple weights and variance-penalized weights, respectively.
 - "**plot_sectional_risk.m**" & "**plot_sectional_bias_var.m**"  & "**plot_sectional_MA_single.m**": Plot the cross-sections of the three-dimensional graph for risk, bias‑variance decomposition, and model averaging vs. single model.
- **`Mallows_Cp_order.m`** – Prioritises regressors via Mallows' Cₚ‑based forward selection.
- **`test_real_data.m`** – Performs predictions on the sorted dataset.
 - "**Large_model_averaging.m**": Core function implementing the LaMA method.
 - "**compare_seven_methods.m**": Generate samples of the numerical examples and perform model averaging using seven methods:
  1. AIC model selection  
  2. Smoothed AIC averaging  
  3. BIC model selection  
  4. Smoothed BIC averaging  
  5. Mallows model averaging  
  6. Jackknife model averaging  
  7. LaMA (Large Model Averaging)


## Workflow

 - Simulation: Run "**compare_seven_methods.m**".
 - Empirical analysis: Modify the path of the data, then run "**Mallows_Cp_order.m**" and "**test_real_data.m**".
