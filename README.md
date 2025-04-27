# Rolling-Forecast-of-Fixed-Income-Portfolio-VaR-and-CVaR
This project builds a dynamic fixed income portfolio risk forecasting system.
I use expanding window rolling forecasting to predict one-day-ahead portfolio risk metrics (VaR and CVaR) over 1000 steps.
Thanks for the insight by github user justeason.

The code has been seperated into different parts(steps).
Methodology
Part A: Data Preparation
1. Load raw Treasury yield curve data (2017–2024).
2. Standardize yields into decimals.
3. Extract maturities from 1M to 30Y.
4. Prepare clean data matrix for modeling.

Part B: Svensson Beta Estimation
1. Fit a Svensson four-factor term structure model to each day's yield curve.
2. Extract daily β₀, β₁, β₂, β₃ coefficients to represent curve dynamics.

Part C: AR(1) Modeling
1. Fit AR(1) time series models individually for each beta series.
2. Extract intercepts, AR(1) coefficients, and residuals.
3. Capture short-term dynamics of yield curve factors.

Part D: Standardizing Residuals
1. Standardize AR(1) residuals to have unit variance.
2. Prepare residuals for DCC-GARCH modeling.

Part E: DCC-GARCH Modeling
1. Fit DCC(1,1)-GARCH(1,1) model to standardized residuals.
2. Model time-varying dynamic correlations between yield curve factors.

Part F: Beta Forecasting
Forecast next day's β₀, β₁, β₂, β₃ using AR(1) models.

Part G: Yield Curve Reconstruction
Reconstruct the forecasted next-day yield curve from predicted betas.

Part H: Portfolio NAV and Log Returns
1. Convert the yield curve into synthetic zero-coupon bond prices.
2. Build an equal-weighted bond portfolio.
3. Compute portfolio NAV and daily log returns.

Part I: Rolling VaR and CVaR Computation
1. Using rolling expanding historical windows:
  - Start with 750 days of history.
  - Expand the training set by 1 day at each step.
  - Forecast one-day-ahead VaR and CVaR.
2. Compute 1%, 2.5%, and 5% VaR and CVaR based on historical returns.

Part J: Expanding Window Rolling Forecast
1. Roll over 1000 steps:
  - Each step extends the training window by 1 day.
  - Forecast one-day-ahead risk metrics.
2. Save rolling VaR/CVaR results into CSV file.
3. Plot rolling risk metrics over steps.

Key Results
1. Successfully tracked rolling forecasts of Portfolio Mean, Portfolio Std, VaR(1%,2.5%,5%), and CVaR(1%,2.5%,5%).
2. Observed smooth dynamics of risk measures over time.
3. Demonstrated the ability of the expanding window model to adapt to evolving market risk.

Folder/File	Description
1. "daily-treasury-rates_2017_2024.csv": Raw input yield curve data
2. "Project.R": Complete R project script
3. "DCC_GARCH_Qt.h5": Saved DCC covariance array (for future Python version)
4. "rolling_CVaR_expanding.csv": Rolling VaR/CVaR results from expanding window forecasting
5. "rolling_CVaR_plots"(Optional): Rolling plots for visualization

Future Work
1. Implement Kupiec and Christoffersen tests to statistically backtest VaR violations.
2. Explore different training window lengths (e.g., 500 days, 1000 days).
3. Extend from 1-step-ahead forecasting to multi-horizon forecasts.
4. Translate project into Python version with a full reproduction.

