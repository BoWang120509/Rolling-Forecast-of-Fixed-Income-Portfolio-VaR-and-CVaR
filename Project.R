# Load necessary libraries
library(rmgarch)
library(fGarch)
library(rhdf5)  

# Read the expanded Treasury yield curve data
# data <- read.csv('C:/Users/19493/Desktop/MSMFT第二学期/MF728/main_project/daily-treasury-rates_2017_2024.csv')
data <- read.csv('C:/Users/19493/Desktop/MSMFT第二学期/MF728/main_project/daily-treasury-rates_2017_2024.csv')

# Check the data structure
head(data)


#######Step 1#######
# Convert Date to Date type
# data$Date <- as.Date(data$Date)
date_vec <- as.Date(data$Date, format="%m/%d/%Y")

# Extract maturities columns
maturity_cols <- c('X1.Mo', 'X3.Mo', 'X6.Mo', 'X1.Yr', 'X2.Yr', 'X3.Yr', 'X5.Yr', 'X7.Yr', 'X10.Yr', 'X20.Yr', 'X30.Yr')
# yield_data <- data[, maturity_cols]
yield_data <- data[, maturity_cols]

if (max(yield_data, na.rm=TRUE) > 1) {
  yield_data <- yield_data / 100
}

data <- yield_data

# Define months for each maturity
months <- c(1, 3, 6, 12, 24, 36, 60, 84, 120, 240, 360)
maturities <- months / 12

colnames(data) <- months


#######Step 2#######
# Define Svensson Yield Function.
# Fir a Svensson model (four-factor term structure model) to each day's yield curve. Extract beta coefficients for each day.
# Reduce the entire yield curve into a few key dynamic factors that can be modeled over time.
svensson_yield <- function(tau, beta0, beta1, beta2, beta3, lambda1, lambda2) {
  term1 <- (1 - exp(-lambda1 * tau)) / (lambda1 * tau)
  term2 <- term1 - exp(-lambda1 * tau)
  term3 <- (1 - exp(-lambda2 * tau)) / (lambda2 * tau) - exp(-lambda2 * tau)
  return(beta0 + beta1 * term1 + beta2 * term2 + beta3 * term3)
}

# Function to fit Svensson for one day's yield curve
fit_svensson <- function(yields, maturities, lambda1 = 0.5, lambda2 = 0.1) {
  
  # Internal function: model to optimize
  svensson_model <- function(params, tau, yields) {
    beta0 <- params[1]
    beta1 <- params[2]
    beta2 <- params[3]
    beta3 <- params[4]
    fitted <- svensson_yield(tau, beta0, beta1, beta2, beta3, lambda1, lambda2)
    return(sum((fitted - yields)^2))
  }
  
  start_params <- c(0.03, -0.02, 0.02, 0.01) # This is the defalt by me.
  
  fit <- optim(
    par = start_params,
    fn = svensson_model,
    tau = maturities,
    yields = yields,
    method = "BFGS"
  )
  
  return(fit$par)
}

# Now apply to every day
beta.df <- as.data.frame(matrix(NA, nrow=nrow(data), ncol=4))
colnames(beta.df) <- c('beta0', 'beta1', 'beta2', 'beta3')

for (i in 1:nrow(data)) {
  tryCatch({
    beta.df[i, ] <- fit_svensson(as.numeric(data[i, ]), maturities)
  }, error = function(e) {
    beta.df[i, ] <- rep(NA, 4)
  })
}

cat("\n Part 2: Svensson Beta estimates preview \n")
print(head(beta.df))
flush.console()




#######Step 3#######
# Fit an AR(1) process to each beta series individually. Extract intercepts, AR(1) coefficients, and residuals.
# Capture the time series dynamics of the yield curve factors.
fit_ar1 <- function(series) {
  x <- series[1:(length(series) - 1)]
  y <- series[2:length(series)]
  
  ar_model <- lm(y ~ x)
  
  intercept <- ar_model$coefficients[1]
  phi <- ar_model$coefficients[2]
  residuals <- ar_model$residuals
  
  return(list(intercept = intercept, phi = phi, residuals = residuals))
}


# Now apply AR(1) fitting for each beta
ar.resid.df <- as.data.frame(matrix(NA, nrow=nrow(beta.df)-1, ncol=4))
colnames(ar.resid.df) <- colnames(beta.df)

ar.coef.df <- as.data.frame(matrix(NA, nrow=2, ncol=4))
rownames(ar.coef.df) <- c('intercept', 'phi')
colnames(ar.coef.df) <- colnames(beta.df)

for (i in 1:ncol(beta.df)) {
  ar_fit <- fit_ar1(beta.df[,i])
  
  ar.resid.df[,i] <- ar_fit$residuals
  ar.coef.df['intercept',i] <- ar_fit$intercept
  ar.coef.df['phi',i] <- ar_fit$phi
}

cat("\n Part 3: AR(1) Residuals preview \n")
print(head(ar.resid.df))
flush.console()

cat("\n Part 3: AR(1) Coefficients \n")
print(ar.coef.df)
flush.console()




#######Step 4#######
# Function to standardize residuals (make variance ~1)
# Prepare the residuals for DCC-GARCH modeling, ensuring consistent scale across factors.
standardize_residuals <- function(residuals_df) {
  
  standardized_resid <- as.data.frame(matrix(NA, nrow=nrow(residuals_df), ncol=ncol(residuals_df)))
  colnames(standardized_resid) <- colnames(residuals_df)
  
  for (i in 1:ncol(residuals_df)) {
    sd_resid <- sd(residuals_df[,i], na.rm=TRUE)
    standardized_resid[,i] <- residuals_df[,i] / sd_resid
  }
  
  return(standardized_resid)
}

# Apply standardization
std.resid.df <- standardize_residuals(ar.resid.df)

cat("\nPart 4: Standardized residuals preview\n")
print(head(std.resid.df))
flush.console()




#######Step 5#######
# Define univariate GARCH(1,1) model for each beta.
# Estimate time-varying conditional correlations among yield curve factors.
# Capture the dynamic covariance structure between yield curve factors over time.
garch11.spec <- ugarchspec(
  mean.model = list(armaOrder = c(0, 0)),
  variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
  distribution.model = "norm"
)

# Define DCC(1,1) specification
dcc.spec <- dccspec(
  uspec = multispec(replicate(4, garch11.spec)),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)

# Fit DCC-GARCH model
dcc.fit <- dccfit(dcc.spec, data = std.resid.df)

cat("\n--- Part E: DCC-GARCH fitted model summary ---\n")
print(dcc.fit)
flush.console()

# Extract last Q covariance matrix
last_Q_matrix <- rcov(dcc.fit)[,,ncol(rcov(dcc.fit))]

cat("\n--- Part E: Last Conditional Covariance Matrix (Qt) ---\n")
print(last_Q_matrix)
flush.console()

# From this part, I try to create a h5 file for python version.
# Just load the file and extract it.
# Skip the DCC part for python (but MUST write mannually DCC in python) and run the rest.
Qt_array <- rcov(dcc.fit)
cat("\n--- Shape of Qt array ---\n")
print(dim(Qt_array))

output_path <- "C:/Users/19493/Desktop/MSMFT第二学期/MF728/main_project/DCC_GARCH_Qt.h5"

h5createFile(output_path)
h5createDataset(output_path, "Qt_array", dims=dim(Qt_array), storage.mode="double")
h5write(Qt_array, output_path, "Qt_array")

cat("\n Qt array saved successfully to: \n", output_path, "\n")





#######Step 6#######
# Function to predict next beta using AR(1) coefficients.
# Generate predicted factor values to reconstruct the forecasted yield curve.
predict_next_beta <- function(last_beta, ar_coef_df) {
  
  beta_pred <- rep(NA, 4)
  names(beta_pred) <- colnames(ar_coef_df)
  
  for (i in 1:length(last_beta)) {
    intercept <- ar_coef_df['intercept', i]
    phi <- ar_coef_df['phi', i]
    beta_pred[i] <- intercept + phi * last_beta[i]
  }
  
  return(beta_pred)
}

# Last observed beta
last_beta <- as.numeric(beta.df[nrow(beta.df), ])

# Predict next beta
beta_pred <- predict_next_beta(last_beta, ar.coef.df)

cat("\n Part 6: Predicted next beta \n")
print(beta_pred)
flush.console()





#######Step 7#######
# Function to generate yield curve using predicted betas.
# Rebuild the next day's yield curve based on predicted beta coefficients using the Svensson model.
# Obtain a forecasted yield curve for portfolio construction and risk analysis.
generate_predicted_yield_curve <- function(beta_pred, maturities, lambda1 = 0.5, lambda2 = 0.1) {
  term1 <- (1 - exp(-lambda1 * maturities)) / (lambda1 * maturities)
  term2 <- term1 - exp(-lambda1 * maturities)
  term3 <- (1 - exp(-lambda2 * maturities)) / (lambda2 * maturities) - exp(-lambda2 * maturities)
  
  yields_pred <- beta_pred[1] + beta_pred[2] * term1 + beta_pred[3] * term2 + beta_pred[4] * term3
  
  return(yields_pred)
}

# Generate predicted yields
mu.pred <- generate_predicted_yield_curve(beta_pred, maturities)

cat("\n Part 7: Predicted Yield Curve \n")
print(mu.pred)
flush.console()






#######Step 8#######
# Convert the predicted yield curve into synthetic bond prices (zero-coupon approximation).
# Construct an equal-weighted bond portfolio. Calculate the portfolio's NAV and daily log returns.
# Translate yield curve forecasts into actual portfolio value changes to assess potential losses. 
num_days <- nrow(beta.df)
num_terms <- length(maturities)

# Yield curve for each period
yield_matrix <- matrix(NA, nrow=num_days, ncol=num_terms)
for (i in 1:num_days) {
  yield_matrix[i, ] <- svensson_yield(
    maturities,
    beta.df[i, 1], beta.df[i, 2], beta.df[i, 3], beta.df[i, 4],
    lambda1 = 0.5, lambda2 = 0.1
  )
}
colnames(yield_matrix) <- paste0(maturities, "Y")

# Using the yield from the previous code to generate the equal weighted bond basket. That's right, I changed to equal weighted.
price_matrix <- matrix(NA, nrow=num_days, ncol=num_terms)
for (j in 1:num_terms) {
  price_matrix[,j] <- exp(-yield_matrix[,j] * maturities[j])
}

# The NAV of the portfolio
weights <- rep(1/num_terms, num_terms)
comb_nav <- as.numeric(price_matrix %*% weights) 
comb_ret <- diff(log(comb_nav))

cat("\n Part 8: Equal-Weight Multi-Period Portfolio's NAV and Daily Yield \n")
print(head(comb_ret))
cat("Mean：", mean(comb_ret, na.rm=TRUE), "\n")
cat("STD：", sd(comb_ret, na.rm=TRUE), "\n")
flush.console()



#######Step 9#######
# Based on rolling historical returns, compute 1%, 2.5%, and 5% VaR and CVaR.
# Quantify portfolio downside risk under different confidence levels.
alphas <- c(0.01, 0.025, 0.05)
VaR <- sapply(alphas, function(a) quantile(comb_ret, a, na.rm=TRUE))
CVaR <- sapply(1:length(alphas), function(i) mean(comb_ret[comb_ret <= VaR[i]], na.rm=TRUE))

names(VaR) <- paste0(alphas*100, "% VaR")
names(CVaR) <- paste0(alphas*100, "% CVaR")

cat("\n Part 9: Equal-Weight Multi-Period Portfolio's VaR/CVaR \n")
print(VaR)
print(CVaR)
flush.console()



##############################Everything above can be consider as functions. The step 10 is the loop####################################


#######Step 10####### 
# Implement an expanding window rolling forecast:
# Start with 750 days of historical data.
# Each day, extend the training window by one additional day.
# Predict the next day's VaR and CVaR, repeating for 500–1001 steps.
# Dynamically track how the risk profile evolves over time as more historical information is accumulated.
results <- data.frame(
  PortfolioMean = rep(NA, 1000),
  PortfolioStd = rep(NA, 1000),
  VaR_1pct = rep(NA, 1000),
  CVaR_1pct = rep(NA, 1000),
  VaR_2_5pct = rep(NA, 1000),
  CVaR_2_5pct = rep(NA, 1000),
  VaR_5pct = rep(NA, 1000),
  CVaR_5pct = rep(NA, 1000)
)

# The initial day for starting
train_days_start <- 750

# This is the loop starts
for (ind in 1:1000) {
  
  train_range <- 1:(train_days_start + ind)
  
  beta_window <- beta.df[train_range, ]
  
  # Yield Matrix
  yield_matrix <- matrix(NA, nrow=nrow(beta_window), ncol=length(maturities))
  for (i in 1:nrow(beta_window)) {
    yield_matrix[i, ] <- svensson_yield(
      maturities, 
      beta_window[i,1], beta_window[i,2], beta_window[i,3], beta_window[i,4],
      lambda1 = 0.5, lambda2 = 0.1
    )
  }
  
  # Bond price Matrix
  price_matrix <- matrix(NA, nrow=nrow(yield_matrix), ncol=length(maturities))
  for (j in 1:length(maturities)) {
    price_matrix[,j] <- exp(-yield_matrix[,j] * maturities[j])
  }
  
  # Port NAV
  weights <- rep(1/length(maturities), length(maturities))
  comb_nav <- as.numeric(price_matrix %*% weights)
  
  # Port log return
  comb_ret <- diff(log(comb_nav))
  
  # VaR/CVaR
  alphas <- c(0.01, 0.025, 0.05)
  VaRs <- sapply(alphas, function(a) quantile(comb_ret, a, na.rm=TRUE))
  CVaRs <- sapply(1:length(alphas), function(i) mean(comb_ret[comb_ret <= VaRs[i]], na.rm=TRUE))
  
  results$Step[ind] <- ind
  results$PortfolioMean[ind] <- mean(comb_ret, na.rm=TRUE)
  results$PortfolioStd[ind] <- sd(comb_ret, na.rm=TRUE)
  results$VaR_1pct[ind] <- VaRs[1]
  results$CVaR_1pct[ind] <- CVaRs[1]
  results$VaR_2_5pct[ind] <- VaRs[2]
  results$CVaR_2_5pct[ind] <- CVaRs[2]
  results$VaR_5pct[ind] <- VaRs[3]
  results$CVaR_5pct[ind] <- CVaRs[3]
  

  if (ind %% 50 == 0) {
    cat("Completed the", ind, "steps\n")
    flush.console()
  }
}

cat("\n Part 10: Expanding Window Finished \n")
print(head(results))
flush.console()

write.csv(results, "C:/Users/19493/Desktop/MSMFT第二学期/MF728/main_project/rolling_CVaR.csv", row.names = FALSE)
cat("\n The Rolling Completed \n")





#######Step 11####### 
# Plotting the curves
results <- read.csv('C:/Users/19493/Desktop/MSMFT第二学期/MF728/main_project/rolling_CVaR_expanding.csv')
steps <- 1:nrow(results)

columns_to_plot <- c("PortfolioMean", "PortfolioStd", 
                     "VaR_1pct", "CVaR_1pct",
                     "VaR_2_5pct", "CVaR_2_5pct",
                     "VaR_5pct", "CVaR_5pct")

neg_columns <- c("VaR_1pct", "CVaR_1pct", "VaR_2_5pct", "CVaR_2_5pct", "VaR_5pct", "CVaR_5pct")

par(mfrow=c(2,4))

for (colname in columns_to_plot) {
  
  if (colname %in% neg_columns) {
    y_values <- -results[[colname]]
  } else {
    y_values <- results[[colname]]
  }
  
  plot(steps, y_values,
       type = "l",
       col = "blue",
       lwd = 2,
       xlab = "Step",
       ylab = colname,
       main = paste("Rolling", colname, "over Steps"),
       xaxt = "n") 
  axis(1, at=seq(0, max(steps), by=50))
  grid()
}

par(mfrow=c(1,1))






