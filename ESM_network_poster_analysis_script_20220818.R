## Comparing two ways of modelling contemporaneous effect between ESM variables
## measured in different timeframes

# By Yong Zhang
# Last updated: August 18, 2022

library(tidyverse)
library(lavaan)
library(semPlot)

# Load the dataset
data_analysis <- readRDS("PA_Phys_data_Yong.rds")

# Lag-1 vector autoregressive "VAR(1)" model: modelling contemporaneous effect as
# residual covariance
VAR1_covariance_model <- '
# measurement models
# regression models
Pos ~ phi_11*Pos_l1 +  phi_21*Phys_l1
Phys ~ phi_12*Pos_l1 + phi_22*Phys_l1
# (co)variances latent variables
# (co)variances observed variables
Pos ~~ r*Phys
'

# VAR(1) model: modelling contemporaneous effect as a directed path
VAR1_directed_path_model <- '
# measurement models
# regression models
Pos ~ phi_11*Pos_l1 +  phi_21*Phys_l1 + beta*Phys
Phys ~ phi_12*Pos_l1 + phi_22*Phys_l1
# (co)variances latent variables
# (co)variances observed variables
'

# Fit both models and compare model fit indices
VAR1_covariance_fit <- lavaan(model = VAR1_covariance_model,
                              data = data_analysis,
                              auto.var = TRUE)
summary(VAR1_covariance_fit, fit.measures = TRUE)

VAR1_directed_path_fit <- lavaan(model = VAR1_directed_path_model,
                                 data = data_analysis,
                                 auto.var = TRUE)
summary(VAR1_directed_path_fit, fit.measures = TRUE)

# Path diagrams
semPaths(VAR1_covariance_fit, style = "ram", nCharNodes = 0, what = "model", 
         edge.label.cex = 1.2, edge.label.position = 0.4, 
         label.cex = 1.2, whatLabels = "estimates", sizeLat = 6)
semPaths(VAR1_directed_path_fit, style = "ram", nCharNodes = 0, what = "model", 
         edge.label.cex = 1.2, edge.label.position = 0.4, 
         label.cex = 1.2, whatLabels = "estimates", sizeLat = 6)

# Cross-validation: 7-block ----
### Calculate MSE for all validation pairs for both models
coeff_cov_model <- matrix(nrow = 10, ncol = 7)
mse_cov_model <- list()
coeff_directed_model <- matrix(nrow = 10, ncol = 7)
mse_directed_model <- list()

for (k in 1:7){
  # Based on the residual covariance model:
  # Extract the path coefficients from the model estimated with the training set
  coeff_cov_model[1:10, k] <- parameterEstimates(lavaan(model = VAR1_covariance_model,
                                                        data = data_analysis[-c((16*k-15): (16*k)),],
                                                        auto.var = TRUE))$est
  # Calculate the mean squared error of the fitted values for the test set
  mse_cov_model[k] <- sum((data_analysis[c((16*k-15): (16*k)), 1] - 
                             (coeff_cov_model[1, k]*data_analysis[c((16*k-15): (16*k)), 3] + 
                                coeff_cov_model[2, k]*data_analysis[c((16*k-15): (16*k)), 4]))^2)/16
  # Repeat the same process for the directed path model:
  coeff_directed_model[1:10, k] <-  parameterEstimates(lavaan(model = VAR1_directed_path_model,
                                                              data = data_analysis[-c((16*k-15): (16*k)),],
                                                              auto.var = TRUE))$est
  mse_directed_model[k] <- sum((data_analysis[c((16*k-15): (16*k)), 1] - 
                                  (coeff_directed_model[1, k]*data_analysis[c((16*k-15): (16*k)), 3] + 
                                     coeff_directed_model[2, k]*data_analysis[c((16*k-15): (16*k)), 4] +
                                     coeff_directed_model[3, k]*data_analysis[c((16*k-15): (16*k)), 2]))^2)/16
}

# Compare the MSE between the two models for all training-test combinations
mse_comparison <- tibble(
  ID_test_block = as.factor(1:7),
  MSE_cov_model = as.numeric(mse_cov_model),
  MSE_directed_model = as.numeric(mse_directed_model)
) %>%
  mutate(MSE_diff = MSE_cov_model - MSE_directed_model,
         # Only compare with 0, not super strong
         Better_model = ifelse(MSE_diff < 0, "correlated residuals", "directed-path"))

View(mse_comparison)

# Reshape the tibble to long format for plotting
mse_long <- mse_comparison %>%
  select(ID_test_block, MSE_cov_model, MSE_directed_model) %>%
  rename(TestBlock = ID_test_block, 
         CorrelatedResiduals = MSE_cov_model,
         DirectedPath = MSE_directed_model) %>%
  gather(Model, MSE, c(CorrelatedResiduals, DirectedPath))

ggplot(data = mse_long, aes(x = TestBlock, y = MSE, group = Model)) +
  geom_line(aes(linetype = Model)) +
  geom_point()