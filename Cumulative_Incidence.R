# Load required packages
library(tibble)
library(dplyr)
library(kableExtra)
library(cmprsk)
library(ggplot2)
library(patchwork)

# Number of observations
n <- 4000

# Modify the creation of the dummy dataset to include 'age' as a significant factor
set.seed(1214223)

# Create a dummy dataset with group-specific event probabilities
dummy_data <- tibble(
  time_to_event = rexp(n, rate = 0.1),  # Generate random survival times
  region = sample(c("Outside London", "Inside London"), size = n, replace = TRUE),
  age = rnorm(n, mean = 8, sd = 10),
  # Adjust the probability to make both age and region significantly affect the event status
  status = ifelse(region == "Inside London",
                  rbinom(n, size = 1, prob = plogis(-1.5 - 0.02 * age)),  # Stronger logistic effect of age for "Inside London"
                  rbinom(n, size = 1, prob = plogis(-0.5 - 0.02 * age)))  # Weaker logistic effect of age for "Outside London"
)

# Create the Fine-Gray model
fg_model <- crr(
  ftime = dummy_data$time_to_event,    # Time to event
  fstatus = dummy_data$status,         # Status (event vs censoring)
  cov1 = model.matrix(~ region + age, data = dummy_data)[, -1]  # Covariates (excluding intercept column)
)

# Summary of the model
summary(fg_model)

# Create cumulative incidence data using 'cuminc' function for uncorrected cumulative incidence
fit <- cuminc(ftime = dummy_data$time_to_event, fstatus = dummy_data$status, group = dummy_data$region)

# Extract and organize uncorrected data into a tidy data frame for ggplot, including cumulative incidence estimates
uncorrected_ci_data <- do.call(rbind, lapply(names(fit), function(region_name) {
  if ("est" %in% names(fit[[region_name]])) {
    data.frame(
      time = fit[[region_name]]$time,
      cuminc = fit[[region_name]]$est,
      region = gsub(" 1", "", region_name)  # Clean the region name
    )
  }
}))

# Extract coefficients from the Fine-Gray model
coefficients <- summary(fg_model)$coef

# Extract coefficients for "regionOutside London" and "age"
outside_london_coef <- coefficients["regionOutside London", "coef"]
age_coef <- coefficients["age", "coef"]

# Calculate the average age in the dataset
mean_age <- mean(dummy_data$age)

# Extract the coefficients vector from the Fine-Gray model
coef_vector <- coefficients[, "coef"]

# Create the design matrix for the covariates (excluding intercept)
covariate_matrix <- model.matrix(~ region + age, data = dummy_data)[, -1]  # Adjust as needed for additional covariates

# Calculate the linear predictor for each observation
# This is the dot product of the covariate values and the corresponding coefficients
linear_predictor <- covariate_matrix %*% coef_vector

# Use the uncorrected cumulative incidence as a baseline and adjust for the effect of all covariates
corrected_ci_data <- uncorrected_ci_data %>%
  mutate(
    # Apply the linear predictor to calculate the corrected cumulative incidence
    corrected_cuminc = cuminc * exp(linear_predictor[match(region, dummy_data$region)]),
    # Calculate standard errors for confidence intervals
    se_cuminc = sqrt(diag(fg_model$var)[1])  # Simplified extraction of standard error for illustration purposes
  ) %>%
  mutate(
    ci_upper = corrected_cuminc + 1.96 * se_cuminc,
    ci_lower = corrected_cuminc - 1.96 * se_cuminc,
    type = "Corrected"
  )

# Add a new column to indicate whether the data is corrected or uncorrected
uncorrected_ci_data <- uncorrected_ci_data %>%
  mutate(type = "Uncorrected")

corrected_ci_data <- corrected_ci_data %>%
  mutate(type = "Corrected")

# Combine both corrected and uncorrected data
combined_ci_data <- bind_rows(
  uncorrected_ci_data %>% rename(cif = cuminc),
  corrected_ci_data %>% rename(cif = corrected_cuminc)
)


# Plot both the corrected (Fine-Gray) and uncorrected cumulative incidence curves for both groups
ggplot(combined_ci_data, aes(x = time, y = cif, color = region, linetype = type)) +
  geom_line(size = 1.2) +
  # Labels and styling
  labs(
    title = "Cumulative Incidence Function: Corrected (Fine-Gray) vs. Uncorrected",
    x = "Time to Event",
    y = "Cumulative Incidence",
    color = "Region",
    linetype = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "bottom"
  )


# Plot corrected cumulative incidence curves with confidence intervals
cuminc_plot_corrected <- ggplot(corrected_ci_data, aes(x = time, y = corrected_cuminc, color = region)) +
  geom_ribbon(aes(x = time, ymin = ci_lower, ymax = ci_upper, fill = region), alpha = 0.2) +
  geom_line(size = 1.2) +
  # Labels and styling
  labs(
    title = "Cumulative Incidence Function: Corrected (Fine-Gray)",
    x = "Time to Event",
    y = "Cumulative Incidence"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "bottom"
  )

# Calculate the number at risk for the unadjusted risk table
time_intervals <- seq(0, max(dummy_data$time_to_event), by = 10)  # Set time intervals (e.g., every 10 units)

risk_data <- lapply(time_intervals, function(t) {
  dummy_data %>%
    filter(time_to_event >= t) %>%
    group_by(region) %>%
    summarise(at_risk = n(), .groups = 'drop') %>%
    mutate(time = t)
}) %>%
  bind_rows()

# Plot the risk table with unadjusted number at risk
risk_table_plot <- ggplot(risk_data, aes(x = time, y = region, label = at_risk)) +
  geom_text(size = 3) +
  labs(
    x = "Time to Event",
    y = "Region",
    title = "Number at Risk"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove grid lines for clarity
    panel.grid.minor = element_blank()
  )

# Combine the corrected cumulative incidence plot and unadjusted risk table using patchwork
final_plot <- cuminc_plot_corrected / risk_table_plot + plot_layout(heights = c(3, 1))

# Print the combined plot
final_plot


#Revised Approach:
#  To maintain statistical validity, the following revised approach is suggested:
#
#  Keep the Number at Risk as an Empirical Count:
#  Retain the original, unadjusted counts of participants at risk over time.
#Use the Fine-Gray Model for Cumulative Incidence:
#  Use the Fine-Gray model coefficients to derive adjusted cumulative incidence curves, showing the covariate effects on survival or cumulative incidence.
#Visualize Covariate Effects through Model Outputs:
#  Present the effects of covariates by plotting adjusted survival or cumulative incidence curves.
#Here's the revised version of the code, where we maintain the empirical number at risk while showcasing the covariate effects through adjusted cumulative incidence:

