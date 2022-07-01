source("single_primary_variable.R")
source("multiple_primary_variable.R")

diabetes = read.csv("diabetes.csv")
diabetes[diabetes$Glucose == 0, "Glucose"] = NA
diabetes[diabetes$BloodPressure == 0, "BloodPressure"] = NA
diabetes[diabetes$SkinThickness == 0, "SkinThickness"] = NA
diabetes[diabetes$Insulin == 0, "Insulin"] = NA
diabetes[diabetes$BMI == 0, "BMI"] = NA

x = diabetes[, c("Glucose", "BloodPressure", "SkinThickness",
            "BMI")]
y = diabetes[, c("Insulin")]
single_data = single_data_preparation(x, y)


# estimate the mean of insulin.
insulin_mean_ra = single_regression_adjustment(single_data)
insulin_mean_ipw = single_ipw(single_data)
insulin_mean_mr = single_multiply_robust(single_data)

insulin_bt_ra = single_bootstrap(single_data, single_regression_adjustment)
insulin_bt_ipw = single_bootstrap(single_data, single_ipw)
insulin_bt_mr = single_bootstrap(single_data, single_multiply_robust)

insulin_cc = mean(data$Insulin, na.rm=T)
cc_position = single_data$A == 1
insulin_cc_sd = sd(data$Insulin, na.rm=T)

# Sensitivity analysis
insulin_mean_sa = single_ipw_sensitivity(single_data, 0.0001)

# transform insulin to binary and again estimate the proportion
indicator <- function(x) as.numeric(x <= 50)
insulin_ind_mean_ra = single_regression_adjustment(single_data,
                      indicator, binary=T)
insulin_ind_mean_ipw = single_ipw(single_data, indicator)
insulin_ind_mean_mr = single_multiply_robust(single_data, indicator, binary=T)

insulin_ind_bt_mean_ra = single_bootstrap(
    single_data, single_regression_adjustment, fun=indicator, binary=T)
insulin_ind_bt_mean_ipw = single_bootstrap(single_data, single_ipw,
                    fun=indicator)
insulin_ind_bt_mean_mr = single_bootstrap(single_data,
                        single_multiply_robust, fun=indicator, binary=T)

insulin_ind_sa_bt = single_bootstrap(single_data,
                        single_ipw_sensitivity, delta=0.001,
                        fun=indicator)

# Multiple data case
# Estimate the averages of Insulin and SkinThickness
x = diabetes[, c("Glucose", "BloodPressure", "BMI")]
y = diabetes[, c("Insulin", "SkinThickness")]
multiple_data = multiple_data_preparation(x, y)

multiple_average = function(x1, x2) 0.5 * (x1 + x2)
# Estimate the averages of Insulin and SkinThickness.
mt_avg_ra_result = multiple_ra_average(multiple_data)
mt_avg_ipw_result = multiple_ipw(multiple_data, multiple_average)
mt_avg_mr_result = multiple_mr_average(multiple_data)

mt_avg_ra_bt = multiple_bootstrap(multiple_data, multiple_ra_average)
mt_avg_ipw_bt = multiple_bootstrap(multiple_data, multiple_ipw,
                        fun=multiple_average)
mt_avg_mr_bt = multiple_bootstrap(multiple_data, multiple_ra_average)


multiple_indicator = function(x1, x2) {
    as.numeric(x1 <= 50) * as.numeric(x2 <= 50)
}
# Estimate P(Insulin <= 50, SkinThickness <= 50).
mt_ind_ra_result = multiple_ra_indicator(multiple_data, 50)
mt_ind_ipw_result = multiple_ipw(multiple_data, multiple_indicator)
mt_ind_mr_result = multiple_mr_indicator(multiple_data, 50)

mt_ind_ra_bt = multiple_bootstrap(multiple_data, multiple_ra_indicator, a=50)
mt_ind_ipw_bt = multiple_bootstrap(multiple_data, multiple_ipw,
                                   fun=multiple_indicator)
mt_ind_mr_bt = multiple_bootstrap(multiple_data, multiple_mr_indicator, a=50)

mt_ind_sa_bt = multiple_bootstrap(multiple_data, multiple_ipw_sensitivity,
                                  delta=0.0001, fun=multiple_indicator)



# Linear regression of Insulin ~ BMI + BloodPressure
# Treat Glucose and SkinThickness as auxillary variables.
x = diabetes[, c("Glucose", "SkinThickness")]
y = diabetes[, c("Insulin", "BMI", "BloodPressure")]
lm_data = multiple_data_preparation(x, y)
lm_weights = ipw_regression(lm_data)
lm_estimate = lm(Insulin ~ BMI + BloodPressure, data=diabetes,
                weights=lm_weights)$coefficients

# Uncertainty estimate
lm_formula = as.formula("Insulin ~ BMI + BloodPressure")
lm_bt = bootstrap_regression(lm_data, lm, lm_formula)
regression_ci = apply(lm_bt, 2, function(x) quantile(x, probs=c(0.025, 0.975)))

# Complete case analysis for comparison.
response = lm_data$y[lm_data$A == 2, 1]
covariates = as.matrix(lm_data$y[lm_data$A == 2, 2:3])
result_lm = lm(Insulin ~ BMI + BloodPressure, data=diabetes)
cc_lm_ci = confint(result_lm)

# Logistic regression of outcome ~ age + diabetespedigreefunction + bmi +
#                               insulin + glucose + bloodpressure # nolint
# Use SkinThickness as auxillary variables.
x = diabetes[, c("SkinThickness", "Glucose")]
y = diabetes[, c("BMI", "BloodPressure", "Insulin",
                 "Age", "Outcome")]
glm_data = multiple_data_preparation(x, y)
glm_weights = ipw_regression(glm_data)
glm_formula = as.formula(paste("Outcome ~ Age + ",
                        " + BMI + Insulin + BloodPressure", sep=""))

glm_estimate = glm(glm_formula, data=diabetes,
                weights=glm_weights, family=binomial)
glm_bt = bootstrap_regression(glm_data, glm, glm_formula, family=binomial)
regression_ci = apply(lm_bt, 2, function(x) quantile(x, probs=c(0.025, 0.975)))
