source("single_primary_variable.R")
source("multiple_primary_variable.R")

diabetes = read.csv("diabetes.csv")
# Glucose, BloodPressure, SkinThickness, Insulin and BMI cannot have values
# being 0. This suggests that these values are in fact missing. For this reason,
# we replace their values by NA. 
diabetes[diabetes$Glucose == 0, "Glucose"] = NA
diabetes[diabetes$BloodPressure == 0, "BloodPressure"] = NA
diabetes[diabetes$SkinThickness == 0, "SkinThickness"] = NA
diabetes[diabetes$Insulin == 0, "Insulin"] = NA
diabetes[diabetes$BMI == 0, "BMI"] = NA

# To estimate the mean of Insulin E[Insulin], we treat Glucose, BloodPressure,
# SkinThickness,and BMI as auxillary variables and treat Insulin as
# primary variables.
x = diabetes[, c("Glucose", "BloodPressure", "SkinThickness",
            "BMI")]
y = diabetes[, c("Insulin")]

# single_data contains data and missing information, like missing patterns
# for both primary and auxillary variables.
single_data = single_data_preparation(x, y)

# Estimate with regression adjustment, ipw and multiply robust estimator.
insulin_mean_ra = single_regression_adjustment(single_data)
insulin_mean_ipw = single_ipw(single_data)
insulin_mean_mr = single_multiply_robust(single_data)

# Bootstrap to obtain 95% confidence intervals. 
insulin_bt_ra = single_bootstrap(single_data, single_regression_adjustment)
insulin_ra_ci = quantile(insulin_bt_ra, probs=c(0.025, 0.975))

insulin_bt_ipw = single_bootstrap(single_data, single_ipw)
insulin_ipw_ci = quantile(insulin_bt_ipw, probs=c(0.025, 0.975))

insulin_bt_mr = single_bootstrap(single_data, single_multiply_robust)
insulin_mr_ci = quantile(insulin_bt_mr, probs=c(0.025, 0.975))

# Sensitivity analysis for the estimate with delta = 0.0001
insulin_mean_sa = single_ipw_sensitivity(single_data, 0.0001)

# Transform Insulin to binary and estimate the proportion P(Insulin <= 50).
indicator <- function(x) as.numeric(x <= 50)
insulin_ind_mean_ra = single_regression_adjustment(single_data,
                      indicator, binary=T)
insulin_ind_mean_ipw = single_ipw(single_data, indicator)
insulin_ind_mean_mr = single_multiply_robust(single_data, indicator, binary=T)

# Bootstrap to estimate the confidence intervals
insulin_ind_bt_mean_ra = single_bootstrap(
    single_data, single_regression_adjustment, fun=indicator, binary=T)
insulin_ind_ra_ci = quantile(insulin_ind_bt_mean_ra, probs=c(0.025, 0.975))

insulin_ind_bt_mean_ipw = single_bootstrap(single_data, single_ipw,
                    fun=indicator)
insulin_ind_ipw_ci = quantile(insulin_ind_bt_mean_ipw, probs=c(0.025, 0.975))

insulin_ind_bt_mean_mr = single_bootstrap(single_data,
                        single_multiply_robust, fun=indicator, binary=T)
insulin_ind_mr_ci = quantile(insulin_ind_bt_mean_mr, probs=c(0.025, 0.975))

# bootstrap for sensitivity analysis.
insulin_ind_sa_bt = single_bootstrap(single_data,
                        single_ipw_sensitivity, delta=0.001,
                        fun=indicator)
insulin_ind_sa_ci = quantile(insulin_ind_sa_bt, probs=c(0.025, 0.975))

# Multiple data case
# Estimate the averages of Insulin and SkinThickness
# i.e., E[0.5 * (Insulin + SkinThickness)]
x = diabetes[, c("Glucose", "BloodPressure", "BMI")]
y = diabetes[, c("Insulin", "SkinThickness")]
multiple_data = multiple_data_preparation(x, y)

multiple_average = function(x1, x2) 0.5 * (x1 + x2)
mt_avg_ra_result = multiple_ra_average(multiple_data)
mt_avg_ipw_result = multiple_ipw(multiple_data, multiple_average)
mt_avg_mr_result = multiple_mr_average(multiple_data)

# Bootstrap to get the 95% confidence intervals.
mt_avg_ra_bt = multiple_bootstrap(multiple_data, multiple_ra_average)
mt_avg_ra_ci = quantile(mt_avg_ra_bt, probs=c(0.025, 0.975))

mt_avg_ipw_bt = multiple_bootstrap(multiple_data, multiple_ipw,
                        fun=multiple_average)
mt_avg_ipw_ci = quantile(mt_avg_ipw_bt, probs=c(0.025, 0.975))

mt_avg_mr_bt = multiple_bootstrap(multiple_data, multiple_ra_average)
mt_avg_mr_ci = quantile(mt_avg_mr_bt, probs=c(0.025, 0.975))

multiple_indicator = function(x1, x2) {
    as.numeric(x1 <= 50) * as.numeric(x2 <= 50)
}
# Estimate P(Insulin <= 50, SkinThickness <= 50).
mt_ind_ra_result = multiple_ra_indicator(multiple_data, 50)
mt_ind_ipw_result = multiple_ipw(multiple_data, multiple_indicator)
mt_ind_mr_result = multiple_mr_indicator(multiple_data, 50)

mt_ind_ra_bt = multiple_bootstrap(multiple_data, multiple_ra_indicator, a=50)
mt_ind_ra_ci = quantile(mt_ind_ra_bt, probs=c(0.025, 0.975))

mt_ind_ipw_bt = multiple_bootstrap(multiple_data, multiple_ipw,
                                   fun=multiple_indicator)
mt_ind_ipw_ci = quantile(mt_ind_ipw_bt, probs=c(0.025, 0.975))

mt_ind_mr_bt = multiple_bootstrap(multiple_data, multiple_mr_indicator, a=50)
mt_ind_mr_ci = quantile(mt_ind_mr_bt, probs=c(0.025, 0.975))

# Sensitivity analysis
mt_ind_sa_bt = multiple_bootstrap(multiple_data, multiple_ipw_sensitivity,
                                  delta=0.0001, fun=multiple_indicator)
mt_ind_sa_ci = quantile(mt_ind_sa_bt, probs=c(0.025, 0.975))


# Linear regression of Insulin ~ BMI + BloodPressure
# Treat Glucose and SkinThickness as auxillary variables.
x = diabetes[, c("Glucose", "SkinThickness")]
y = diabetes[, c("Insulin", "BMI", "BloodPressure")]
lm_data = multiple_data_preparation(x, y)
# Use weights to estimate the regression parameters
lm_weights = ipw_regression(lm_data)
lm_formula = as.formula("Insulin ~ BMI + BloodPressure")
lm_estimate = lm(lm_formula, data=diabetes,
                weights=lm_weights)$coefficients

# Use bootstrap to obtain the confidence intervals.
lm_bt = bootstrap_regression(lm_data, lm, lm_formula)
lm_est_ci = apply(lm_bt, 2, function(x) quantile(x, probs=c(0.025, 0.975)))

# Logistic regression of outcome ~ age + bmi +
#                               insulin + bloodpressure # nolint
# Use SkinThickness and glucose as auxillary variables.
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
glm_ci = apply(glm_bt, 2, function(x) quantile(x, probs=c(0.025, 0.975)))
colnames(glm_ci) = c("Intercept", attr(terms(glm_formula), "term.labels"))
