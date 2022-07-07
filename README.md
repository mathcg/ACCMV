# Available Complete-Case Missing Value (ACCMV)

Code for paper "Handling Nonmonotone Missing Data with Available Complete-Case Missing Value Assumption". 

- Paper reference: [Handling Nonmonotone Missing Data with Available Complete-Case Missing Value Assumption](https://arxiv.org/abs/2207.02289)

- Contact: mathchenggang@gmail.com

## single_primary_variabe.R
This R script contains functions when there is a single primary variable. 

Functions include:
- `single_data_preparation`: this function prepares the data and finds the missing patterns in the data. 
- `single_regression_adjustment`: this function implements the regression adjustment estimator.
- `single_ipw`: this function implements the inverse probability weighting (IPW) estimator. 
- `single_ipw_sensitivity`: this function implements the sensitivity analysis with exponential tilting for IPW estimator. 
- `single_multiply_robust`: this function implements the multiply robust estimator.
- `single_bootstrap`: this function implements the bootstrap to estimate the confidence intervals for the above estimators. 

## multiple_primary_variables.R
This R script contains functions when there are multiple primary variable. If not specified, the function below assumes that there are two primary variables. 

Functions include:
- `multiple_data_preparation`: similar to the `single_data_preparation` function above. 
- `multiple_ra_average`: this function implements the regression adjustment estimator when the outcomes are continuous.
- `multiple_ra_indicator`: this function implements the regression adjustment estimator when the outcomes are binary.
- `multiple_ipw`: this function implements the inverse probability weighting (IPW) estimator. 
- `multiple_ipw_sensitivity`: this function implements the sensitivity analysis with exponential tilting for IPW estimator. 
- `multiple_mr_average`: this function implements the multiply robust estimator when the outcomes are continuous. 
- `multiple_mr_indicator`: this function implements the multiply robust estimator when the outcomes are binary. 
- `multiple_bootstrap`: this function implements the bootstrap to estimate the confidence intervals for the above estimators. 
- `ipw_regression`: this function output weights for IPW and the weights can be used for estimating the regression parameters. 
- `bootstrap_regression`: this function implements bootstrap for estimating regression parameters. 

## diabetes_demo.R
An R script for demonstrating how to use the functions mentioned above. We apply our estimators to the Pima Indians Diabetes Dataset. 

## diabetes.csv
The Pima Indians Diabetes Dataset. It is downloaded from https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database. This dataset is to diagnostically predict whether or not a patient has diabetes, based on certain diagnostic measurements included in the dataset. 


