library(mgcv) # uniquecombs

single_data_preparation <- function(x, y) {
    # Prepare the data when there is single primary variable y.
    # Returns the missing patterns and number of possible missing patterns for
    # secondary variables x.
    # Args:
    #    x: secondary variables.
    #    y: primary variable.
    # Returns:
    #    data: a list containing the data, missing indicators
    #          and number of possible missing patterns for the
    #          secondary and primary variables.
    n = length(y)
    missing_pattern_x = uniquecombs(!is.na(x))
    num_pattern_x = nrow(missing_pattern_x)
    missing_indicator_x = attr(missing_pattern_x, "index")

    missing_indicator_y = rep(1, n)
    missing_indicator_y[is.na(y)] = 0

    data = list()
    data$y = y
    data$x = as.matrix(x)
    attr(missing_pattern_x, "index") = NULL
    rownames(missing_pattern_x) = 1 : num_pattern_x

    data$missing_pattern_x = missing_pattern_x
    data$R = missing_indicator_x
    data$A = missing_indicator_y
    data$num_pattern_x = num_pattern_x

    return(data)
}

find_available_case_patterns <- function(data, num_observed_variables,
                                        observed_variables) {
    # Find available cases for the auxillary variables x.
    if (num_observed_variables == 0) {
        more_x_observed = 1 : data$num_pattern_x
    } else {
        more_x_observed = which(apply(
            data$missing_pattern_x[, observed_variables, drop=FALSE],
            1, function(v) sum(v) == num_observed_variables))
    }
    return(more_x_observed)
}

single_regression_adjustment <- function(data, fun=identity, binary=FALSE) {
    # Use regression adjustment to estimate the expectation of fun(data$y).
    # By default, linear regression are used when binary is False; otherwise
    # logistic regression are used when binary is True.
    # Args:
    #     data: a list containing data and missing information.
    #     fun:  a function applied to the primary variables.
    #     binary: if fun transforms the primary variables into binary(0/1)
    #             variables, then this argument should be set to TRUE.
    # Returns:
    #     theta_estimate: estimate by regression adjustment.
    n = length(data$y)
    theta_estimate = 0

    # Create a data frame for glm.
    data_x_y = data.frame(data$x, y=fun(data$y))
    colnames(data_x_y) = c(colnames(data$x), "y")

    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(fun(data$y[position_A_1])) / n

    # Logistic regression for binary responses
    family_glm = if (binary) binomial else gaussian
  
    for (i in 1 : data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)
        new_data_position = data$A == 0 & data$R == i

        # If there is no data from pattern i that misses primary variables,
        # then skip this missing pattern as there is no need to do regression
        # adjustment. This might occur in bootstrap.
        if (sum(new_data_position) == 0) {
            next
        }

        if (num_observed_variables > 0) {
            # lm will use available cases to fit the linear regression.
            ra_fm = as.formula(paste("y ~ ", 
                                paste(colnames(data$x)[observed_variables],
                                collapse = '+'), sep=""))
            model_matrix_x = model.matrix.lm(ra_fm, data=data_x_y,
                                             na.action=na.pass)
            ra_fit = glm(ra_fm, data=data_x_y, family=family_glm)
            new_data = data.frame(
                matrix(model_matrix_x[new_data_position, ],
                        nrow=sum(new_data_position)))
            colnames(new_data) = names(ra_fit$coefficients)
            if (!binary) {
                new_data_fit = predict(ra_fit, new_data)
            } else {
                new_data_fit = predict(ra_fit, new_data, type="response")
            }
            theta_est = sum(new_data_fit) / n
        } else {
            # when there are no secondary variables observed.
            ra_fit = mean(fun(data$y[position_A_1]))
            theta_est = ra_fit * sum(new_data_position) / n
        }
        theta_estimate = theta_estimate + theta_est
    }
    return(theta_estimate)
}

single_ipw <- function(data, fun=identity) {
    # Use IPW to estimate the expectation of fun(data$y).
    # By default, logistic regression is used to estimate the odds function.
    # Args:
    #     data: a list containing data and missing information.
    #     fun:  a function applied to the primary variables.
    # Returns:
    #     theta_estimate: estimate by IPW.
    n = length(data$y)
    theta_estimate = 0

    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(fun(data$y[position_A_1])) / n

    for (i in 1:data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)

        # find available cases pattern number for missing pattern i.
        # available cases should have observed_variables observed.
        more_x_observed = find_available_case_patterns(data,
                num_observed_variables, observed_variables)

        # create binary response for estimating the odds function.
        position_0 = data$A == 1 & data$R %in% more_x_observed
        position_1 = data$A == 0 & data$R == i

        # If there is no data from pattern i that misses primary variables,
        # then skip this missing pattern. This might occur in bootstrap.
        if (sum(position_1) == 0) {
            next
        }
        position = position_0 | position_1
        response = rep(0, n)
        response[position_1] = 1

        if (num_observed_variables > 0) {
            x_part = data$x[position, observed_variables]
            ipw_fit = glm(response[position] ~ x_part, family=binomial)
            new_data = cbind(rep(1, sum(position_0)),
                            data$x[position_0, observed_variables])
            odds = exp(new_data %*% ipw_fit$coefficients)
            theta_est = sum(odds * fun(data$y[position_0])) / n
        } else {
            ipw_fit = glm(response[position] ~ 1, family=binomial)
            odds = exp(ipw_fit$coefficients)
            theta_est = sum(odds * fun(data$y[position_0])) / n
        }
        theta_estimate = theta_estimate + theta_est
    }
    return(theta_estimate)
}

single_ipw_sensitivity <- function(data, delta, fun=identity) {
    # Sensitivity analysis with IPW with exponential tilting.
    # Args:
    #     data: a list containing data and missing information.
    #     delta: parameter for sensitivity analysis.
    #     fun: a function applied to the primary variables.
    # Returns:
    #     theta_estimate: estimate by IPW.
    n = length(data$y)
    theta_estimate = 0
    weights = rep(0, n)

    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(fun(data$y[position_A_1]))
    weights[position_A_1] = weights[position_A_1] + 1

    for (i in 1:data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)

        more_x_observed = find_available_case_patterns(data,
                num_observed_variables, observed_variables)
        position_0 = data$A == 1 & data$R %in% more_x_observed
        position_1 = data$A == 0 & data$R == i
        if (sum(position_1) == 0) {
            next
        }
        position = position_0 | position_1
        response = rep(0, n)
        response[position_1] = 1

        # Use exponential tilting for sensitivity analysis.
        odds_part_unobserved = exp(data$y[position_0] * delta)

        if (num_observed_variables > 0) {
            x_part = data$x[position, observed_variables]
            ipw_fit = glm(response[position] ~ x_part, family=binomial)
            new_data = cbind(rep(1, sum(position_0)),
                            data$x[position_0, observed_variables])
            odds = exp(new_data %*% ipw_fit$coefficients)
            odds_sensitivity = odds * odds_part_unobserved
            theta_est_unnorm = sum(odds_sensitivity * fun(data$y[position_0]))
        } else {
            ipw_fit = glm(response[position] ~ 1, family=binomial)
            odds = exp(ipw_fit$coefficients)
            odds_sensitivity = odds * odds_part_unobserved
            theta_est_unnorm = sum(odds_sensitivity * fun(data$y[position_0]))
        }
        weights[position_0] = weights[position_0] + odds_sensitivity
        theta_estimate = theta_estimate + theta_est_unnorm
    }
    # normalize with weight.
    theta_estimate = theta_estimate / sum(weights)
    return(theta_estimate)
}

single_multiply_robust <- function(data, fun=identity, binary=F) {
    # Multiply robust estimation for expectation of fun(data$y).
    # By default, for regression function estimation, linear regression
    # are used when binary is False; otherwise logistic regression are
    # used when binary is True.
    # For odds function, logistic regression are being used.
    # Args:
    #     data: a list containing data and missing information.
    #     fun:  a function applied to the the primary variables.
    #     binary: if fun transform the primary variables into binary, then this
    #             argument should be set to T.
    # Returns:
    #     theta_estimate: estimate by multiply robust estimator.
    n = length(data$y)
    theta_estimate = 0

    # Create a data frame for regression adjustment.
    data_x_y = data.frame(data$x, y=fun(data$y))
    colnames(data_x_y) = c(colnames(data$x), "y")

    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(fun(data$y[position_A_1])) / n

    family_glm = if (binary) binomial else gaussian

    for (i in 1:data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)

        more_x_observed = find_available_case_patterns(data,
                num_observed_variables, observed_variables)
        position_0 = data$A == 1 & data$R %in% more_x_observed
        position_1 = data$A == 0 & data$R == i

        if (sum(position_1) == 0) {
            next
        }

        position = position_0 | position_1
        response = rep(0, n)
        response[position_1] = 1

        if (num_observed_variables > 0) {
            x_part = data$x[position, observed_variables]
            ipw_fit = glm(response[position] ~ x_part, family=binomial)
            new_data_ipw = cbind(rep(1, sum(position_0)),
                            data$x[position_0, observed_variables])
            odds = exp(new_data_ipw %*% ipw_fit$coefficients)

            ra_fm = as.formula(paste("y ~ ",
                                paste(colnames(data$x)[observed_variables],
                                collapse = '+'), sep=""))
            model_matrix_x = model.matrix.lm(ra_fm, data=data_x_y,
                                                na.action=na.pass)
            ra_fit = glm(ra_fm, data=data_x_y, family=family_glm)
            new_data_ra = data.frame(
                matrix(model_matrix_x[position_1, ],
                        nrow=sum(position_1)))
            colnames(new_data_ra) = names(ra_fit$coefficients)
            if (!binary) {
                new_data_fit = predict(ra_fit, new_data_ra)
            } else {
                new_data_fit = predict(ra_fit, new_data_ra, type="response")
            }
            temp = rep(0, n)
            temp[position_0] = (fun(data$y[position_0]) - fitted(ra_fit)) * odds
            temp[position_1] = new_data_fit
            theta_est = mean(temp)
        } else {
            ipw_fit = glm(response[position] ~ 1, family=binomial)
            odds = exp(ipw_fit$coefficients)

            ra_fit = mean(fun(data$y[position_0]))
            temp = rep(0, n)
            temp[position_0] = (fun(data$y[position_0]) - ra_fit) * odds
            temp[position_1] = ra_fit
            theta_est = mean(temp)
        }
        pattern_est[i] = theta_est
        theta_estimate = theta_estimate + theta_est
    }
    return(theta_estimate)
}

single_bootstrap <- function(data, method, n_B=1000, ...) {
    # Bootstrap to obtain the confidence intervals for regression adjustment,
    # IPW and multiply robust estimator.
    # Args:
    #     data: a list containing data and missing information.
    #     method: function for the estimators to perform bootstrap.
    #     n_B: number of bootstrap samples.
    #     ...: optional arguments for method.
    # Returns:
    #     bt_estimate: an array of bootstrap estimates.
    # Example usage:
    #   bt_ra = single_bootstrap(data, single_regression_adjustment,
    #                           fun=indicator, binary=TRUE)
    #   bt_ipw = single_bootstrap(data, single_ipw)
    #   bt_mr = single_bootstrap(data, single_multiply_robust, n_B=500)
    n = length(data$y)
    bt_estimate = rep(0, n_B)

    for (i in 1:n_B) {
        if (i %% 100 == 0) {
            cat("now is the ", i, " iteration\n")
        }
        bt_index = sample(n, size=n, replace=T)
        data_bt = list()
        data_bt$x = data$x[bt_index, ]
        data_bt$y = data$y[bt_index]
        data_bt$A = data$A[bt_index]
        data_bt$R = data$R[bt_index]
        data_bt$missing_pattern_x = data$missing_pattern_x
        data_bt$num_pattern_x = data$num_pattern_x
        bt_estimate[i] = method(data_bt, ...)
    }
    return(bt_estimate)
}