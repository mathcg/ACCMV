library(mgcv)

multiple_data_preparation <- function(x, y) {
    # Prepare the data when there are multiple primary variables.
    # Returns the missing patterns and number of possible missing patterns for
    # secondary variables x.
    # Args:
    #    x: secondary variables.
    #    y: primary variables.
    # Returns:
    #    data: a list containing the data, missing indicators
    #          and number of possible missing patterns for the
    #          secondary and primary variables.
    missing_pattern_x = uniquecombs(!is.na(x))
    num_pattern_x = nrow(missing_pattern_x)
    missing_indicator_x = attr(missing_pattern_x, "index")

    missing_pattern_y = uniquecombs(!is.na(y))
    num_pattern_y = nrow(missing_pattern_y)
    missing_indicator_y = attr(missing_pattern_y, "index")

    data = list()
    data$y = as.matrix(y)
    data$x = as.matrix(x)
    attr(missing_pattern_x, "index") = NULL
    rownames(missing_pattern_x) = 1 : num_pattern_x
    data$missing_pattern_x = missing_pattern_x

    attr(missing_pattern_y, "index") = NULL
    rownames(missing_pattern_y) = 1 : num_pattern_y
    data$missing_pattern_y = missing_pattern_y

    data$R = missing_indicator_x
    data$A = missing_indicator_y
    data$num_pattern_x = num_pattern_x
    data$num_pattern_y = num_pattern_y
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

multiple_ra_average <- function(data) {
    # Use regression adjustment to estimate the averages of
    # primary variables.
    # We assume that there are two primary variables and this
    # function can be easily modified to extend to more primary
    # variables.
    # Args:
    #   data: a list containing data and missing information.
    # Returns:
    #   theta_estimate: estimate of the averages.
    n = length(data$A)
    theta_estimate = 0

    # pattern number for complete case
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    theta_estimate = theta_estimate + 
        sum(0.5 * (data$y[position_A_1, 1] + data$y[position_A_1, 2])) / n

    for (i in 1 : data$num_pattern_y) {
        # Skip this missing pattern for y if y is fully observed.
        if (i == cp_num) {
            next
        }

        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        missed_variables_y = which(missing_pattern_y == F)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            # Find out all available cases for the auxillary variables x.
            more_x_observed = find_available_case_patterns(data,
                num_observed_variables_x, observed_variables_x)

            # Complete cases for primary variables and available cases
            # for auxillary variable.
            position = data$A == cp_num & data$R %in% more_x_observed

            new_data_position = data$A == i & data$R == j

            # There are no data that has missing pattern i for y and missing
            # pattern j for x, skip to next.
            if (sum(new_data_position) == 0) {
                next
            }

            if (num_observed_variables_y == 0) {
                response =  0.5 * (data$y[position, 1] + data$y[position, 2])
                if (num_observed_variables_x > 0) {
                    x_part = data$x[position, observed_variables_x]
                    ra_fit = lm(response ~ x_part)
                    new_data = cbind(rep(1, sum(new_data_position)),
                        data$x[new_data_position, observed_variables_x,
                            drop=FALSE])
                    new_data_fit = new_data %*% ra_fit$coefficients
                    theta_est = sum(new_data_fit) / n
                } else {
                    temp = mean(response)
                    theta_est = temp * sum(new_data_position) / n
                }
            } else {
                # we fit linear regression with the unobserved y variable
                # as response and observed x and observed y as covariates.
                response = data$y[position, missed_variables_y]
                x_part = cbind(data$x[position, observed_variables_x],
                            data$y[position, observed_variables_y])
                ra_fit = lm(response ~ x_part)
                new_data = cbind(rep(1, sum(new_data_position)),
                    data$x[new_data_position, observed_variables_x, drop=FALSE],
                    data$y[new_data_position, observed_variables_y, drop=FALSE])
                new_data_fit = new_data %*% ra_fit$coefficients
                theta_est = sum(0.5 * (new_data_fit + 
                    data$y[new_data_position, observed_variables_y])) / n
            }
            theta_estimate = theta_estimate + theta_est
        }
    }
    return(theta_estimate)
}

multiple_ra_indicator <- function(data, a) {
    # Use regression adjustment to estimate P(y_1 <= a, y_2 <= a).
    # We assume that there are two primary variables and this
    # function can be easily modified to extend to more primary
    # variables.
    # Args:
    #   data: a list containing data and missing information.
    #   a: a constant for converting y into binary variables.
    # Returns:
    #   theta_estimate: estimate of the proportion.
    n = length(data$A)
    theta_estimate = 0

    # pattern number for complete case
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    response_indicator = as.numeric(data$y[position_A_1, 1] <= a) *
                        as.numeric(data$y[position_A_1, 2] <= a)
    theta_estimate = theta_estimate + sum(response_indicator) / n

    for (i in 1 : data$num_pattern_y) {
        if (i == cp_num) {
            next
        }

        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        missed_variables_y = which(missing_pattern_y == F)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            more_x_observed = find_available_case_patterns(data,
                        num_observed_variables_x, observed_variables_x)

            position = data$A == cp_num & data$R %in% more_x_observed
            new_data_position = data$A == i & data$R == j

            if (sum(new_data_position) == 0) {
                next
            }

            if (num_observed_variables_y == 0) {
                response = as.numeric(data$y[position, 1] <= a) *
                    as.numeric(data$y[position, 2] <= a)

                if (num_observed_variables_x > 0) {
                    x_part = data$x[position, observed_variables_x]
                    ra_fit = glm(response ~ x_part, family=binomial)
                    new_data = cbind(rep(1, sum(new_data_position)),
                            data$x[new_data_position, observed_variables_x,
                                drop=FALSE])
                    new_data_fit = exp(-new_data %*% ra_fit$coefficients)
                    new_data_fit = 1 / (1 + new_data_fit)
                    theta_est = sum(new_data_fit) / n
                } else {
                    ra_fit = glm(response ~ 1, family=binomial)
                    new_data_fit = 1 / (1 + exp(-ra_fit$coefficients))
                    theta_est = new_data_fit * sum(new_data_position) / n
                }
            } else {
                # we observe one y, now fit a logistic regression for the
                # missed y with observed x and y as covariates.
                response = as.numeric(data$y[position, missed_variables_y] <= a)
                x_part = cbind(data$x[position, observed_variables_x],
                            data$y[position, observed_variables_y])
                ra_fit = glm(response ~ x_part, family=binomial)

                new_data = cbind(rep(1, sum(new_data_position)),
                    data$x[new_data_position, observed_variables_x, drop=FALSE],
                    data$y[new_data_position, observed_variables_y, drop=FALSE])
                new_data_fit = exp(-new_data %*% ra_fit$coefficients)
                new_data_fit = 1 / (1 + new_data_fit)
                observed_response = as.numeric(data$y[new_data_position,
                            observed_variables_y] <= a)
                theta_est = sum(new_data_fit * observed_response) / n
            }
            theta_estimate = theta_estimate + theta_est
        }
    }
    return(theta_estimate)
}

multiple_ipw <- function(data, fun) {
    # Use IPW to estimate the expectation of fun(data$y).
    # By default, logistic regression is used to estimate the odds function.
    # We assume that there are two primary variables.
    # Args:
    #     data: a list containing data and missing information.
    #     fun:  a function applied to the primary variables that
    #           takes two arguments.
    # Returns:
    #     theta_estimate: estimate by IPW.
    n = length(data$A)
    theta_estimate = 0

    # pattern number for complete case
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    theta_estimate = theta_estimate +
        sum(fun(data$y[position_A_1, 1], data$y[position_A_1, 2])) / n

    for (i in 1 : data$num_pattern_y) {
        if (i == cp_num) {
            next
        }
        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            more_x_observed = find_available_case_patterns(data,
                        num_observed_variables_x, observed_variables_x)

            position_0 = data$A == cp_num & data$R %in% more_x_observed
            position_1 = data$A == i & data$R == j
            if (sum(position_1) == 0) {
                next
            }
            position = position_0 | position_1

            # create binary response for the logistic regression to estimate
            # the odds function.
            response_ipw = rep(0, n)
            response_ipw[position_1] = 1

            response = fun(data$y[position_0, 1], data$y[position_0, 2])

            if (num_observed_variables_x + num_observed_variables_y == 0) {
                ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
                odds = exp(ipw_fit$coefficients)
                theta_est = odds * sum(response) / n
            } else {
                x_part = cbind(data$x[position, observed_variables_x],
                            data$y[position, observed_variables_y])
                ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)

                new_data = cbind(rep(1, sum(position_0)),
                        data$x[position_0, observed_variables_x, drop=FALSE],
                        data$y[position_0, observed_variables_y, drop=FALSE])
                odds = exp(new_data %*% ipw_fit$coefficients)
                theta_est = sum(odds * response) / n
            }
            theta_estimate = theta_estimate + theta_est
        }
    }
    return(theta_estimate)
}

multiple_ipw_sensitivity <- function(data, delta, fun) {
    # Sensitivity analysis with IPW with exponential tilting.
    # We assume that there are two primary variables.
    # Args:
    #     data: a list containing data and missing information.
    #     delta: parameter for sensitivity analysis.
    #     fun: a function applied to the primary variables that
    #          takes two arguments.
    # Returns:
    #     theta_estimate: estimate by IPW.
    n = length(data$A)
    theta_estimate = 0
    weights = rep(0, n)

    # pattern number for complete cases of primary variables.
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    theta_estimate = theta_estimate + 
            sum(fun(data$y[position_A_1, 1], data$y[position_A_1, 2]))
    weights[position_A_1] = 1

    for (i in 1 : data$num_pattern_y) {
        if (i == cp_num) {
            next
        }

        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        missed_variables_y = which(missing_pattern_y == F)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            more_x_observed = find_available_case_patterns(data,
                        num_observed_variables_x, observed_variables_x)

            position_0 = data$A == cp_num & data$R %in% more_x_observed
            position_1 = data$A == i & data$R == j
            if (sum(position_1) == 0) {
                next
            }
            position = position_0 | position_1
            response_ipw = rep(0, n)
            response_ipw[position_1] = 1

            response = fun(data$y[position_0, 1], data$y[position_0, 2])
            odds_depend_unobserved =
                exp(0.5 * (data$y[position_0, 1] +
                    data$y[position_0, 2]) * delta)

            if (num_observed_variables_y == 0) {
                if (num_observed_variables_x > 0) {
                    x_part = data$x[position, observed_variables_x]
                    ipw_fit = 
                        glm(response_ipw[position] ~ x_part, family=binomial)
                    new_data = cbind(rep(1, sum(position_0)),
                                    data$x[position_0, observed_variables_x])
                    odds = exp(new_data %*% ipw_fit$coefficients)
                    odds_sensitivity = odds * odds_depend_unobserved
                    theta_est = sum(odds_sensitivity * response)
                } else {
                    ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
                    odds = exp(ipw_fit$coefficients)
                    odds_sensitivity = odds * odds_depend_unobserved
                    theta_est = sum(odds_sensitivity * response)
                }
            } else {
                # we observe one y
                x_part = cbind(data$x[position, observed_variables_x],
                            data$y[position, observed_variables_y])
                ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)

                new_data = cbind(rep(1, sum(position_0)),
                                data$x[position_0, observed_variables_x],
                                data$y[position_0, observed_variables_y])
                odds_depend_on_unobserved =
                    exp(data$y[position_0, missed_variables_y] * delta)
                odds = exp(new_data %*% ipw_fit$coefficients)
                odds_sensitivity = odds_depend_on_unobserved * odds
                theta_est = sum(odds_sensitivity * response)
            }
            weights[position_0] = weights[position_0] + odds_sensitivity
            theta_estimate = theta_estimate + theta_est
        }
    }
    theta_estimate = theta_estimate / sum(weights)
    return(theta_estimate)
}

multiple_mr_average <- function(data) {
    # Use multiply robust to estimate the averages of primary variables.
    # Linear regressions are used for estimating the regression functions.
    # and logistic regressions are used for estimating the odds functions.
    # We assume that there are two primary variables and this
    # function can be easily modified to extend to more primary
    # variables.
    # Args:
    #   data: a list containing data and missing information.
    # Returns:
    #   theta_estimate: estimate of the averages.
    n = length(data$A)
    theta_estimate = 0

    # pattern number for complete case
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    theta_estimate = theta_estimate +
            sum(0.5 * (data$y[position_A_1, 1] + data$y[position_A_1, 2])) / n

    for (i in 1 : data$num_pattern_y) {
        if (i == cp_num) {
            next
        }
        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        missed_variables_y = which(missing_pattern_y == F)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            more_x_observed = find_available_case_patterns(data,
                    num_observed_variables_x, observed_variables_x)

            position_0 = data$A == cp_num & data$R %in% more_x_observed
            position_1 = data$A == i & data$R == j

            if (sum(position_1) == 0) {
                next
            }

            position = position_0 | position_1
            response_ipw = rep(0, n)
            response_ipw[position_1] = 1
            response = 0.5 * (data$y[position_0, 1] + data$y[position_0, 2])

            if (num_observed_variables_y == 0) {
                if (num_observed_variables_x > 0) {
                    x_part = data$x[position, observed_variables_x]
                    ipw_fit =
                        glm(response_ipw[position] ~ x_part, family=binomial)
                    new_data = cbind(rep(1, sum(position_0)),
                        data$x[position_0, observed_variables_x, drop=FALSE])
                    odds = exp(new_data %*% ipw_fit$coefficients)

                    x_part = data$x[position_0, observed_variables_x]
                    ra_fit = lm(response ~ x_part)
                    new_data = cbind(rep(1, sum(position_1)),
                        data$x[position_1, observed_variables_x, drop=FALSE])
                    temp = rep(0, n)
                    temp[position_0] = (response - fitted(ra_fit)) * odds
                    temp[position_1] = new_data %*% ra_fit$coefficients
                    theta_est = mean(temp)
                } else {
                    ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
                    odds = exp(ipw_fit$coefficients)
                    ra_fit = mean(response)
                    temp = rep(0, n)
                    temp[position_0] = (response - ra_fit) * odds
                    temp[position_1] = ra_fit
                    theta_est = mean(temp)
                }
            } else {
                # we observe one y
                x_part = cbind(data$x[position, observed_variables_x],
                    data$y[position, observed_variables_y, drop=FALSE])
                ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
                new_data = cbind(rep(1, sum(position_0)),
                        data$x[position_0, observed_variables_x, drop=FALSE],
                        data$y[position_0, observed_variables_y, drop=FALSE])
                odds = exp(new_data %*% ipw_fit$coefficients)

                response = data$y[position_0, missed_variables_y]
                x_part =
                    cbind(data$x[position_0, observed_variables_x, drop=FALSE],
                    data$y[position_0, observed_variables_y, drop=FALSE])

                ra_fit = lm(response ~ x_part)
                new_data = cbind(rep(1, sum(position_1)),
                        data$x[position_1, observed_variables_x, drop=FALSE],
                        data$y[position_1, observed_variables_y, drop=FALSE])
                new_data_fit = new_data %*% ra_fit$coefficients
                temp = rep(0, n)
                temp[position_0] =
                    (0.5 * (data$y[position_0, 1] + data$y[position_0, 2]) -
                        0.5 * (fitted(ra_fit) +
                        data$y[position_0, observed_variables_y])) * odds
                temp[position_1] =
                    0.5 * (data$y[position_1, observed_variables_y] +
                             new_data_fit)
                theta_est = mean(temp)
            }
            theta_estimate = theta_estimate + theta_est
        }
    }
    return(theta_estimate)
}

multiple_mr_indicator <- function(data, a) {
    # Use multiply robust to estimate P(y_1 <= a, y_2 <= a).
    # Logistic regressions are used to estimate both the regression and
    # odds functions.
    # We assume that there are two primary variables and this
    # function can be easily modified to extend to more primary
    # variables.
    # Args:
    #   data: a list containing data and missing information.
    #   a: a constant for converting y into binary variables.
    # Returns:
    #   theta_estimate: estimate of the averages.
    n = length(data$A)
    theta_estimate = 0

    # pattern number for complete case
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    response_indicator = as.numeric(data$y[position_A_1, 1] <= a) *
            as.numeric(data$y[position_A_1, 2] <= a)
    theta_estimate = theta_estimate + sum(response_indicator) / n

    for (i in 1 : data$num_pattern_y) {
        if (i == cp_num) {
            next
        }
        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        missed_variables_y = which(missing_pattern_y == F)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            more_x_observed = find_available_case_patterns(data,
                num_observed_variables_x, observed_variables_x)

            position_0 = data$A == cp_num & data$R %in% more_x_observed
            position_1 = data$A == i & data$R == j
            if (sum(position_1) == 0) {
                next
            }
            position = position_0 | position_1
            response_ipw = rep(0, n)
            response_ipw[position_1] = 1

            if (num_observed_variables_y == 0) {
                response_indicator = as.numeric(data$y[position_0, 1] <= a) *
                        as.numeric(data$y[position_0, 2] <= a)
                if (num_observed_variables_x > 0) {
                    x_part = data$x[position, observed_variables_x]
                    ipw_fit =
                        glm(response_ipw[position] ~ x_part, family=binomial)
                    new_data = cbind(rep(1, sum(position_0)),
                        data$x[position_0, observed_variables_x, drop=FALSE])
                    odds = exp(new_data %*% ipw_fit$coefficients)

                    x_part = data$x[position_0, observed_variables_x]
                    ra_fit = glm(response_indicator ~ x_part, family=binomial)
                    new_data = cbind(rep(1, sum(position_1)),
                        data$x[position_1, observed_variables_x, drop=FALSE])
                    temp = rep(0, n)
                    temp[position_0] =
                        (response_indicator - fitted(ra_fit)) * odds
                    temp[position_1] = exp(-new_data %*% ra_fit$coefficients)
                    temp[position_1] = 1 / (1 + temp[position_1])
                    theta_est = mean(temp)
                } else {
                    ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
                    odds = exp(ipw_fit$coefficients)

                    ra_fit = glm(response_indicator ~ 1, family=binomial)
                    ra_fit_val = 1 / (1 + exp(-ra_fit$coefficients))
                    temp = rep(0, n)
                    temp[position_0] =
                        (response_indicator - fitted(ra_fit)) * odds
                    temp[position_1] = ra_fit_val
                    theta_est = mean(temp)
                }
            } else {
                # we observe one y
                x_part =
                    cbind(data$x[position, observed_variables_x, drop=FALSE],
                        data$y[position, observed_variables_y, drop=FALSE])
                ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
                new_data = cbind(rep(1, sum(position_0)),
                    data$x[position_0, observed_variables_x, drop=FALSE],
                    data$y[position_0, observed_variables_y, drop=FALSE])
                odds = exp(new_data %*% ipw_fit$coefficients)

                response_indicator =
                    as.numeric(data$y[position_0, missed_variables_y] <= a)
                x_part = cbind(data$x[position_0, observed_variables_x],
                            data$y[position_0, observed_variables_y])

                ra_fit = glm(response_indicator ~ x_part, family=binomial)
                new_data = cbind(rep(1, sum(position_1)),
                        data$x[position_1, observed_variables_x, drop=FALSE],
                        data$y[position_1, observed_variables_y, drop=FALSE])
                new_data_fit = exp(-new_data %*% ra_fit$coefficients)
                new_data_fit = 1 / (1 + new_data_fit)
                temp = rep(0, n)

                response_ind = as.numeric(data$y[position_0, 1] <= a) *
                        as.numeric(data$y[position_0, 2] <= a)
                temp[position_0] = (response_ind - fitted(ra_fit) *
                    as.numeric(
                        data$y[position_0, observed_variables_y] <= a)) * odds
                temp[position_1] =
                    as.numeric(data$y[position_1, observed_variables_y] <= a) *
                        new_data_fit
                theta_est = mean(temp)
            }
            theta_estimate = theta_estimate + theta_est
        }
    }
    return(theta_estimate)
}

multiple_bootstrap <- function(data, method, n_B=1000, ...) {
    # Bootstrap to obtain the confidence intervals for regression adjustment,
    # IPW and multiply robust estimator.
    # Args:
    #     data: a list containing data and missing information.
    #     method: function for the estimators to perform bootstrap.
    #     n_B: number of bootstrap samples.
    #     ...: optional arguments for the method function.
    # Returns:
    #     bt_estimate: an array of bootstrap estimates.
    # Example usage:
    #   bt_ra_avg = multiple_bootstrap(data, multiple_ra_average)
    #   bt_ra_ind = multiple_bootstrap(data, multiple_ra_indicator, a=50)
    #   bt_ipw_avg = multiple_bootstrap(data, multiple_ipw,
    #                     fun=average)
    #   bt_ipw_ind = multiple_bootstrap(data, multiple_ipw,
    #                     fun=multiple_indicator)
    #   bt_mr = multiple_bootstrap(data, multiple_mr_average)
    n = nrow(data$y)
    theta_estimate = rep(0, n_B)

    for (i in 1:n_B) {
        if (i %% 100 == 0) {
            cat("now is the ", i, " iteration\n")
        }
        bt_index = sample(n, size=n, replace=T)
        data_bt = list()
        data_bt$x = data$x[bt_index, ]
        data_bt$y = data$y[bt_index, ]
        data_bt$A = data$A[bt_index]
        data_bt$R = data$R[bt_index]
        data_bt$missing_pattern_x = data$missing_pattern_x
        data_bt$num_pattern_x = data$num_pattern_x
        data_bt$missing_pattern_y = data$missing_pattern_y
        data_bt$num_pattern_y = data$num_pattern_y
        theta_estimate[i] = method(data_bt, ...)
    }
    return(theta_estimate)
}

ipw_regression <- function(data) {
    # Output weights for IPW, this can be used for estimating regression
    # parameters.
    # By default, logistic regression is used to estimate the odds function.
    # Args:
    #     data: a list containing data and missing information.
    # Return:
    #     weight: IPW weights for each data point.
    n = length(data$A)
    weight = rep(0, n)

    # pattern number for complete case
    cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

    # estimate the mean when A = 11
    position_A_1 = data$A == cp_num
    weight[position_A_1] = 1

    for (i in 1 : data$num_pattern_y) {
        if (i == cp_num) {
            next
        }
        missing_pattern_y = data$missing_pattern_y[i, ]
        observed_variables_y = which(missing_pattern_y == T)
        num_observed_variables_y = length(observed_variables_y)

        for (j in 1 : data$num_pattern_x) {
            missing_pattern_x = data$missing_pattern_x[j, ]
            observed_variables_x = which(missing_pattern_x == T)
            num_observed_variables_x = length(observed_variables_x)

            more_x_observed = find_available_case_patterns(data,
                        num_observed_variables_x, observed_variables_x)

            position_0 = data$A == cp_num & data$R %in% more_x_observed
            position_1 = data$A == i & data$R == j
            if (sum(position_1) == 0) {
                next
            }
            position = position_0 | position_1
            relative_position = data$A[position] == cp_num &
                data$R[position] %in% more_x_observed

            response_ipw = rep(0, n)
            response_ipw[position_1] = 1

            if (num_observed_variables_y == 0) {
                if (num_observed_variables_x > 0) {
                    x_part = data$x[position, observed_variables_x]
                    ipw_fit = 
                        glm(response_ipw[position] ~ x_part, family=binomial)
                    new_data = model.matrix(ipw_fit)[relative_position, ]
                    odds = exp(new_data %*% ipw_fit$coefficients)
                    weight[position_0] = weight[position_0] + odds
                } else {
                    ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
                    odds = exp(ipw_fit$coefficients)
                    weight[position_0] = weight[position_0] + odds
                }
            } else {
                # we observe at least one y
                x_part = as.matrix(
                    cbind(data$x[position, observed_variables_x],
                    data$y[position, observed_variables_y]))
                ipw_fit =
                    glm(response_ipw[position] ~ x_part, family=binomial)
                new_data = model.matrix(ipw_fit)[relative_position, ]
                odds = exp(new_data %*% ipw_fit$coefficients)
                weight[position_0] = weight[position_0] + odds
            }
        }
    }
    return(weight)
}

bootstrap_regression <- function(data, method, fm, n_B=1000, ...) {
    # Bootstrap to obtain the confidence intervals for regression parameters.
    # Args:
    #    data: a list containing data and missing information.
    #    method: function for performing regression, for example, lm or glm.
    #    fm: formula for the regression.
    #    n_B: number of bootstrap times.
    #    ...: optional parameters for the method function.
    # Returns:
    #    beta_estimate: a matrix of bootstrap estimate of parameters.
    n = length(data$A)
    num_beta = length(attr(terms(fm), "term.labels")) + 1
    beta_estimate = matrix(0, nrow=n_B, ncol=num_beta)

    for (i in 1:n_B) {
        if (i %% 100 == 0) {
            cat("now is the ", i, " iteration\n")
        }
        bt_index = sample(n, size=n, replace=T)
        data_bt = list()
        data_bt$x = as.matrix(reg_data$x[bt_index, ])
        data_bt$y = as.matrix(reg_data$y[bt_index, ])
        data_bt$A = reg_data$A[bt_index]
        data_bt$R = reg_data$R[bt_index]
        data_bt$missing_pattern_x = reg_data$missing_pattern_x
        data_bt$num_pattern_x = reg_data$num_pattern_x
        data_bt$missing_pattern_y = reg_data$missing_pattern_y
        data_bt$num_pattern_y = reg_data$num_pattern_y

        weight = ipw_regression(data_bt)
        y_df = as.data.frame(data_bt$y)
        y_df$weight = weight
        beta_estimate[i, ] =
            method(formula=fm, data=y_df, weights=weight, ...)$coefficients
    }
    return(beta_estimate)
}
