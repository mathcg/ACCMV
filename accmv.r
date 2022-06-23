library(mgcv)
setwd('P:/CCW_Local_Diabetes_Statistics_Zhao/Analysis_Data_Deidentified/Zhao_Diabetes_Methods/gang/')
a1c_data = read.csv("a1c_dat.csv")[, -1]
a1c_first_year = a1c_data[, 1:6]
n = nrow(a1c_first_year)

single_data_preparation <- function(x, y) {
    # Prepare the data when there is single primary variable y.
    # Args:
    #    x: secondary variables.
    #    y: primary variables.
    # Returns:
    #    data:
    #    x: x
    #    y: y
    #    R: missing_indicator_x
    #    A: missing_indicator_y
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

single_data = single_data_preparation(a1c_first_year[, 2:5], 
                                      a1c_first_year[, 6])

single_regression_adjustment_original <- function(data) {
    n = length(data$y)
    theta_estimate = 0

    data_x_y = data.frame(data$x, y=data$y)
    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(data$y[position_A_1]) / n
  
    for (i in 1 : data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)

        if (num_observed_variables == 0) {
        more_x_observed = 1 : data$num_pattern_x
        } else if (num_observed_variables == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables])
        } else {
        more_x_observed = which(apply(
            data$missing_pattern_x[, observed_variables], 
            1, function(v) sum(v) == num_observed_variables))
        }
        position = data$A == 1 & data$R %in% more_x_observed
        new_data_position = data$A == 0 & data$R == i
        if (num_observed_variables > 0) {
            # lm will automatically use available cases to fit the linear regression.
            ra_fm = as.formula(paste("y ~ ", 
                                paste(colnames(data$x)[observed_variables], 
                                collapse = '+'), sep=""))
            model_matrix_x = model.matrix.lm(ra_fm, data=data_x_y,
                                             na.action=na.pass)
            ra_fit = lm(as.formula(ra_fm), data=data_x_y, subset=position)
            new_data = cbind(rep(1, sum(new_data_position)),
                            data$x[new_data_position, observed_variables])
            new_data_fit = new_data %*% ra_fit$coefficients
            theta_est = sum(new_data_fit) / n
        } else {
            temp = mean(data$y[position])
            theta_est = temp * sum(new_data_position) / n
        }
        theta_estimate = theta_estimate + theta_est
    }
    return(theta_estimate)
}

single_indicator <- function(x) as.numeric(x <= 7)
single_regression_adjustment <- function(data, fun=identity, binary=F) {
    n = length(data$y)
    theta_estimate = 0

    # Create a data frame for glm.
    data_x_y = data.frame(data$x, y=fun(data$y))
    colnames(data_x_y) = c(colnames(data$x), "y")

    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(fun(data$y[position_A_1])) / n

    family = if (binary) binomial else gaussian
  
    pattern_est = rep(0, data$num_pattern_x)
    for (i in 1 : data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)
        new_data_position = data$A == 0 & data$R == i

        if (num_observed_variables > 0) {
            # lm will use available cases to fit the linear regression.
            ra_fm = as.formula(paste("y ~ ", 
                                paste(colnames(data$x)[observed_variables], 
                                collapse = '+'), sep=""))
            model_matrix_x = model.matrix.lm(ra_fm, data=data_x_y,
                                             na.action=na.pass)
            ra_fit = glm(ra_fm, data=data_x_y, family=family)
            new_data = data.frame(model_matrix_x[new_data_position, ])
            colnames(new_data) = names(ra_fit$coefficients)
            if (!binary) {
                new_data_fit = predict(ra_fit, new_data)
            } else {
                new_data_fit = predict(ra_fit, new_data, type="response")
            }
            theta_est = sum(new_data_fit) / n
        } else {
            temp = mean(fun(data$y[position_A_1]))
            theta_est = temp * sum(new_data_position) / n
        }
        pattern_est[i] = theta_est
        theta_estimate = theta_estimate + theta_est
    }
    return(theta_estimate)
}

single_ipw <- function(data, fun=identity) {
  n = length(data$y)
  theta_estimate = 0
  
  # estimate the mean when A = 1
  position_A_1 = data$A == 1
  theta_estimate = theta_estimate + sum(fun(data$y[position_A_1])) / n
  
  for (i in 1:data$num_pattern_x) {
    missing_pattern_x = data$missing_pattern_x[i, ]
    observed_variables = which(missing_pattern_x == T)
    num_observed_variables = length(observed_variables)

    # Find available cases
    if (num_observed_variables == 0) {
      more_x_observed = 1:data$num_pattern_x
    } else if (num_observed_variables == 1) {
      more_x_observed = which(data$missing_pattern_x[, observed_variables])
    } else {
      more_x_observed = which(apply(
        data$missing_pattern_x[, observed_variables], 
        1, function(v) sum(v) == num_observed_variables))
    }
    position_0 = data$A == 1 & data$R %in% more_x_observed
    position_1 = data$A == 0 & data$R == i
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
    if (num_observed_variables == 0) {
      more_x_observed = 1:data$num_pattern_x
    } else if (num_observed_variables == 1) {
      more_x_observed = which(data$missing_pattern_x[, observed_variables])
    } else {
      more_x_observed = which(apply(
        data$missing_pattern_x[, observed_variables], 
        1, function(v) sum(v) == num_observed_variables))
    }
    position_0 = data$A == 1 & data$R %in% more_x_observed
    position_1 = data$A == 0 & data$R == i
    position = position_0 | position_1
    response = rep(0, n)
    response[position_1] = 1
    odds_part_unobserved = exp((data$y[position_0] - 7) * delta)
    
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
  theta_estimate = theta_estimate / sum(weights)
  return(theta_estimate)
}

single_multiply_robust <- function(data, fun=identity, binary=F) {
    n = length(data$y)
    theta_estimate = 0
  
    # Create a data frame for regression adjustment.
    data_x_y = data.frame(data$x, y=fun(data$y))
    colnames(data_x_y) = c(colnames(data$x), "y")
    
    # estimate the mean when A = 1
    position_A_1 = data$A == 1
    theta_estimate = theta_estimate + sum(fun(data$y[position_A_1])) / n
  
    family = if (binary) binomial else gaussian
    for (i in 1:data$num_pattern_x) {
        missing_pattern_x = data$missing_pattern_x[i, ]
        observed_variables = which(missing_pattern_x == T)
        num_observed_variables = length(observed_variables)

        if (num_observed_variables == 0) {
            more_x_observed = 1:data$num_pattern_x
        } else if (num_observed_variables == 1) {
            more_x_observed = which(data$missing_pattern_x[, observed_variables])
        } else {
            more_x_observed = which(apply(
                data$missing_pattern_x[, observed_variables], 
                1, function(v) sum(v) == num_observed_variables))
        }
        position_0 = data$A == 1 & data$R %in% more_x_observed
        position_1 = data$A == 0 & data$R == i
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
            ra_fit = glm(ra_fm, data=data_x_y, family=family)
            new_data_ra = data.frame(model_matrix_x[position_1, ])
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
        theta_estimate = theta_estimate + theta_est
    }
    return(theta_estimate)
}

single_bootstrap <- function(data, method, n_B=1000, ...) {
    # Bootstrap to obtain the confidence intervals.
    args = list(...)
    n = length(data$y)
    theta_estimate = rep(0, n_B)
    if (!is.null(args$fun)) {
      fun = args$fun
    } else {
      fun = identity
    }

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
        if (!is.null(args$binary)) {
          theta_estimate[i] = method(data_bt, fun, args$binary)
        } else {
          theta_estimate[i] = method(data_bt, fun)
        }
    }
    return(theta_estimate)
}

single_bootstrap_sensitivity <- function(data, delta, n_B = 1000, 
                                        fun=identity) {
  n = length(data$y)
  theta_estimate = rep(0, n_B)
  
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
    theta_estimate[i] = single_ipw_sensitivity(data_bt, delta,
                        fun=fun)
  }
  return(theta_estimate)
}

# sensitivity analysis for E[Y_4]
cc_num = sum(single_data$A == 1)
theta_estimate_cc = mean(single_data$y[single_data$A == 1])
cc_sd = sd(single_data$y[single_data$A == 1]) / sqrt(cc_num)
delta = seq(-0.5, 0.5, 0.1)
theta_estimate_ipw_sa = rep(0, length(delta))
theta_estimate_ipw_sa_sd = rep(0, length(delta))
theta_estimate_ipw_sa_l = rep(0, length(delta))
theta_estimate_ipw_sa_u = rep(0, length(delta))
for (i in 1 : length(delta)) {
  theta_estimate_ipw_sa[i] = single_ipw_sensitivity(single_data, delta[i])
  theta_estimate_ipw_sa_bt = single_bootstrap_sensitivity(single_data, 
                          delta[i])
  temp_ci = quantile(theta_estimate_ipw_sa_bt, probs = c(0.025, 0.975))
  theta_estimate_ipw_sa_sd[i] = sd(theta_estimate_ipw_sa_bt)
  theta_estimate_ipw_sa_l[i] = temp_ci[1]
  theta_estimate_ipw_sa_u[i] = temp_ci[2]
}

delta_0 = c(delta, 0)
sa_estimate = c(theta_estimate_ipw_sa, theta_estimate_cc)
sa_ci_l = c(theta_estimate_ipw_sa_l, theta_estimate_cc - 1.96 * cc_sd)
sa_ci_u = c(theta_estimate_ipw_sa_u, theta_estimate_cc + 1.96 * cc_sd)
method = c(rep("IPW perturbed", length(delta)),  "complete-case")

theta_df = data.frame(delta=delta_0, sa_estimate = sa_estimate,
                      sa_ci_l = sa_ci_l,
                      sa_ci_u = sa_ci_u,
                      method=method)


# multiple primary variables
# now compute E[Y_4 Y_5]

multiple_data_preparation <- function(x, y) {
    # Prepare the data when there are multiple primary variables y.
    # Args:
    #   x: secondary variables
    #   y: primary variables
    # Returns:
    #   data:
    missing_pattern_x = uniquecombs(!is.na(x))
    num_pattern_x = nrow(missing_pattern_x)
    missing_indicator_x = attr(missing_pattern_x, "index")

    missing_pattern_y = uniquecombs(!is.na(y))
    num_pattern_y = nrow(missing_pattern_y)
    missing_indicator_y = attr(missing_pattern_y, "index")

    data = list()
    data$y = y
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

multiple_data = multiple_data_preparation(a1c_first_year[, 2:4],
                                      a1c_first_year[, 5:6])

# Estimate E[Y_3 + Y_4] / 2
multiple_regression_adjustment_average <- function(data) {
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
      
      # Find out all available cases for the auxillary variables x.
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      # Complete cases for primary variables and available cases for auxillary variable.
      position = data$A == cp_num & data$R %in% more_x_observed
      new_data_position = data$A == i & data$R == j
      
      if (num_observed_variables_y == 0) {
        response =  0.5 * (data$y[position, 1] + data$y[position, 2])
        
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ra_fit = lm(response ~ x_part)
          new_data = cbind(rep(1, sum(new_data_position)),
                           data$x[new_data_position, observed_variables_x])
          new_data_fit = new_data %*% ra_fit$coefficients
          theta_est = sum(new_data_fit) / n
        } else {
          temp = mean(response)
          theta_est = temp * sum(new_data_position) / n
        }
      } else {
        # we observe one y variable and fit linear regression on the unobserved y variable with
        # observed x and observed y as covariates.
        response = data$y[position, missed_variables_y]
        x_part = cbind(data$x[position, observed_variables_x], 
                       data$y[position, observed_variables_y])
        ra_fit = lm(response ~ x_part)
        
        new_data = cbind(rep(1, sum(new_data_position)),
                         data$x[new_data_position, observed_variables_x],
                         data$y[new_data_position, observed_variables_y])
        new_data_fit = new_data %*% ra_fit$coefficients
        theta_est = sum(0.5 * (new_data_fit + data$y[new_data_position, observed_variables_y])) / n
      }
      theta_estimate = theta_estimate + theta_est
    }
  }
  return(theta_estimate)
}

# estimate P(Y_3 <= 7, Y_4 <= 7)
multiple_regression_adjustment_indicator <- function(data) {
  n = length(data$A)
  theta_estimate = 0
  
  # pattern number for complete case
  cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

  # estimate the mean when A = 11
  position_A_1 = data$A == cp_num
  response_indicator = as.numeric(data$y[position_A_1, 1] <= 7) * 
                      as.numeric(data$y[position_A_1, 2] <= 7)
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
      
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      position = data$A == cp_num & data$R %in% more_x_observed
      new_data_position = data$A == i & data$R == j
      
      if (num_observed_variables_y == 0) {
        response = as.numeric(data$y[position, 1] <= 7) * as.numeric(data$y[position, 2] <= 7)
        
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ra_fit = glm(response ~ x_part, family=binomial)
          new_data = cbind(rep(1, sum(new_data_position)),
                           data$x[new_data_position, observed_variables_x])
          new_data_fit = exp(-new_data %*% ra_fit$coefficients)
          new_data_fit = 1 / (1 + new_data_fit)
          theta_est = sum(new_data_fit) / n
        } else {
          ra_fit = glm(response ~ 1, family=binomial)
          new_data_fit = 1 / (1 + exp(-ra_fit$coefficients))
          theta_est = new_data_fit * sum(new_data_position) / n
        }
      } else {
        # we observe one y, now fit a logistic regression for the missed y with 
        # observed x and y as covariates.
        response = as.numeric(data$y[position, missed_variables_y] <= 7)
        x_part = cbind(data$x[position, observed_variables_x], 
                       data$y[position, observed_variables_y])
        ra_fit = glm(response ~ x_part, family=binomial)
        
        new_data = cbind(rep(1, sum(new_data_position)),
                         data$x[new_data_position, observed_variables_x],
                         data$y[new_data_position, observed_variables_y])
        new_data_fit = exp(-new_data %*% ra_fit$coefficients)
        new_data_fit = 1 / (1 + new_data_fit)
        observed_response = as.numeric(data$y[new_data_position, observed_variables_y] <= 7)
        theta_est = sum(new_data_fit * observed_response) / n
      }
      theta_estimate = theta_estimate + theta_est
    }
  }
  return(theta_estimate)
}

average <- function(x1, x2) 0.5 * (x1 + x2)
indicator <- function(x1, x2) {
  as.numeric(x1 <= 7) * as.numeric(x2 <= 7)
}

multiple_ipw <- function(data, fun) {
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
    missed_variables_y = which(missing_pattern_y == F)
    num_observed_variables_y = length(observed_variables_y)
    
    for (j in 1 : data$num_pattern_x) {
      missing_pattern_x = data$missing_pattern_x[j, ]
      observed_variables_x = which(missing_pattern_x == T)
      num_observed_variables_x = length(observed_variables_x)
      
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      position_0 = data$A == cp_num & data$R %in% more_x_observed
      position_1 = data$A == i & data$R == j
      position = position_0 | position_1

      # create binary response for the logistic regression to estimate the odds function
      response_ipw = rep(0, n)
      response_ipw[position_1] = 1

      response = fun(data$y[position_0, 1], data$y[position_0, 2])
      
      if (num_observed_variables_y == 0) {
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
          new_data = cbind(rep(1, sum(position_0)),
                           data$x[position_0, observed_variables_x])
          odds = exp(new_data %*% ipw_fit$coefficients)
          theta_est = sum(odds * response) / n
        } else {
          ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
          odds = exp(ipw_fit$coefficients)
          theta_est = odds * sum(response) / n
        }
      } else {
        # we observe one y
        x_part = cbind(data$x[position, observed_variables_x], 
                       data$y[position, observed_variables_y])
        ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
        
        new_data = cbind(rep(1, sum(position_0)),
                         data$x[position_0, observed_variables_x],
                         data$y[position_0, observed_variables_y])
        odds = exp(new_data %*% ipw_fit$coefficients)
        theta_est = sum(odds * response) / n
      }
      theta_estimate = theta_estimate + theta_est
    }
  }
  return(theta_estimate)
}

multiple_ipw_sensitivity_analysis <- function(data, delta, fun) {
  n = length(data$A)
  theta_estimate = 0
  weights = rep(0, n)
  
  # pattern number for complete case
  cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

  # estimate the mean when A = 11
  position_A_1 = data$A == cp_num
  theta_estimate = theta_estimate + 
        sum(fun(data$y[position_A_1, 1], data$y[position_A_1, 2]))
        # sum(0.5 * (data$y[position_A_1, 1] + data$y[position_A_1, 2]))
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
      
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      position_0 = data$A == 1 & data$R %in% more_x_observed
      position_1 = data$A == i & data$R == j
      position = position_0 | position_1
      response_ipw = rep(0, n)
      response_ipw[position_1] = 1
      # response = 0.5 * (data$y[position_0, 1] + data$y[position_0, 2])
      response = fun(data$y[position_0, 1], data$y[position_0, 2])
      odds_depend_unobserved = exp((data$y[position_0, 1] + data$y[position_0, 2] - 14) * delta)
      
      if (num_observed_variables_y == 0) {
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
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
        odds_depend_on_unobserved = exp((data$y[position_0, missed_variables_y] - 7) * delta)
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

multiple_multiply_robust_average <- function(data) {
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
      
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      position_0 = data$A == cp_num & data$R %in% more_x_observed
      position_1 = data$A == i & data$R == j
      position = position_0 | position_1
      response_ipw = rep(0, n)
      response_ipw[position_1] = 1
      response = 0.5 * (data$y[position_0, 1] + data$y[position_0, 2])
      
      if (num_observed_variables_y == 0) {
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
          new_data = cbind(rep(1, sum(position_0)),
                           data$x[position_0, observed_variables_x])
          odds = exp(new_data %*% ipw_fit$coefficients)
          
          x_part = data$x[position_0, observed_variables_x]
          ra_fit = lm(response ~ x_part)
          new_data = cbind(rep(1, sum(position_1)),
                           data$x[position_1, observed_variables_x])
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
                       data$y[position, observed_variables_y])
        ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
        new_data = cbind(rep(1, sum(position_0)),
                         data$x[position_0, observed_variables_x],
                         data$y[position_0, observed_variables_y])
        odds = exp(new_data %*% ipw_fit$coefficients)
        
        response = data$y[position_0, missed_variables_y]
        x_part = cbind(data$x[position_0, observed_variables_x], 
                       data$y[position_0, observed_variables_y])
        
        ra_fit = lm(response ~ x_part)
        new_data = cbind(rep(1, sum(position_1)),
                         data$x[position_1, observed_variables_x],
                         data$y[position_1, observed_variables_y])
        new_data_fit = new_data %*% ra_fit$coefficients
        temp = rep(0, n)
        temp[position_0] = (0.5 * (data$y[position_0, 1] + data$y[position_0, 2]) - 
                            0.5 * (fitted(ra_fit) + data$y[position_0, observed_variables_y])) * odds
        temp[position_1] = 0.5 * (data$y[position_1, observed_variables_y] + new_data_fit)
        theta_est = mean(temp)
      }
      theta_estimate = theta_estimate + theta_est
    }
  }
  return(theta_estimate)
}

multiple_multiply_robust_indicator <- function(data) {
  n = length(data$A)
  theta_estimate = 0
  
  # pattern number for complete case
  cp_num = which(apply(data$missing_pattern_y, 1, prod) == 1)

  # estimate the mean when A = 11
  position_A_1 = data$A == cp_num
  response_indicator = as.numeric(data$y[position_A_1, 1] <= 7) * 
        as.numeric(data$y[position_A_1, 2] <= 7)
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
      
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      position_0 = data$A == cp_num & data$R %in% more_x_observed
      position_1 = data$A == i & data$R == j
      position = position_0 | position_1
      response_ipw = rep(0, n)
      response_ipw[position_1] = 1
      
      if (num_observed_variables_y == 0) {
        response_indicator = as.numeric(data$y[position_0, 1] <= 7) * 
                  as.numeric(data$y[position_0, 2] <= 7)
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
          new_data = cbind(rep(1, sum(position_0)),
                           data$x[position_0, observed_variables_x])
          odds = exp(new_data %*% ipw_fit$coefficients)
          
          x_part = data$x[position_0, observed_variables_x]
          ra_fit = glm(response_indicator ~ x_part, family=binomial)
          new_data = cbind(rep(1, sum(position_1)),
                           data$x[position_1, observed_variables_x])
          temp = rep(0, n)
          temp[position_0] = (response_indicator - fitted(ra_fit)) * odds
          temp[position_1] = exp(-new_data %*% ra_fit$coefficients)
          temp[position_1] = 1 / (1 + temp[position_1])
          theta_est = mean(temp)
        } else {
          ipw_fit = glm(response_ipw[position] ~ 1, family=binomial)
          odds = exp(ipw_fit$coefficients)
          
          ra_fit = glm(response_indicator ~ 1, family=binomial)
          ra_fit_val = 1 / (1 + exp(-ra_fit$coefficients))
          temp = rep(0, n)
          temp[position_0] = (response_indicator - fitted(ra_fit)) * odds
          temp[position_1] = ra_fit_val
          theta_est = mean(temp)
        }
      } else {
        # we observe one y
        x_part = cbind(data$x[position, observed_variables_x], 
                       data$y[position, observed_variables_y])
        ipw_fit = glm(response_ipw[position] ~ x_part, family=binomial)
        new_data = cbind(rep(1, sum(position_0)),
                         data$x[position_0, observed_variables_x],
                         data$y[position_0, observed_variables_y])
        odds = exp(new_data %*% ipw_fit$coefficients)
        
        response_indicator = as.numeric(data$y[position_0, missed_variables_y] <= 7)
        x_part = cbind(data$x[position_0, observed_variables_x], 
                       data$y[position_0, observed_variables_y])
        
        ra_fit = glm(response_indicator ~ x_part, family=binomial)
        new_data = cbind(rep(1, sum(position_1)),
                         data$x[position_1, observed_variables_x],
                         data$y[position_1, observed_variables_y])
        new_data_fit = exp(-new_data %*% ra_fit$coefficients)
        new_data_fit = 1 / (1 + new_data_fit)
        temp = rep(0, n)
        
        response_ind = as.numeric(data$y[position_0, 1] <= 7) *
                   as.numeric(data$y[position_0, 2] <= 7)
        temp[position_0] = (response_ind - fitted(ra_fit) * 
                              as.numeric(data$y[position_0, observed_variables_y] <= 7)) * odds
        temp[position_1] = as.numeric(data$y[position_1, observed_variables_y] <= 7) * new_data_fit
        theta_est = mean(temp)
      }
      theta_estimate = theta_estimate + theta_est
    }
  }
  return(theta_estimate)
}

multiply_bootstrap <- function(data, f, n_B=1000, ...) {
  args = list(...)
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
    if (!is.null(args$fun)) {
      theta_estimate[i] = f(data_bt, args$fun)
    } else {
      theta_estimate[i] = f(data_bt)
    }
  }
  return(theta_estimate)
}

start_time = Sys.time()
result = multiply_bootstrap(data, multiple_ipw, fun=average)
end_time = Sys.time()

start_time = Sys.time()
result = multiply_bootstrap(data, multiple_regression_adjustment_average)
end_time = Sys.time()

start_time = Sys.time()
result = multiply_bootstrap(data, multiple_regression_adjustment_indicator)
end_time = Sys.time()

start_time = Sys.time()
result = multiply_bootstrap(data, multiple_ipw, fun=indicator)
end_time = Sys.time()

start_time = Sys.time()
result = multiply_bootstrap(data, multiple_multiply_robust_average)
end_time = Sys.time()

start_time = Sys.time()
result = multiply_bootstrap(data, multiple_multiply_robust_indicator)
end_time = Sys.time()

multiple_bootstrap_sensitivity_analysis <- function(data, delta, fun, n_B=1000) {
  n = length(data$A)
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
    theta_estimate[i] = multiple_ipw_sensitivity_analysis(data_bt, delta, fun)
  }
  return(theta_estimate)
}

# complete case analysis for P(Y_3 <= 7, Y_4 <= 7)
cp_num = which(apply(multiple_data$missing_pattern_y, 1, prod) == 1)
cc_position = multiple_data$A == cp_num
mt_response_indicator = as.numeric(multiple_data$y[cc_position, 1] <= 7) * 
            as.numeric(multiple_data$y[cc_position, 2] <= 7)
prop_mt_estimate_cc = mean(mt_response_indicator)
prop_mt_estimate_cc_sd = sd(mt_response_indicator) / sqrt(length(mt_response_indicator))
prop_mt_ci = c(prop_mt_estimate_cc - 1.96 * prop_mt_estimate_cc_sd, prop_mt_estimate_cc + 1.96 * prop_mt_estimate_cc_sd)

delta = seq(-0.5, 0.5, 0.1)
prop_mt_estimate_ipw_sa = rep(0, length(delta))
prop_mt_estimate_ipw_sa_sd = rep(0, length(delta))
prop_mt_estimate_ipw_sa_l = rep(0, length(delta))
prop_mt_estimate_ipw_sa_u = rep(0, length(delta))
for (i in 1 : length(delta)) {
  prop_mt_estimate_ipw_sa[i] = multiple_ipw_sensitivity_analysis(
            multiple_data, delta[i], indicator)
  prop_mt_estimate_ipw_sa_bt = multiple_bootstrap_sensitivity_analysis(
            multiple_data, delta[i], indicator)
  temp_ci = quantile(prop_mt_estimate_ipw_sa_bt, probs = c(0.025, 0.975))
  prop_mt_estimate_ipw_sa_sd[i] = sd(prop_mt_estimate_ipw_sa_bt)
  prop_mt_estimate_ipw_sa_l[i] = temp_ci[1]
  prop_mt_estimate_ipw_sa_u[i] = temp_ci[2]
}

delta_0 = c(delta, 0)
sa_estimate = c(prop_mt_estimate_ipw_sa, prop_mt_estimate_cc)
sa_ci_l = c(prop_mt_estimate_ipw_sa_l, prop_mt_ci[1])
sa_ci_u = c(prop_mt_estimate_ipw_sa_u, prop_mt_ci[2])
method = c(rep("IPW perturbed", length(delta)),  "complete-case")

prop_mt_theta_df = data.frame(delta=delta_0, sa_estimate = sa_estimate,
                      sa_ci_l = sa_ci_l,
                      sa_ci_u = sa_ci_u,
                      method=method)

# estimate 0.5 * (E[Y_3 + Y_4])
mt_response_average = 0.5 * (multiple_data$y[cc_position, 1] + 
          multiple_data$y[cc_position, 2])
average_mt_estimate_cc = mean(mt_response_average)
average_mt_estimate_cc_sd = sd(mt_response_average) / sqrt(length(mt_response_average))
average_mt_ci = c(average_mt_estimate_cc - 1.96 * average_mt_estimate_cc_sd, average_mt_estimate_cc + 1.96 * average_mt_estimate_cc_sd)

delta = seq(-0.5, 0.5, 0.1)
average_mt_estimate_ipw_sa = rep(0, length(delta))
average_mt_estimate_ipw_sa_sd = rep(0, length(delta))
average_mt_estimate_ipw_sa_l = rep(0, length(delta))
average_mt_estimate_ipw_sa_u = rep(0, length(delta))
for (i in 1:length(delta)) {
  average_mt_estimate_ipw_sa[i] = multiple_ipw_sensitivity_analysis(
                  multiple_data, delta[i], average)
  average_mt_estimate_ipw_sa_bt = multiple_bootstrap_sensitivity_analysis(
                  multiple_data, delta[i], average)
  temp_ci = quantile(average_mt_estimate_ipw_sa_bt, probs = c(0.025, 0.975))
  average_mt_estimate_ipw_sa_sd[i] = sd(average_mt_estimate_ipw_sa_bt)
  average_mt_estimate_ipw_sa_l[i] = temp_ci[1]
  average_mt_estimate_ipw_sa_u[i] = temp_ci[2]
}

delta_0 = c(delta, 0)
sa_estimate = c(average_mt_estimate_ipw_sa, average_mt_estimate_cc)
sa_ci_l = c(average_mt_estimate_ipw_sa_l, average_mt_ci[1])
sa_ci_u = c(average_mt_estimate_ipw_sa_u, average_mt_ci[2])
method = c(rep("IPW perturbed", length(delta)),  "complete-case")

average_mt_theta_df = data.frame(delta=delta_0, sa_estimate = sa_estimate,
                      sa_ci_l = sa_ci_l,
                      sa_ci_u = sa_ci_u,
                      method=method)

cbp2 <- c("red", "black")
g2_single_average <- ggplot(theta_df) + 
  geom_point(aes(x=delta, y=sa_estimate, color=method)) + 
  geom_errorbar(aes(ymax=sa_ci_u, ymin=sa_ci_l, x=delta, color=method), width=0.08) + 
  scale_colour_manual(values=cbp2) + 
  xlab(expression(delta)) + 
  ylab(expression(paste("E[", Y[4], "]", sep=""))) + 
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.2)) + 
  ggtitle(expression(paste("(a) Sensitivity analysis for E[", Y[4], "]", sep=""))) + 
  theme_bw() + 
  theme(
      plot.title = element_text(hjust=0.5, size=20),
      legend.text = element_text(size=18),
      legend.title = element_text(size=20),
      axis.text=element_text(size=16),
      axis.title=element_text(size=18, face="bold"),
      legend.position="none"
  )

g2_mt_prop <- ggplot(prop_mt_theta_df) + 
  geom_point(aes(x=delta, y=sa_estimate, color=method)) + 
  geom_errorbar(aes(ymax=sa_ci_u, ymin=sa_ci_l, x=delta, color=method), 
        width=0.08) + 
  scale_colour_manual(values=cbp2) + 
  xlab(expression(delta)) + 
  ylab(expression(paste("P(", Y[3], " <= 7,", Y[4], " <= 7)", sep=""))) + 
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.2)) + 
  ggtitle(expression(paste("(b) Sensitivity analysis for P(", Y[3], " <= 7,", Y[4], " <= 7)", sep=""))) + 
  theme_bw() + 
  theme(
      plot.title = element_text(hjust=0.5, size=20),
      legend.text = element_text(size=18),
      legend.title = element_text(size=20),
      axis.text=element_text(size=16),
      axis.title=element_text(size=18, face="bold"),
      legend.position="none"
  )

g2_mt_mean <- ggplot(average_mt_theta_df) + 
  geom_point(aes(x=delta, y=sa_estimate, color=method)) + 
  geom_errorbar(aes(ymax=sa_ci_u, ymin=sa_ci_l, x=delta, color=method), 
        width=0.08) + 
  scale_colour_manual(values=cbp2) + 
  xlab(expression(delta)) + 
  ylab(expression(paste("(E[", Y[3], "] + E[", Y[4], "]) / 2", sep=""))) + 
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.2)) + 
  ggtitle(expression(paste("(c) Sensitivity analysis for (E[", Y[3], "] + E[", Y[4], "]) / 2", sep=""))) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust=0.5, size=20),
    axis.text=element_text(size=16),
    axis.title=element_text(size=18, face="bold"),
    legend.text = element_text(size=18),
    legend.title = element_text(size=20),
    # legend.position = c(0.2, 0.9)
    legend.position = "bottom"
  )

ggarrange(g2_single_average, g2_mt_prop, g2_mt_mean, ncol=3, 
  common.legend=T, legend="bottom")
ggsave(file=paste("P:/CCW_Local_Diabetes_Statistics_Zhao/",
      "Analysis_Data_Deidentified/Zhao_Diabetes_Methods/gang",
      "/Sensitivity_analysis_summary_measure.pdf", sep=""), width=20, height=6)
# Marginal parametric models

x = a1c_first_year[, 2:3]
y = a1c_first_year[, 4:6]
mpm_data = multiple_data_preparation(x, y)

# marginal parametric model, now try to fit a regression model with Y_5 ~ Y_4 + Y_3
ipw_marginal_parametric_model <- function(data) {
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
    missed_variables_y = which(missing_pattern_y == F)
    num_observed_variables_y = length(observed_variables_y)
    
    for (j in 1 : data$num_pattern_x) {
      missing_pattern_x = data$missing_pattern_x[j, ]
      observed_variables_x = which(missing_pattern_x == T)
      num_observed_variables_x = length(observed_variables_x)
      
      if (num_observed_variables_x == 0) {
        more_x_observed = 1 : data$num_pattern_x
      } else if (num_observed_variables_x == 1) {
        more_x_observed = which(data$missing_pattern_x[, observed_variables_x])
      } else {
        more_x_observed = which(apply(
          data$missing_pattern_x[, observed_variables_x], 
          1, function(v) sum(v) == num_observed_variables_x))
      }
      
      position_0 = data$A == cp_num & data$R %in% more_x_observed
      position_1 = data$A == i & data$R == j
      position = position_0 | position_1
      relative_position = data$A[position] == cp_num & data$R[position] %in% more_x_observed

      response = rep(0, n)
      response[position_1] = 1
      
      if (num_observed_variables_y == 0) {
        if (num_observed_variables_x > 0) {
          x_part = data$x[position, observed_variables_x]
          ipw_fit = glm(response[position] ~ x_part, family=binomial)
          new_data = model.matrix(ipw_fit)[relative_position, ]
          odds = exp(new_data %*% ipw_fit$coefficients)
          weight[position_0] = weight[position_0] + odds
        } else {
          ipw_fit = glm(response[position] ~ 1, family=binomial)
          odds = exp(ipw_fit$coefficients)
          weight[position_0] = weight[position_0] + odds
        }
      } else {
        # we observe at least one y
        x_part = as.matrix(cbind(data$x[position, observed_variables_x], 
                       data$y[position, observed_variables_y]))
        ipw_fit = glm(response[position] ~ x_part, family=binomial)
        new_data = model.matrix(ipw_fit)[relative_position, ]
        odds = exp(new_data %*% ipw_fit$coefficients)
        weight[position_0] = weight[position_0] + odds
      }
    }
  }
  
  response = data$y[data$A == cp_num, 3]
  covariate = as.matrix(data$y[data$A == cp_num, 1:2])
  result_lm = lm(response ~ covariate, weights=weight[data$A == cp_num])
  beta_estimate = result_lm$coefficients
  return(beta_estimate)
}


bootstrap_regression <- function(data, n_B) {
  n = length(data$A)
  beta_estimate = matrix(0, nrow=n_B, ncol=3)
  
  for (i in 1:n_B) {
    if (i %% 100 == 0) {
      cat("now is the ", i, " iteration\n")
    }
    bt_index = sample(n, size=n, replace=T)
    data_bt = list()
    data_bt$X = data$X[bt_index, ]
    data_bt$Y = data$Y[bt_index, ]
    data_bt$A = data$A[bt_index]
    data_bt$R = data$R[bt_index]
    beta_estimate[i, ] = ipw_marginal_parametric_model(data_bt)
  }
  return(beta_estimate)
  
}

regression_lm = ipw_marginal_parametric_model(mpm_data)
regression_bt = bootstrap_regression(mpm_data, n_B)
regression_ci = apply(regression_bt, 2, function(x) quantile(x, probs=c(0.025, 0.975)))


response = mpm_data$Y[mpm_data$A == 1, 3]
covariates = as.matrix(mpm_data$Y[mpm_data$A == 1, 1:2])
result_lm = lm(response ~ covariates)
cc_lm_ci = confint(result_lm)

ipw_marginal_parametric_model_sensitivity_analysis <- function(data, delta) {
  n = length(data$A)
  weight = rep(0, n)
  
  # estimate the mean when A = 1
  position_A_1 = data$A == 1
  weight[position_A_1] = 1
  
  for (i in 2:8) {
    missing_pattern_Y = mpm_possible_missing_pattern_for_Y[i, ]
    observed_variables_Y = which(missing_pattern_Y == T)
    missed_variables_Y = which(missing_pattern_Y == F)
    num_observed_variables_Y = length(observed_variables_Y)
    
    for (j in 1:4) {
      missing_pattern_X = mpm_possible_missing_pattern_for_X[j, ]
      observed_variables_X = which(missing_pattern_X == T)
      num_observed_variables_X = length(observed_variables_X)
      
      if (num_observed_variables_X == 0) {
        more_X_observed = 1:4
      } else if (num_observed_variables_X == 1) {
        more_X_observed = which(mpm_possible_missing_pattern_for_X[, observed_variables_X])
      } else {
        more_X_observed = which(apply(
          mpm_possible_missing_pattern_for_X[, observed_variables_X], 
          1, function(v) sum(v) == num_observed_variables_X))
      }
      
      position_0 = data$A == 1 & data$R %in% more_X_observed
      position_1 = data$A == i & data$R == j
      position = position_0 | position_1
      relative_position = data$A[position] == 1 & data$R[position] %in% more_X_observed
      response = rep(0, n)
      response[position_1] = 1
      
      if (num_observed_variables_Y == 0) {
        if (num_observed_variables_X > 0) {
          X_part = data$X[position, observed_variables_X]
          ipw_fit = glm(response[position] ~ X_part, family=binomial)
          new_data = model.matrix(ipw_fit)[relative_position, ]
          odds_depend_on_unobserved = exp((data$Y[position_0, 1] + data$Y[position_0, 2] + 
                                           data$Y[position_0, 3] - 21) * delta)
          odds = exp(new_data %*% ipw_fit$coefficients)
          odds_sensitivity = odds * odds_depend_on_unobserved
          # weight[position_0] = weight[position_0] + odds
          weight[position_0] = weight[position_0] + odds_sensitivity
        } else {
          ipw_fit = glm(response[position] ~ 1, family=binomial)
          odds = exp(ipw_fit$coefficients)
          odds_depend_on_unobserved = exp((data$Y[position_0, 1] + data$Y[position_0, 2] + 
                                           data$Y[position_0, 3] - 21) * delta)
          odds_sensitivity = odds * odds_depend_on_unobserved
          # weight[position_0] = weight[position_0] + odds
          weight[position_0] = weight[position_0] + odds_sensitivity
        }
      } else {
        # we observe one Y
        X_part = as.matrix(cbind(data$X[position, observed_variables_X], 
                       data$Y[position, observed_variables_Y]))
        ipw_fit = glm(response[position] ~ X_part, family=binomial)
        new_data = model.matrix(ipw_fit)[relative_position, ]
        odds = exp(new_data %*% ipw_fit$coefficients)
        if (num_observed_variables_Y == 2) {
          odds_depend_on_unobserved = exp((data$Y[position_0, -observed_variables_Y] - 7) * delta)
        } else {
          odds_depend_on_unobserved = exp(rowSums((data$Y[position_0, -observed_variables_Y] - 7) * delta))
        }
        odds_sensitivity = odds_depend_on_unobserved * odds
        # weight[position_0] = weight[position_0] + odds
        weight[position_0] = weight[position_0] + odds_sensitivity
      }
    }
  }
  
  response = data$Y[data$A == 1, 3]
  covariate = as.matrix(data$Y[data$A == 1, 1:2])
  result_lm = lm(response ~ covariate, weights=weight[data$A == 1])
  beta_estimate = result_lm$coefficients
  return(beta_estimate)
}

marginal_parametric_bootstrap_sensitivity_analysis <- function(data, f, n_B=1000, delta) {
  n = length(data$A)
  theta_estimate = matrix(0, nrow=n_B, ncol=3)
  
  for (i in 1:n_B) {
    if (i %% 100 == 0) {
      cat("now is the ", i, " iteration\n")
    }
    bt_index = sample(n, size=n, replace=T)
    data_bt = list()
    data_bt$X = data$X[bt_index, ]
    data_bt$Y = data$Y[bt_index, ]
    data_bt$A = data$A[bt_index]
    data_bt$R = data$R[bt_index]
    theta_estimate[i, ] = f(data_bt, delta)
  }
  return(theta_estimate)
}

delta = seq(-0.5, 0.5, 0.1)
n_B = 1000
marginal_parametric_ipw_sa = matrix(0, nrow=length(delta), ncol=3)
marginal_parametric_ipw_sa_sd = matrix(0, nrow=length(delta), ncol=3)
marginal_parametric_ipw_sa_l = matrix(0, nrow=length(delta), ncol=3)
marginal_parametric_ipw_sa_u = matrix(0, nrow=length(delta), ncol=3)

for (i in 1:length(delta)) {
  marginal_parametric_ipw_sa[i, ] = ipw_marginal_parametric_model_sensitivity_analysis(mpm_data, delta[i])
  marginal_parametric_ipw_sa_bt = marginal_parametric_bootstrap_sensitivity_analysis(mpm_data,
                                  ipw_marginal_parametric_model_sensitivity_analysis,
                                  n_B, delta[i]
  )
  betas_ci = apply(marginal_parametric_ipw_sa_bt, 2, function(v) quantile(v, probs=c(0.025, 0.975)))
  marginal_parametric_ipw_sa_sd[i, ] = apply(marginal_parametric_ipw_sa_bt, 2, sd)
  marginal_parametric_ipw_sa_l[i, ] = betas_ci[1, ]
  marginal_parametric_ipw_sa_u[i, ] = betas_ci[2, ]
}


delta_0 = c(delta, 0)
beta_0_sa_estimate = c(marginal_parametric_ipw_sa[, 1], result_lm$coefficients[1])
beta_0_sa_ci_l = c(marginal_parametric_ipw_sa_l[, 1], cc_lm_ci[1, 1])
beta_0_sa_ci_u = c(marginal_parametric_ipw_sa_u[, 1], cc_lm_ci[1, 2])
method = c(rep("IPW perturbed", length(delta)),  "complete-case")

beta_0_theta_df = data.frame(delta=delta_0, sa_estimate = beta_0_sa_estimate,
                      sa_ci_l = beta_0_sa_ci_l,
                      sa_ci_u = beta_0_sa_ci_u,
                      method=method)
library(ggplot2)
cbp2 <- c("red", "black")
g2 <- ggplot(beta_0_theta_df) + 
  geom_point(aes(x=delta, y=sa_estimate, color=method)) + 
  geom_errorbar(aes(ymax=sa_ci_u, ymin=sa_ci_l, x=delta, color=method), width=0.08) + 
  scale_colour_manual(values=cbp2) + 
  xlab(expression(delta)) + 
  ylab(expression(beta[0])) + 
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1)) + 
  ggtitle(expression(paste("Sensitivity analysis for ", beta[0], sep=""))) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust=0.5, size=18),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16, face="bold"),
    legend.position="none"
  )

plot(g2)

# delta_0 = c(delta, 0)
beta_1_sa_estimate = c(marginal_parametric_ipw_sa[, 2], result_lm$coefficients[2])
beta_1_sa_ci_l = c(marginal_parametric_ipw_sa_l[, 2], cc_lm_ci[2, 1])
beta_1_sa_ci_u = c(marginal_parametric_ipw_sa_u[, 2], cc_lm_ci[2, 2])
method = c(rep("IPW perturbed", length(delta)),  "complete-case")

beta_1_theta_df = data.frame(delta=delta_0, sa_estimate = beta_1_sa_estimate,
                      sa_ci_l = beta_1_sa_ci_l,
                      sa_ci_u = beta_1_sa_ci_u,
                      method=method)

g2 <- ggplot(beta_1_theta_df) + 
  geom_point(aes(x=delta, y=sa_estimate, color=method)) + 
  geom_errorbar(aes(ymax=sa_ci_u, ymin=sa_ci_l, x=delta, color=method), width=0.08) + 
  scale_colour_manual(values=cbp2) + 
  xlab(expression(delta)) + 
  ylab(expression(beta[1])) + 
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1)) + 
  ggtitle(expression(paste("Sensitivity analysis for ", beta[1], sep=""))) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust=0.5, size=18),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16, face="bold"),
    legend.position="none"
  )

plot(g2)

beta_2_sa_estimate = c(marginal_parametric_ipw_sa[, 3], result_lm$coefficients[3])
beta_2_sa_ci_l = c(marginal_parametric_ipw_sa_l[, 3], cc_lm_ci[3, 1])
beta_2_sa_ci_u = c(marginal_parametric_ipw_sa_u[, 3], cc_lm_ci[3, 2])
method = c(rep("IPW perturbed", length(delta)),  "complete-case")

beta_2_theta_df = data.frame(delta=delta_0, sa_estimate = beta_2_sa_estimate,
                      sa_ci_l = beta_2_sa_ci_l,
                      sa_ci_u = beta_2_sa_ci_u,
                      method=method)

g2 <- ggplot(beta_2_theta_df) + 
  geom_point(aes(x=delta, y=sa_estimate, color=method)) + 
  geom_errorbar(aes(ymax=sa_ci_u, ymin=sa_ci_l, x=delta, color=method), width=0.08) + 
  scale_colour_manual(values=cbp2) + 
  xlab(expression(delta)) + 
  ylab(expression(beta[1])) + 
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1)) + 
  ggtitle(expression(paste("Sensitivity analysis for ", beta[2], sep=""))) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust=0.5, size=18),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16, face="bold"),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    legend.position = c(0.2, 0.9)
  )

plot(g2)
