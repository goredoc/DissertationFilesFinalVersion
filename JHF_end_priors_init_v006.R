# Saved as Priors.R

# Creates priors, start points, boundaries and fixes parameters



#' Updates lists for prior, upper and lower bounds
#'
#' Appends new elements to prior, upper, and lower bounds lists
#'
#' @param param_name - parameter name to append to the lists
#' @param mu_mean - mean of mu hyperlevel prior
#' @param mu_sd - standard deviation of mu hyperlevel prior
#' @param sigma_shape - gamma shape of sigma hyperlevel prior
#' @param sigma_rate - gamma rate of sigma hyperlevel prior
#' @param lb - lower bound of the parameter to be added
#' @param ub - upper bound of the parameter to be added
#' @param param_names - names of all parameters to estimate
#' @param prior_list - list of priors to be appended
#' @param lower_bound_list - list of lower bounds to be appended
#' @param upper_bound_list - list of upper bounds to be appended
#'
#' @example
#'
#' prior_list <- list()
#' upper_bound_list <- list()
#' lower_bound_list <- list()
#'
#' out <- update_prior_bounds_lists(param_name = "v",
#'                                  mu_mean = .1, mu_sd = 1,
#'                                  sigma_shape = 1, sigma_rate = 1,
#'                                  lb = -Inf, ub = Inf,
#'                                  param_names = param_names,
#'                                  prior_list = prior_list,
#'                                  lower_bound_list = lower_bound_list,
#'                                  upper_bound_list = upper_bound_list)
#' out$prior_list
#' out$upper_bound_list
#' out$lower_bound_list
update_prior_bounds_lists <- function(param_name,
                                      mu_mean,
                                      mu_sd,
                                      sigma_shape,
                                      sigma_rate,
                                      lb,
                                      ub,
                                      param_names,
                                      prior_list,
                                      lower_bound_list,
                                      upper_bound_list) {
  tmp <- grep(param_name, param_names, value = TRUE)
  for (n in 1:length(tmp)) {
    tmp2 <- tmp[n]
    prior_list[[tmp2]] <- list(mu = c(mu_mean, mu_sd), sigma = c(sigma_shape, sigma_rate))
    lower_bound_list[[tmp2]] <- lb
    upper_bound_list[[tmp2]] <- ub
  }
  list(prior_list = prior_list, lower_bound_list = lower_bound_list, upper_bound_list = upper_bound_list)
}


# -------------------------------------------------------- Parameters
A <- 1
v1_A <- 4
v1_B <- 3
v1_C <- 2
v2 <- 1.5
b_acc <- 3
b_speed <- 1
rho <- .5

param_names <- c("v1_A", "v1_B", "v1_C", "v2", "b_acc", "b_speed", "rho")#
param_names_fixed <- c("A")
params_fixed <- c(A)

# -------------------------------------------------------- Priors

prior_list <- list()
upper_bound_list <- list()
lower_bound_list <- list()

#change line below to use the transformation function
start_points <- mu_mean_vec <- c(log(v1_A), log(v1_B), log(v1_C), log(v2), log(b_acc), log(b_speed), log(rho/(1-rho)))
start_sds <- mu_sd_vec <-      c(        1,         1,         1,       1,          1,            1,              1)
sigma_shape_vec        <-      c(   .574e5,    .574e5,    .574e5,  .574e5,    2.692e5,      2.692e5,        0.324e5) 
sigma_rate_vec         <- rep(100000,length(param_names))
#sigma_shape_vec  <- rep(1, length(param_names))
#sigma_rate_vec   <- rep(5, length(param_names))
lower_bounds_vec <- rep(-Inf, length(param_names))
upper_bounds_vec <- rep(Inf, length(param_names))

for (i in seq_along(param_names)) {
  out <- update_prior_bounds_lists(
    param_name = param_names[i],
    mu_mean = mu_mean_vec[i],
    mu_sd = mu_sd_vec[i],
    sigma_shape = sigma_shape_vec[i],
    sigma_rate = sigma_rate_vec[i],
    lb = lower_bounds_vec[i],
    ub = upper_bounds_vec[i],
    param_names = param_names,
    prior_list = prior_list,
    lower_bound_list = lower_bound_list,
    upper_bound_list = upper_bound_list
  )
  prior_list <- out$prior_list
  upper_bound_list <- out$upper_bound_list
  lower_bound_list <- out$lower_bound_list
}
