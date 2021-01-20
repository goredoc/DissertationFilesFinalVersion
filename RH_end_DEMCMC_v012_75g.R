# Saved as src/DEMCMC_functions.R

# Functions used to perform DEMCMC posterior estimation

#' Computes Metropolis-Hastings Step
#' 
#' Determines whether to accept or reject the proposal based on weights of new and old proposals
#' 
#' @param new_weight - new weight of proposal to determine whether to accept
#' @param old_weight - old weight of former proposal value
#' 
#' @examples 
#' to_test <- c(-Inf, Inf, NULL, NA, NaN, 0, 1, -1, 99, -99)
#' to_test_pairs <- expand.grid(to_test, to_test)
#' 
#' for(p in 1:nrow(to_test_pairs)){
#'   test_pair <- unname(unlist(to_test_pairs[p, ]))
#'   print("Pair: ")
#'   print(test_pair)
#'   print(mh_step(test_pair[1], test_pair[2]))
#' }
mh_step <- function(new_weight, old_weight) {
  out <- "reject"
  
  if (is.null(new_weight)) new_weight <- -Inf
  if (!is.na(new_weight) & !is.finite(new_weight) & new_weight < 0) new_weight <- -Inf
  if (is.na(new_weight)) new_weight <- -Inf
  
  if (is.null(old_weight)) old_weight <- -Inf
  if (!is.na(old_weight) & !is.finite(old_weight) & old_weight < 0) old_weight <- -Inf
  if (is.na(old_weight)) old_weight <- -Inf
  
  
  mh <- exp(new_weight - old_weight)
  if (is.finite(mh)) {
    if (runif(1) < mh) {
      out <- "accept"
    }
  }else if (is.infinite(mh) & mh > 1){
    out <- "accept"
  }
  out
}


# ----------------------------------------------------------------------- Individual distribution

#' Calculates crossover step in DEMCMC
#' 
#' @param i - chain to be updated
#' @param param_inds - indices of parameters to perform crossover step
#' @param use_theta - matrix of theta to update
#' @param use_like - weight for theta
#' @param data - one subject's list of data: Resp, Time, Cond, Stim
#' @param hyper - matrix of hyper parameters
#' @param param_names - vector of names of indivudal level parameters to estimate
#' @param curr_iter - current iteration of samples, used to redraw from likelihood
#' 
#' @example 
#' See tests/recover_conjugate.R
crossover <- function(i, 
                      param_inds, 
                      use_theta, 
                      use_like, 
                      data, 
                      hyper, 
                      param_names,
                      params_fixed, 
                      param_names_fixed,
                      mc,
                      curr_iter,
                      CURRENT_TO_BEST) {
  if (curr_iter %% 5 == 0) use_like[i] <- log_dens_like(params = use_theta[, i], 
                                                        data = data, 
                                                        param_names = param_names, 
                                                        params_fixed = params_fixed, 
                                                        param_names_fixed = param_names_fixed,
                                                        mc = mc
                                                        )
  
  use_weight <- use_like[i] + log_dens_prior(use_theta[, i], hyper[, i])
  
  gamma <- .5
  index <- sample(c(1:num_chains)[-i], 2, replace = F)
  theta <- use_theta[, i]
  theta[param_inds] <- use_theta[param_inds, i] + 
    gamma * (use_theta[param_inds, index[1]] - use_theta[param_inds, index[2]]) + 
    runif(1, -b, b)
  
  if(CURRENT_TO_BEST){
    max_ind <- which.max(use_like)
    theta[param_inds] <- use_theta[param_inds, i] + 
      gamma * (- use_theta[param_inds, i] + use_theta[param_inds, max_ind]) + 
      runif(1, -b, b)
  }
  
  # print(theta)
  
  prior_like <- log_dens_prior(theta, hyper[, i])
  if (is.finite(prior_like) & prior_like > -Inf) {
    like <- try(log_dens_like(params = theta, 
                          data = data, 
                          param_names = param_names, 
                          params_fixed = params_fixed, 
                          param_names_fixed = param_names_fixed,
                          mc = mc
                          ))
    like <- ifelse(is.character(like), -Inf, like)
  } else {
    like <- -Inf
  }
  weight <- like + prior_like
  
  mh <- mh_step(new_weight = weight, old_weight = use_weight)
  if(mh == "accept"){
    use_theta[, i] <- theta
    use_like[i] <-  weight
  }
  c(use_like[i], use_theta[, i])
}

crossover_jones <- function(i, 
                      param_inds, 
                      use_theta, 
                      use_like, 
                      data, 
                      auxdata,
                      hyper, 
                      param_names,
                      params_fixed, 
                      param_names_fixed,
                      mc,
                      curr_iter,
                      CURRENT_TO_BEST) {
  if (curr_iter %% 5 == 0) use_like[i] <- log_dens_like_jones(params = use_theta[, i], 
                                                        data = data,
                                                        auxdata = auxdata,
                                                        param_names = param_names, 
                                                        params_fixed = params_fixed, 
                                                        param_names_fixed = param_names_fixed,
                                                        mc = mc
  )
  
  use_weight <- use_like[i] + log_dens_prior(use_theta[, i], hyper[, i])
  
  gamma <- .5
  index <- sample(c(1:num_chains)[-i], 2, replace = F)
  theta <- use_theta[, i]
  theta[param_inds] <- use_theta[param_inds, i] + 
    gamma * (use_theta[param_inds, index[1]] - use_theta[param_inds, index[2]]) + 
    runif(1, -b, b)
  
  if(CURRENT_TO_BEST){
    max_ind <- which.max(use_like)
    theta[param_inds] <- use_theta[param_inds, i] + 
      gamma * (- use_theta[param_inds, i] + use_theta[param_inds, max_ind]) + 
      runif(1, -b, b)
  }
  
  # print(theta)
  
  prior_like <- log_dens_prior(theta, hyper[, i])
  if (is.finite(prior_like) & prior_like > -Inf) {
    like <- try(log_dens_like_jones(params = theta, 
                              data = data,
                              auxdata = auxdata,
                              param_names = param_names, 
                              params_fixed = params_fixed, 
                              param_names_fixed = param_names_fixed,
                              mc = mc
    )) # like is log prob with no information about group
    # store like for each individual to compute WAIS
    like <- ifelse(is.character(like), -Inf, like)
  } else {
    like <- -Inf
  }
  weight <- like + prior_like # prior at group plus like individual
  
  mh <- mh_step(new_weight = weight, old_weight = use_weight)
  if(mh == "accept"){
    use_theta[, i] <- theta
    use_like[i] <-  weight
  }
  c(use_like[i], use_theta[, i]) # use_like is like + prior = weight
  # Ill want to put in just likelihood
  # this comes out then goes to hierarchical fit script

}


#' Performs migration step of DEMCMC
#' 
#' @param params_inds - indices of parameters to estimate
#' @param use_theta - matrix of theta estimates chains x parameters
#' @param use_like - vector of current likelihood evaluations
#' @param data - one subject's list of data: Resp, Time, Cond, Stim
#' @param hyper - matrix of hyper distribution
#' @param param_names - vector of names of parameters to be freely estimated
#' 
#' @example 
#' See tests/recover_conjugate.R
migration_crossover <- function(param_inds, use_theta, use_like, data, hyper, param_names) {
  n_migration_chains <- ceiling(runif(1, 0, num_chains))
  use_chains <- sample(1:num_chains, n_migration_chains)
  migration_use_weight <- rep(NA, n_migration_chains)
  migration_weight <- rep(NA, n_migration_chains)
  
  for (mi in 1:n_migration_chains) {
    migration_use_weight[mi] <- use_like[use_chains[mi]] + 
      log_dens_prior(use_theta[param_inds, use_chains[mi]], hyper[, use_chains[mi]])
    new_chain <- mi - 1
    
    if (new_chain == 0) new_chain <- n_migration_chains
    
    migration_weight[mi] <- use_like[use_chains[new_chain]] + 
      log_dens_prior(use_theta[param_inds, use_chains[new_chain]], hyper[, use_chains[new_chain]])
    
    
    mh <- mh_step(new_weight = migration_weight[mi], old_weight = migration_use_weight[mi])
    if(mh == "accept"){
      use_theta[, use_chains[mi]] <- use_theta[, use_chains[new_chain]]
      use_like[use_chains[mi]] <- use_like[use_chains[new_chain]]
    }
  }
  rbind(use_like, use_theta)
}

# ----------------------------------------------------------------------- Hyper distribution

#' Calculates crossover step in DEMCMC for hyper distribution
#' 
#' @param i - chain to be updated
#' @param param_inds - indices of parameters to perform crossover step
#' @param use_theta - matrix of theta
#' @param use_like - weight for theta
#' @param use_phi - matrix of phi to be updated
#' @param prior - list of priors for this parameters mean and standard deviation hyperdistributions
#' @param p - current index or name of parameter to update
#' 
#' @example 
#' See tests/recover_conjugate.R
crossover_hyper <- function(i, 
                            param_inds, 
                            use_theta, 
                            use_phi, 
                            use_weight, 
                            prior, 
                            p, 
                            CURRENT_TO_BEST) {
  use_weight_i <- log_dens_hyper(use_theta[, i], use_phi[param_inds, i], prior, p)
  gamma <- .5
  index <- sample(c(1:num_chains)[-i], 2, replace = F)
  phi <- use_phi[, i]
  
  phi[param_inds] <- use_phi[param_inds, i] + 
    gamma * (use_phi[param_inds, index[1]] - use_phi[param_inds, index[2]]) + 
    runif(1, -b, b)
  
  if(CURRENT_TO_BEST){
    max_ind <- which.max(use_weight)
    phi[param_inds] <- use_phi[param_inds, i] + 
      gamma * (use_phi[param_inds, max_ind] - use_phi[param_inds, i]) + 
      runif(1, -b, b)
  }
  
  weight <- log_dens_hyper(use_theta[, i], phi[param_inds], prior, p)
  
  mh <- mh_step(new_weight = weight, old_weight = use_weight_i)
  if(mh == "accept"){
    use_phi[, i] <- phi
    use_weight[i] <- weight
  } else {
    use_weight[i] <- use_weight_i
  }
  c(use_weight[i], use_phi[, i])
}

crossover_hyper_one_sub <- function(i, param_inds, use_theta, use_phi, prior, p) {
  # use_weight <- log_dens_hyper(use_theta[, i], use_phi[param_inds, i], prior, p)
  use_weight <- log_dens_hyper(use_theta[i], use_phi[param_inds, i], prior, p)
  gamma <- .5
  index <- sample(c(1:num_chains)[-i], 2, replace = F)
  phi <- use_phi[, i]
  
  phi[param_inds] <- use_phi[param_inds, i] + 
    gamma * (use_phi[param_inds, index[1]] - use_phi[param_inds, index[2]]) + 
    runif(1, -b, b)
  
  # weight <- log_dens_hyper(use_theta[, i], phi[param_inds], prior, p)
  weight <- log_dens_hyper(use_theta[i], phi[param_inds], prior, p)
  
  mh <- mh_step(new_weight = weight, old_weight = use_weight)
  if(mh == "accept"){
    use_phi[, i] <- phi
  }
  use_phi[, i]
}
#' Performs migration step of DEMCMC for hyperdistribution
#' 
#' @param params_inds - indices of parameters to estimate
#' @param use_theta - matrix of theta estimates chains x parameters
#' @param use_phi - vector of current likelihood evaluations
#' @param use_phi - matrix of phi to be updated
#' @param prior - list of priors for this parameters mean and standard deviation hyperdistributions
#' @param p - current index or name of parameter to update
#' 
#' @example 
#' See tests/recover_conjugate.R
migration_crossover_hyper <- function(param_inds, use_theta, use_phi, prior, p) {
  n_migration_chains <- ceiling(runif(1, 0, num_chains))
  use_chains <- sample(1:num_chains, n_migration_chains)
  migration_use_weight <- rep(NA, n_migration_chains)
  migration_weight <- rep(NA, n_migration_chains)
  
  for (mi in 1:n_migration_chains) {
    migration_use_weight[mi] <- log_dens_hyper(use_theta[, use_chains[mi]], 
                                               use_phi[param_inds, use_chains[mi]], 
                                               prior, 
                                               p = p)
    new_chain <- mi - 1
    if (new_chain == 0) new_chain <- n_migration_chains
    migration_weight[mi] <- log_dens_hyper(use_theta[, use_chains[mi]], 
                                           use_phi[param_inds, use_chains[new_chain]], 
                                           prior, 
                                           p = p)
    
    mh <- mh_step(new_weight = migration_weight[mi], old_weight = migration_use_weight[mi])
    if (mh == "accept") {
      use_phi[param_inds, use_chains[mi]] <- use_phi[param_inds, use_chains[new_chain]]
    }
  }
  use_phi
}

# ---------------------------------------------------------------- Helpers

#' Combines matrices into arrays along the third dim
#' 
#' @param x - matrix to be combined
#' @param ... - other matrices to be combined
#' 
#' @example
#' bind_3(cars, cars)
bind_3 <- function(x, ...){
  abind::abind(x, ..., along = 3)
}

#' lengthens containers for DEMCMC
#' 
#' @param container - containers for weights, theta, phi
#' @param modifying_dim - dimension to extend
#' @param modifying_dim_length - new length of modified dimension
#' 
#' @example 
#' (not run)
lengthen_container <- function(container, modifying_dim, modifying_dim_length) {
  container_old <<- container
  new_dim <- dim(container)
  new_dim[modifying_dim] <- modifying_dim_length
  container <- array(NA, dim = new_dim)
  container
}