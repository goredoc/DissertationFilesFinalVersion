transform_to_normal_space = function(params){
  primes <- params
  primes <- log(params)
  # primes["b_speed"] <- -log(1 / params["b_speed"] - 1)
  # primes["v1_med"] <- -log(1 / params["v1_med"] - 1)
  # primes["v1_hard"] <- -log(1 / params["v1_hard"] - 1)
  primes["A"] <- -log(1 / params["A"] - 1)
  primes["rho"] <- -log((1 / params["rho"]) - 1)
  primes
}

transform_to_normal_space_probit = function(params){
  primes <- params
  primes <- log(params)
  # primes["b_speed"] <- -log(1 / params["b_speed"] - 1)
  # primes["v1_med"] <- -log(1 / params["v1_med"] - 1)
  # primes["v1_hard"] <- -log(1 / params["v1_hard"] - 1)
  primes["A"]   <- qnorm(params["A"])
  primes["rho"] <- qnorm(params["rho"])
  primes
}

transform_to_balba_space = function(primes){
  params <- primes
  params <- exp(primes)
  # params["b_speed"] <- 1 / (1 + exp(-primes["b_speed"]))
  # params["v1_med"] <- 1 / (1 + exp(-primes["v1_med"]))
  # params["v1_hard"] <- 1 / (1 + exp(-primes["v1_hard"]))
  params["A"] <- 1 / (1 + exp(-primes["A"]))
  params["rho"] <- 1 / (1 + exp(-primes["rho"]))
  params
}

transform_to_balba_space_probit = function(primes){
  params <- primes
  params <- exp(primes)
  # params["b_speed"] <- 1 / (1 + exp(-primes["b_speed"]))
  # params["v1_med"] <- 1 / (1 + exp(-primes["v1_med"]))
  # params["v1_hard"] <- 1 / (1 + exp(-primes["v1_hard"]))
  params["A"] <- pnorm(primes["A"])
  params["rho"] <- pnorm(primes["rho"])
  params
}


rbalba_handler <- function(n, v1, v2, A, alpha, beta, b) {
  t <- numeric(n)
  r <- numeric(n)
  dat <- NULL
  for(i in 1:n){
    out <- c(-1, -1)
    while(out[1] < 0){
      out <- rbalba(1, 
                    v1 = v1, 
                    v2 = v2, 
                    A = A, 
                    alpha = alpha[i], 
                    beta = beta[i], 
                    b = b, 
                    seed = get_seed())
    }
    dat <- rbind(dat, out)
  }
  # print(dat)
  return(list(rt = unname(dat[, 1]), resp = unname(dat[, 2])))
}

get_seed <- function() {
  sample.int(500, 1)
}


# -----------------------------------------------------
### pnum needs to be added
simulate_exp <- function(pnum,
                         n, 
                         num_blocks, 
                         cond_vec, 
                         b_acc, 
                         b_speed, 
                         v1_A, 
                         v1_B, 
                         v1_C, 
                         v2, 
                         rho, 
                         A){
  cond_vec_list <- list()
  t_list <- list()
  resp_list <- list()
  for(bl in 1:num_blocks){
    if(bl %% 2 == 1){
      b <- b_acc
    }else{
      b <- b_speed
    }
    print(bl)
    stim_vec <- (cond_vec > 3) + 1
    t <- resp <-  numeric(n) 
    
    alpha <- beta <- numeric(n)
    alpha[1] <- beta[1] <- 1
    for(i in 2:n){
      #      alpha[i] <- rho * alpha[i-1] + 2^(2 - stim_vec[i-1])
      #      beta[i] <- rho * beta[i-1] + 2^(stim_vec[i-1] - 1)
      alpha[i] <- 1
      beta[i]  <- 1
      
    }
    
    print("alpha")
    for(cond_i in 1:6){
      print(c("cond ", cond_i))
      # rights are correct
      if(cond_i == 1){
        v1 <- v1_A
        logic_tmp <- cond_vec == cond_i
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        out <- rbalba_handler(length(logic_tmp), v1, v2, A, alpha_vec, beta_vec, b)
        resp[logic_tmp] <- out$resp
        t[logic_tmp] <- out$rt
        
      }else if(cond_i == 2){
        v1 <- v1_B
        logic_tmp <- cond_vec == cond_i
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        out <- rbalba_handler(length(logic_tmp), v1, v2, A, alpha_vec, beta_vec, b)
        resp[logic_tmp] <- out$resp
        t[logic_tmp] <- out$rt
        
      }else if(cond_i == 3){
        v1 <- v1_C
        logic_tmp <- cond_vec == cond_i
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        out <- rbalba_handler(length(logic_tmp), v1, v2, A, alpha_vec, beta_vec, b)
        resp[logic_tmp] <- out$resp
        t[logic_tmp] <- out$rt
      }else if(cond_i == 4){ # left is correct
        v1 <- v1_A
        logic_tmp <- cond_vec == cond_i
        alpha_vec <- beta[logic_tmp]
        beta_vec <- alpha[logic_tmp]
        out <- rbalba_handler(length(logic_tmp), v1, v2, A, alpha_vec, beta_vec, b)
        resp[logic_tmp] <- out$resp
        t[logic_tmp] <- out$rt
      }else if(cond_i == 5){
        v1 <- v1_B
        logic_tmp <- cond_vec == cond_i
        alpha_vec <- beta[logic_tmp]
        beta_vec <- alpha[logic_tmp]
        out <- rbalba_handler(length(logic_tmp), v1, v2, A, alpha_vec, beta_vec, b)
        resp[logic_tmp] <- out$resp
        t[logic_tmp] <- out$rt
      }else if(cond_i == 6){
        v1 <- v1_C
        logic_tmp <- cond_vec == cond_i
        alpha_vec <- beta[logic_tmp]
        beta_vec <- alpha[logic_tmp]
        out <- rbalba_handler(length(logic_tmp), v1, v2, A, alpha_vec, beta_vec, b)
        resp[logic_tmp] <- out$resp
        t[logic_tmp] <- out$rt
      }
      
    }
    # print(summary(t))
    # cond_vec_list[[bl]] <- cond_vec
    t_list[[bl]] <- t
    resp_list[[bl]] <- resp
  }
  list(t_list = t_list, resp_list = resp_list)
}


log_dens_like <- function(params, 
                          data, 
                          param_names, 
                          params_fixed, 
                          param_names_fixed,
                          mc){
  dens <- 0
  params <- transform_to_balba_space(params)
  names(params) <- param_names
  params[param_names_fixed] <- params_fixed
  
  for(bl in 1:num_blocks){
    if(bl %% 2 ==1 ){
      params["b"] <- params["b_acc"]
    }else{
      
      params["b"] <- params["b_speed"]
    }
    
    t <- data$t_list[[bl]]
    resp <- data$resp_list[[bl]]
    
    stim_vec <- (cond_vec > 3) + 1 
    alpha <- beta <- length(cond_vec)
    alpha[1] <- beta[1] <- 1
    for(i in 2:n){
      # alpha[i] <- params["rho"] * alpha[i-1] + 2^(2 - stim_vec[i-1])
      # beta[i] <- params["rho"] * beta[i-1]  + 2^(stim_vec[i-1] - 1)
      alpha[i] <- 1
      beta[i]  <- 1 
    }
    for(cond_i in 1:6){
      
      # rights are correct
      if(cond_i == 1){
        logic_tmp <- cond_vec == cond_i
        params["v1"] <- params["v1_A"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, length(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = alpha_vec, 
                                     beta = beta_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
        
      }else if(cond_i == 2){
        logic_tmp <- cond_vec == cond_i
        params["v1"] <- params["v1_B"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, length(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = alpha_vec, 
                                     beta = beta_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
        
      }else if(cond_i == 3){
        logic_tmp <- cond_vec == cond_i
        params["v1"] <- params["v1_C"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, length(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = alpha_vec, 
                                     beta = beta_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }else if(cond_i == 4){ # left is correct
        logic_tmp <- cond_vec == cond_i
        params["v1"] <- params["v1_A"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, length(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = beta_vec, 
                                     beta = alpha_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }else if(cond_i == 5){
        logic_tmp <- cond_vec == cond_i
        params["v1"] <- params["v1_B"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, length(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = beta_vec, 
                                     beta = alpha_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }else if(cond_i == 6){
        logic_tmp <- cond_vec == cond_i
        params["v1"] <- params["v1_C"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, length(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = beta_vec, 
                                     beta = alpha_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }
      
    }
  }
  dens
}

log_dens_prior <- function(params, hyper_list) {
  out <- 0
  for (p in param_names) {
    tmp <- grep((paste(p, "_", sep = "")), phi_names, value = TRUE)
    out <- out + dnorm(params[p], 
                       hyper_list[[tmp[1]]], 
                       hyper_list[[tmp[2]]],  
                       log = TRUE
    )
  }
  out
}

#--------------------------------------------------- Hyper functions
#' Calulates the log density of hyper distribution likelihood and prior
#' 
#' @param theta - vector of theta estimates for individuals
#' @param phi - vector of phi estimate, mean and standard deviation
#' @param p - name or index of current parameter
#' 
#' @example  
#' theta <- start_points[1]
#' phi <- c(start_points[1], start_sds[1])
#' prior <- prior_list[[1]]
#' p <- "v"
#' 
#' log_dens_hyper(theta, phi, prior, p)
log_dens_hyper <- function(theta, phi, prior, p) {
  sum(dnorm(theta, phi[1], phi[2], log = TRUE) +
        dnorm(phi[1], prior$mu[1], prior$mu[2], log = TRUE) +
        dgamma(phi[2], prior$sigma[1], prior$sigma[2], log = TRUE)
  )
}


log_dens_like_jones <- function(params, 
                          data, 
                          auxdata,
                          pnum,
                          param_names, 
                          params_fixed, 
                          param_names_fixed,
                          mc){
  dens <- 0
  params <- transform_to_balba_space(params)
  names(params) <- param_names
  params[param_names_fixed] <- params_fixed
  
  for(bl in 1:num_blocks){
    #if(bl %% 2 ==1 ){
    #  params["b"] <- params["b_acc"]
    #}else{
      
    #  params["b"] <- params["b_speed"]
    #}

    sac = auxdata$speedaccuracyYN[[bl]][1]
    if (sac==1) {params["b"] <- params["b_acc"]}
    else {params["b"] <- params["b_acc"]} #?
    # T0 = 100
    # t <- data$t_list[[bl]] - 100
    t <- data$t_list[[bl]]
    resp <- data$resp_list[[bl]]
    censorYN <- auxdata$censorYN[[bl]]
    
    cond_vec <- auxdata$freqCondition[[bl]]
    
    stim_vec <- (cond_vec > 3) + 1
    alpha <- beta <- length(cond_vec)
    alpha[1] <- beta[1] <- 1
    for(i in 2:n){
      # alpha[i] <- params["rho"] * alpha[i-1] + 2^(2 - stim_vec[i-1])
      # beta[i]  <- params["rho"] * beta[i-1] + 2^(stim_vec[i-1] - 1)
       alpha[i] <- 1
       beta[i]  <- 1
    }
    for(cond_i in 1:6){
      
      # rights are correct
      if(cond_i == 1){
        #logic_tmp <- cond_vec == cond_i
        logic_tmp <- (cond_vec == cond_i) & (censorYN == 0)
        params["v1"] <- params["v1_A"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, sum(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = alpha_vec, 
                                     beta = beta_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
        
      }else if(cond_i == 2){
        #logic_tmp <- cond_vec == cond_i
        logic_tmp <- (cond_vec == cond_i) & (censorYN == 0)
        params["v1"] <- params["v1_B"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, sum(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = alpha_vec, 
                                     beta = beta_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
        
      }else if(cond_i == 3){
        #logic_tmp <- cond_vec == cond_i
        logic_tmp <- (cond_vec == cond_i) & (censorYN == 0)
        params["v1"] <- params["v1_C"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, sum(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = alpha_vec, 
                                     beta = beta_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }else if(cond_i == 4){ # left is correct
        #logic_tmp <- cond_vec == cond_i
        logic_tmp <- (cond_vec == cond_i) & (censorYN == 0)
        params["v1"] <- params["v1_A"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, sum(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = beta_vec, 
                                     beta = alpha_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }else if(cond_i == 5){
        #logic_tmp <- cond_vec == cond_i
        logic_tmp <- (cond_vec == cond_i) & (censorYN == 0)
        params["v1"] <- params["v1_B"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, sum(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = beta_vec, 
                                     beta = alpha_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }else if(cond_i == 6){
        #logic_tmp <- cond_vec == cond_i
        logic_tmp <- (cond_vec == cond_i) & (censorYN == 0)
        params["v1"] <- params["v1_C"]
        alpha_vec <- alpha[logic_tmp]
        beta_vec <- beta[logic_tmp]
        dens <- dens + logdens_balba(rt = t[logic_tmp], 
                                     resp = resp[logic_tmp], 
                                     stim = rep(1, sum(logic_tmp)),
                                     v1 = params["v1"],
                                     v2 = params["v2"],
                                     A = params["A"], 
                                     alpha = beta_vec, 
                                     beta = alpha_vec, 
                                     b = params["b"], 
                                     mc_size = mc, 
                                     seed = get_seed())
      }
      
    }
  }
  dens
}

