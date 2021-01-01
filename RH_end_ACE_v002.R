# purpose of RH_end_ACE_v002.R 
# Due to poor fit return to the old alpha, beta update rules

# purpose of RH_end_ACE_v001.R 
# make sure that revised code to recover parameters reasonably well
# changed updating rule for alpha to
# alpha[i] = rho*alpha[i-1] + 1 (if last stim is word)
# this did not fit the rho parameter well at all

stopCluster(cl)
rm(list = ls())
set.seed(1)
# Saved as recovery_hierarchical.R

# Recovery of Jones experiment

# import libraries
# setwd("C:/Users/Chief/Desktop/Re__BALBA_Rho_Recovery")
setwd("C:/Users/Robert Gore/Desktop/dissertation/20200930")

library(Rcpp)
library(tidyverse)

sourceCpp("BALBA.cpp")
source("RH_end_BALBA2_v002.R")
source("RH_end_priors_init_v002.R")

##########################################################

hypers_means <- NULL
for (p in seq_along(prior_list)) {
  hypers_means <- c(hypers_means, prior_list[[p]][[1]][1], prior_list[[p]][[1]][2])
}
names(hypers_means) <- rep(param_names, each = 2)

N <- 5
theta_prior_list <- list()
for (p in param_names) {
  theta_prior <- numeric(N)
  for (i in 1:N) {
    # Draw Mu
    Mu <- rnorm(1, mean = prior_list[[p]]$mu[1], sd = prior_list[[p]]$mu[2])
    # Draw Sigma
    Sigma <- rgamma(1, prior_list[[p]]$sigma[1], prior_list[[p]]$sigma[2])
    # Draw theta
    theta_prior[i] <- rnorm(1, Mu, Sigma)
  }
  theta_prior_list[[p]] <- theta_prior
}

# params_true <- params_true[1, ]
params_true <- theta_prior_list %>%
  unlist() %>%
  matrix(ncol = length(param_names), byrow = F)
colnames(params_true) <- param_names

params_trans <- params_true
for (s in 1:nrow(params_true)) {
  params_trans[s, ] <- transform_to_balba_space(params_true[s, ])[-c(8)]
}
params_trans

n <- 96
num_blocks <- 20

cond_vec <- sample(rep(1:6, 96/6), size = 96, replace = F)
Synth_data <- list()
for(s in 1:N){
  print(s)
  #debugonce(simulate_exp)
  synth_data <- simulate_exp(1,
                             n, 
                             num_blocks, 
                             cond_vec, 
                             b_acc   = params_trans[s, "b_acc"],
                             b_speed = params_trans[s, "b_speed"],
                             v1_A    = params_trans[s, "v1_A"],
                             v1_B    = params_trans[s, "v1_B"],
                             v1_C    = params_trans[s, "v1_C"], 
                             v2      = params_trans[s, "v2"], 
                             rho     = params_trans[s, "rho"],
                             A       = 1)
  print(s)
  Synth_data[[s]] <- synth_data
}


rm(list = setdiff(ls(), list("params_true",
                             "params_trans",
                             "Synth_data",
                             "cond_vec"
)
)
)

# -------------------------------------------------------- Set up Estimation
# import libraries
#setwd("C:/Users/Chief/Desktop/Re__BALBA_Rho_Recovery")
setwd("C:/Users/Robert Gore/Desktop/dissertation/20200930")

library(Rcpp)
library(tidyverse)
library(foreach)
library(doParallel)

sourceCpp("BALBA.cpp")
source("RH_end_BALBA2_v002.R")
source("RH_end_priors_reduced_v002.R")
source("RH_end_DEMCMC_v002.R")

Data <- Synth_data
S <- length(Data)

num_blocks <- 20
n <- 96

# Parameters for DEMCMC
num_pars <- length(param_names)
num_chains <- length(param_names) * 3
num_monte_carlo <- 1000
# number of monte carlo draws for the posterior - particles

CTB_vec <- seq(1, 25, by = 1)
b <- .001

# Kernel size for MC_Size
# set to 5 for initial troubleshooting
# 500 is a lot
# 50 may be good
# this is the number of draws from the beta start point
# this is the monte carlo integration across start points 
num_samples <- 5

# -------------------------------------------------------- Storage
theta <- array(NA, c(num_pars, S, num_chains, num_monte_carlo))
rownames(theta) <- param_names

phi <- array(NA, c(num_pars * 2, num_chains, num_monte_carlo))
phi_names <- paste(rep(param_names, each = 2), c("mu", "sigma"), sep = "_")
rownames(phi) <- phi_names

weight <- array(NA, c(S, num_chains, num_monte_carlo))
hyper_weight <- array(NA, c(num_chains, num_monte_carlo))

stopImplicitCluster()
registerDoParallel(3)
# -------------------------------------------------------- Initialization
for (i in 1:num_chains) {
  print(i)
  #i <- 1 #recommented with Noah 10/5/2020
  # s <- 1 # rm
  out <- foreach(s = 1:S, .combine = cbind, .packages = c("Rcpp"), 
                 .noexport = c("logdens_balba")) %dopar% {
    # s <- 1
    current_weight <- -Inf
    while (current_weight == -Inf) {
      sourceCpp("BALBA.cpp")
      current_thetas <- rnorm(n = num_pars, mean = start_points, sd = start_sds)
      current_weight <- log_dens_like(params = current_thetas, 
                                      data = Data[[s]], 
                                      param_names = param_names, 
                                      params_fixed = params_fixed, 
                                      param_names_fixed = param_names_fixed, 
                                      mc = num_samples)
    }
    c(current_weight, current_thetas)
  }
  weight[, i, 1] <- out[1, ]
  theta[, , i, 1] <- out[1+1:num_pars, ]
  # weight[, i, 1] <- out[1]
  # theta[, , i, 1] <- out[1+1:num_pars]
}

(theta[1,,, 1])

for (c in 1:num_chains) {
  hyper_weight[c, 1] <- -Inf
  while (!is.finite(hyper_weight[c, 1])) {
    for(p in param_names){
      tmp <- paste(p, "sigma", sep = "_")
      phi[tmp, c, 1] <- rgamma(1, prior_list[[p]]$sigma[1], prior_list[[p]]$sigma[2])
      
      tmp <- (paste(p, "mu", sep = "_"))
      phi[tmp, c, 1] <- rnorm(1, prior_list[[p]]$mu[1], prior_list[[p]]$mu[2])
    }
    
    new_weight <- 0
    for (p in param_names) {
      print(new_weight)
      print(p)
      which_theta <- match(x = p, table = param_names)
      which_phi <- match(x = paste(p, c("mu", "sigma"), sep = "_"), table = phi_names)
      print(phi[which_phi, c, 1])
      new_weight <- new_weight + log_dens_hyper(theta = theta[which_theta,, c, 1], 
                                                phi = phi[which_phi, c, 1], 
                                                prior = prior_list[[p]], 
                                                p = p)
    }
    hyper_weight[c, 1] <- new_weight
  }
}
print("Initialized")

I <- 3
stopImplicitCluster()
registerDoParallel(3)
# -------------------------------------------------------- Sample
print("Sampling")
tic <- tictoc::tic()
for (i in (I - 1):num_monte_carlo) {
  # i <- 2
  if(i %in% CTB_vec){
    CURRENT_TO_BEST <- T
  } else {
    CURRENT_TO_BEST <- F
  }
  
  cat("\n ", i, "  ")
  if (i %% 100 == 0) save.image("RH_end_v002_centsave.RData")
  # Hyper level
  phi[,, i] <- phi[,, i - 1]
  hyper_weight[, i] <- hyper_weight[, i - 1]
  for (p in param_names) {
    which_theta <- match(x = p, table = param_names)
    which_phi <- match(x = paste(p, c("mu", "sigma"), sep = "_"), table = phi_names)
    # if (i %in% CURRENT_TO_BEST) {
    # phi[,,i] = migration_crossover_hyper(param_inds = which_phi, 
    #                                      use_theta = theta[which_theta,,, i - 1], 
    #                                      use_phi = phi[,, i], 
    #                                      prior = prior_list[[p]],
    #                                      p = p,
    #                                      )
    # } else {
    # debugonce(crossover_hyper)
    out <- sapply(1:num_chains,
                  crossover_hyper,
                  param_inds = which_phi,
                  use_theta = theta[which_theta,,, i - 1],
                  use_phi = phi[,, i],
                  use_weight = hyper_weight[, i],
                  prior = prior_list[[p]],
                  p = p, 
                  CURRENT_TO_BEST = F
    )
    hyper_weight[, i] <- out[1, ]
    phi[,, i] <- out[1+1:(num_pars * 2), ]
  }
  # }
  
  # Individual level
  # if (i %in% CURRENT_TO_BEST) {
  #   temp = foreach(s = 1:S) %dopar% {
  #     library(here)
  #     library(msm)
  #     if(!Linux)dyn.load(here("Code/DMC_v2.dll"))
  #     if(Linux)dyn.load(here("Code/DMC_v2.so"))
  #     debugonce(migration_crossover)
  #     migration_crossover(param_inds = 1:num_pars, 
  #                         use_theta = theta[, s,, i - 1], 
  #                         use_like = weight[s,, i - 1], 
  #                         data = Data[[s]], 
  #                         hyper = phi[,, i], 
  #                         param_names = param_names
  #                         )
  #   }
  # } 
  
out_array <- foreach(s = 1:S, 
                     .combine = bind_3, 
                     .packages = c("Rcpp"), 
                     .noexport = c("logdens_balba")) %dopar% {
    sourceCpp("BALBA.cpp")
    # debugonce(crossover)
    out <- sapply(1:num_chains,
                  crossover,
                  param_inds = 1:num_pars,
                  use_theta = theta[, s,, i - 1],
                  use_like = weight[s,, i - 1],
                  data = Data[[s]],
                  hyper = phi[,, i],
                  param_names = param_names, 
                  params_fixed = params_fixed, 
                  param_names_fixed = param_names_fixed,
                  mc = num_samples,
                  curr_iter = i,
                  CURRENT_TO_BEST = CURRENT_TO_BEST
    )
  }
  
  out_array <- aperm(out_array, c(1, 3, 2))
  weight[,, i] <- out_array[1,,]
  theta[,,, i] <- out_array[1+1:num_pars,,]
  # weight[,, i] <- out_array[1, ]
  # theta[,,, i] <- out_array[1+1:num_pars, ]
}

#stuff to do when stopped (optional)
save.image("RH_end_v002_A.RData")
rho1a = phi[13,1,1:i]
rho1c = exp(rho1a)/(1+exp(rho1a))
mrho1a = mean(phi[13,,1:i])
mrho1c = exp(mrho1a)/(1+exp(mrho1a))
mrho1c
#

I <- i
toc <- tictoc::toc()

theta[, , 1, 25]


### burnin done by eye to get to ergodic sequence
burnin <- 1
keep <- burnin:I
matplot(t(weight[1,, keep]))
matplot(t(hyper_weight[, keep]))

save.image("RH_end_v002_B.RData")
# ---------------------------------------- Setup for Posterior Analysis
library(psych)
library(tidyverse)
prior_draws <- 10000

theta_prior_list <- list()
for(p in names(prior_list)){
  draws_vec <- numeric(prior_draws)
  for(d in 1:prior_draws){
    mu <- rnorm(1, prior_list[[p]]$mu[1], prior_list[[p]]$mu[2])
    sigma <- rgamma(1, prior_list[[p]]$sigma[1], prior_list[[p]]$sigma[2])
    draws_vec[d] <- rnorm(1, mu, sigma)
  }
  theta_prior_list[[p]] <- draws_vec
}


# ---------------------------------------- Phi plotting
hypers_means <- NULL
for (p in param_names) {
  hypers_means <- c(hypers_means, prior_list[[p]][[1]][1])
  hypers_means <- c(hypers_means, prior_list[[p]][[1]][2])
}
names(hypers_means) <- phi_names


phi_scatter <- array(dim = c(length(phi_names), num_chains * length(keep)))
for (p in 1:length(phi_names)) {
  phi_scatteri <- NULL
  for (c in 1:num_chains) {
    phi_scatteri <- c(phi_scatteri, phi[p, c, keep])
  }
  phi_scatter[p, ] <- phi_scatteri
}
dimnames(phi_scatter)[[1]] <- phi_names
pairs.panels(t(phi_scatter[c(1,3,5,7,9,11,13),])) # just the means

for (p in seq_along(phi_names)) {
  matplot(t(phi[p, , keep]), 
          main = paste(phi_names[p], "full"), type = "l")
  
  
  if (p %% 2 == 1) {
    abline(h = hypers_means[p], col = "green", cex = 2)
    print(prior_list[[ceiling(p / 2)]][[1]][2])
    print(sd(t(phi[, p, keep])))
    # lines(x=x, y=dnorm(x, unlist(prior)[prior_indx[1]], unlist(prior)[prior_indx[2]]), col="red")
  } else {
    abline(h = 1, col = "green", cex = 2)
    # lines(x=x, y=dgamma(x, unlist(prior)[prior_indx[1]], unlist(prior)[prior_indx[2]]), col="red")
  }
  
  xmin <- min(phi[p,, keep], na.rm = T)
  xmax <- max(phi[p,, keep], na.rm = T)
  x <- seq(xmin, xmax, by = (xmax - xmin) / 100)
  prior_indx <- grep(p, names(unlist(prior_list)))
  hist((phi[p,, keep]), main = paste(phi_names[p], "full"), col = "black", freq = F)
  # abline(v=hypers_means[l], col=c("hot pink"))
  if (p %% 2 == 1) {
    abline(v = hypers_means[p], col = "green", cex = 2)
    lines(x=x, y=dnorm(x, prior_list[[ceiling(p / 2)]]$mu[1], prior_list[[ceiling(p / 2)]]$mu[2]), col="red")
  } else {
    abline(v = 1, col = "green", cex = 2)
    lines(x=x, y=dgamma(x, prior_list[[ceiling(p / 2)]]$sigma[1], prior_list[[ceiling(p / 2)]]$sigma[2]), col="red")
  }
}

#--------------------------- Prediction
# Synth_data <- list()
# s <- 1
# Synth_data <- predict_from_chain(
#   num_posts = num_posts, 
#   Synth_data = Synth_data, 
#   num_per = 1, 
#   MEDIAN = T,
#   MEAN = T, 
#   MAP = T, 
#   PHI = T,
#   s = s
# )
# 
# Data_mean_subj <- as_tibble(Data[[1]])
# for (d in 2:length(Data)) {
#   tb <- as_tibble(Data[[d]])
#   Data_mean_subj <- rbind(Data_mean_subj, tb)
# }

# ------------------------------------Plot CDFs
# plot_all_CDFs(data = Data_mean_subj, 
#               Synth_data = Synth_data)


# ---------------------------------------- theta
s <- 1
theta_scatter <- array(dim = c(num_pars, num_chains * length(keep)))
for (p in 1:num_pars) {
  theta_scatteri <- NULL
  for (c in 1:num_chains) {
    theta_scatteri <- c(theta_scatteri, theta[p, s, c, keep])
  }
  theta_scatter[p, ] <- theta_scatteri
}
dimnames(theta_scatter)[[1]] <- param_names
pairs.panels(t(theta_scatter))

#### decide which subject to examine
s <- 1
####
for (p in 1:num_pars) {
  matplot(t(theta[p, s,,]),
          main = paste(param_names[p], "full"), type = "l",
          col = c("black", "hot pink")
  )
  abline(h = params_true[s, p], col = "green", cex = 2)
  print(params_true[s, p])
  
  # lines(y=rep(median(theta_prior_list[[p]]), length(1:max(keep))), x=1:max(keep), col="green")
  hist(theta[p, s,,], 
       main = paste(param_names[p], "full"), 
       col = "black", 
       freq = F)
  abline(v = params_true[s, p], col = "green", cex = 2)
  lines(density(theta_prior_list[[p]]), col="red")
  
  matplot(t(theta[p, s,, keep]),
          main = paste(param_names[p], "samples"), 
          type = "l",
          col = c("black", "hot pink")
  )
  abline(h = params_true[s, p], col = "green", cex = 2)
  # lines(y=rep(median(theta_prior_list[[p]]), length(1:max(keep))), x=1:max(keep), col="green")
  hist(theta[p, s,, keep], 
       main = paste(param_names[p], "samples"), 
       col = "black", 
       freq = F)
  abline(v = params_true[s, p], col = "green", cex = 2)
  lines(density(theta_prior_list[[p]]), col="red")
  
  invisible(readline(prompt = "Press [enter] to continue"))
}
#  invisible(readline(prompt="Press [enter] to continue"))

save.image("RH_end_v002_C.RData")
#--------------------------- Prediction
# Synth_data <- list()
# s <- 15
# Synth_data <- predict_from_chain(num_posts = num_posts, 
#                                  Synth_data = Synth_data,
#                                  num_per = 1, 
#                                  MEDIAN = T, 
#                                  MEAN = T, 
#                                  MAP = F, 
#                                  PHI = F,
#                                  s = s)
# 
# # ------------------------------------Plot CDFs
# plot_all_CDFs(data = Data[[s]], Synth_data = Synth_data)
# 
# print(paste("Subject", s))
# #-------------------------------------------- delta Plots
# make_delta_plots(
#   num_posts = num_posts, 
#   Synth_data = Synth_data, 
#   data = Data[[s]], 
#   time_scale = time_scale,
#   s = s, 
#   MEAN = F, 
#   MEDIAN = F, 
#   MAP = F, 
#   PREDICT = T
)

# invisible(readline(prompt = "Press [enter] to continue"))
# }
