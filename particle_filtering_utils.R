# Parameters

# Functions

get_transition_moments <- function(x0_, delta_, method_, obs_ = NULL){
  if(method_ == "euler"){
    mean <- x0_ - rho * (x0_ - mu) * delta_
    sd <- sigma * sqrt(delta_)
  }
  else if(method_ == "ou"){
    mean <- exp(-rho * delta_) * (x0_ - mu) + mu
    sd <- sigma * sqrt(.5   / rho * (1 - exp(-2 * rho * delta)))
  }
  else if(stringr::str_detect(method_, "proposal")){
    # Dynamics moments
    if(method_ == "proposal_ou"){
      mean_dyn <- exp(-rho * delta_) * (x0_ - mu) + mu # Mean (OU)
      var_dyn <- .5 * sigma^2 / rho * (1 - exp(-2 * rho * delta)) # Variance (OU)
    }
    else if(method_ == "proposal_euler"){
      mean_dyn <- x0_ - rho * (x0_ - mu) * delta_ # Mean (Euler)
      var_dyn <- sigma^2 * delta_ # Variance (Euler)
    }
    else{
      print(paste0("Your method is ", method_))
      stop("Unknown method")
    }
    # Observation moments
    mean_obs <- obs_ # Observation model mean
    var_obs <- sigma_obs * sigma_obs # Observation model variance
    # Proposal moments, obtained by conjugation
    var_proposal <- var_dyn * var_obs / (var_dyn + var_obs) 
    mean <- var_proposal * (mean_dyn / var_dyn + mean_obs / var_obs)
    sd <- sqrt(var_proposal)
  }
  else
    stop("Unknown method, must be `ou`, `euler` or `proposal`")
  list(mean = mean, sd = sd)
}

get_transition_density <- function(xF_, x0_, delta_, method_, obs_ = NULL,
                                   log_ = FALSE){
  moments <- get_transition_moments(x0_, delta_, method_, obs_)
  dnorm(xF_, moments[["mean"]], moments[["sd"]], log = log_)
}

get_ou_samples <- function(times_, x0_ = m0){
  n_obs <- length(times_)
  x <- rep(NA, n_obs)
  x[1] <- x0_
  for(t in 2:n_obs){
    delta <- times_[t] - times_[t - 1]
    moments <- get_transition_moments(x0_ = x[t - 1], delta_ = delta, method_ = "ou")
    x[t] <- rnorm(n = 1, mean =  moments[["mean"]], sd = moments[["sd"]])
  }
  x
}

get_proposal_sample <- function(x0_, delta_, obs_, method_){
  moments <- get_transition_moments(x0_, delta_, obs_, method_ = method_)
  rnorm(length(x0_), moments[["mean"]], moments[["sd"]])
}

get_initial_moments <- function(method_, obs_ = NULL){
  if(method_ == "initial"){
    mean = m1
    sd = sigma1
  }
  else if(method_ == "proposal"){
    mean_dyn <- m1
    var_dyn <- sigma1 * sigma1
    mean_obs <- obs_
    var_obs <- sigma_obs * sigma_obs
    var_proposal <- var_dyn * var_obs / (var_dyn + var_obs) 
    mean <- var_proposal * (mean_dyn / var_dyn + mean_obs / var_obs)
    sd <- sqrt(var_proposal)
  }
  list(mean = mean, 
       sd = sd)
}

get_initial_density <- function(x0_, method_, obs_ = NULL, log_ = FALSE){
  moments <-  get_initial_moments(method_, obs_)
  dnorm(x0_, moments[["mean"]], moments[["sd"]], log = log_)
}

get_initial_proposal_samples <- function(n_, obs_){
  moments <- get_initial_moments(method = "proposal", obs_ = obs_)
  rnorm(n_, moments[["mean"]], moments[["sd"]])
}

# perform_particle_filter <- function(data_, n_particles, n_backward,
#                                     density_method){
#   observations <- pull(data_, y)
#   obs_times <- pull(data_, t)
#   n_obs <- nrow(data_)
#   # Creation of needed quantities
#   particles_set <- weights_set  <- matrix(nrow = n_particles, ncol = n_obs)
#   taus_set <- array(0, dim = c(n_particles, 2, n_obs))
#   predictions_table <- map_dfr(1:n_obs,
#                                function(t){
#                                  tibble(t_state = t,
#                                         t_obs = t:n_obs,
#                                         Prediction = rep(NA, n_obs - t + 1))
#                                })
#   # Initialization
# 
#   particles_set[, 1] <- get_initial_proposal_samples(n_particles, observations[1])
#   weights_set[, 1] <- (get_initial_density(x0_ = particles_set[, 1],
#                                            method_ = "initial") * dnorm(observations[1], particles_set[, 1], sigma_obs) /
#                          get_initial_density(x0_ = particles_set[, 1],
#                                              method_ = "proposal",
#                                              obs_ = observations[1])) %>%
#     {. / sum(.)}
#   predictions_table$Prediction[1] <- sum(particles_set[, 1] * weights_set[, 1])
#   for(t in 2:n_obs){
#     # Particle filtering
# 
#     delta_t <- obs_times[t] - obs_times[t - 1]
#     selected_indexes <- sample(1:n_particles,
#                                size = n_particles, replace = TRUE,
#                                prob = weights_set[, t - 1])
#     selected_particles <- particles_set[selected_indexes, t - 1]
#     new_particles <- get_proposal_sample(x0_ = selected_particles, delta_ = delta_t,
#                                          obs_ = observations[t])
#     weights_set[, t] <- (get_transition_density(xF_ = new_particles,
#                                                 x0_ = selected_particles,
#                                                 delta_ = delta_t,
#                                                 method = density_method) *
#                            dnorm(observations[t], new_particles, sigma_obs) /
#                            get_transition_density(xF_ = new_particles,
#                                                   x0_ = selected_particles,
#                                                   delta_ = delta_t,
#                                                   method = "proposal",
#                                                   obs = observations[t])) %>%
#       {. / sum(.)}
#     particles_set[, t] <- new_particles
#     predictions_table[predictions_table$t_state == t &
#                         predictions_table$t_obs == t, "Prediction"] <- sum(particles_set[,t] * weights_set[, t])
# 
#     # backward sampling
#     for(i in 1:n_particles){
#       sigma_plus <- max(get_transition_density(xF_ = particles_set[i, t],
#                                                x0_ = particles_set[, t - 1],
#                                                delta_ = delta_t,
#                                                method = density_method))
#       for(j in 1:n_backward){
#         accepted <- FALSE
#         while(!accepted){
#           proposed_ancestor <- sample(1:n_particles,
#                                       size = 1,
#                                       prob = weights_set[,t - 1])
#           accept_ratio <- get_transition_density(xF_ = particles_set[i, t],
#                                                  x0_ = particles_set[proposed_ancestor, t - 1],
#                                                  delta_ = delta_t,
#                                                  method = density_method) / sigma_plus
#           accepted <- (runif(1) < accept_ratio)
#         }
#         taus_set[i, 2, 0:(t - 1)] <- taus_set[i, 2, 0:(t - 1)] +
#           (taus_set[proposed_ancestor, 1, 0:(t - 1)] +
#              ((1:(t-1))==(t-1)) * particles_set[proposed_ancestor, t - 1]) / n_backward
#       }
#     }
#     taus_set[,1,1:(t-1)] <- taus_set[,2,1:(t-1)]
#     taus_set[,2,1:(t-1)] <- 0
#     for(ell in 1:(t-1)){
#       predictions_table[predictions_table$t_state == ell &
#                           predictions_table$t_obs == t,
#                         "Prediction"] <- sum(taus_set[, 1, ell] * weights_set[, t])
#     }
#   }
#   list(particles_set = particles_set,
#        weights_set = weights_set,
#        predictions = predictions_table)
# }

online_paris_smoothing <- function(data_, n_particles, n_backward,
                                   density_method,
                                   obs_H_ = 1,
                                   smoothing_lag = 1){
  proposal <- paste0("proposal_", density_method)
  observations <- pull(data_, y)
  obs_times <- pull(data_, t)
  n_obs <- nrow(data_)
  smoothing_indexes <- seq(n_obs, 2, by = -smoothing_lag)
  # Creation of needed quantities
 
  particles_set <- weights_set  <- matrix(nrow = n_particles, ncol = n_obs)
  taus_old <- taus_new <- rep(0, n_particles)
  predictions_table <- data.frame(obs_index = c(1, rev(smoothing_indexes))) %>% 
    mutate(t = obs_times[obs_index]) %>% 
    mutate(Prediction_somme = NA)
  # Initialization
  
  particles_set[, 1] <- get_initial_proposal_samples(n_particles, observations[1])
  weights_set[, 1] <- (get_initial_density(x0_ = particles_set[, 1],
                                           method_ = "initial") * 
                         dnorm(observations[1], obs_H_ * particles_set[, 1], sigma_obs) /
                         get_initial_density(x0_ = particles_set[, 1],
                                             method_ = "proposal",
                                             obs_ = observations[1])) %>%
    {. / sum(.)}
  predictions_table$Prediction_somme[1] <- sum(particles_set[, 1] * weights_set[, 1])
  for(t in 2:n_obs){
    # Particle filtering
    delta_t <- obs_times[t] - obs_times[t - 1]
    selected_indexes <- sample(1:n_particles,
                               size = n_particles, replace = TRUE,
                               prob = weights_set[, t - 1])
    selected_particles <- particles_set[selected_indexes, t - 1]
    new_particles <- get_proposal_sample(x0_ = selected_particles, delta_ = delta_t,
                                         obs_ = observations[t], 
                                         method_ = proposal)
    weights_set[, t] <- (get_transition_density(xF_ = new_particles,
                                                x0_ = selected_particles,
                                                delta_ = delta_t,
                                                method = density_method) *
                           dnorm(observations[t], obs_H_ * new_particles, sigma_obs) /
                           get_transition_density(xF_ = new_particles,
                                                  x0_ = selected_particles,
                                                  delta_ = delta_t,
                                                  method = proposal,
                                                  obs = observations[t])) %>%
      {. / sum(.)}
    particles_set[, t] <- new_particles
    # backward sampling
    for(i in 1:n_particles){
      sigma_plus <- max(get_transition_density(xF_ = particles_set[i, t],
                                               x0_ = particles_set[, t - 1],
                                               delta_ = delta_t,
                                               method = density_method))
      for(j in 1:n_backward){
        accepted <- FALSE
        while(!accepted){
          proposed_ancestor <- sample(1:n_particles,
                                      size = 1,
                                      prob = weights_set[,t - 1])
          accept_ratio <- get_transition_density(xF_ = particles_set[i, t],
                                                 x0_ = particles_set[proposed_ancestor, t - 1],
                                                 delta_ = delta_t,
                                                 method = density_method) / sigma_plus
          accepted <- (runif(1) < accept_ratio)
        }
        taus_new[i] <- taus_new[i] + 
          (taus_old[proposed_ancestor] + particles_set[proposed_ancestor, t - 1]) / n_backward
      }
    } # End loops over particles for backward sampling
    taus_old <- taus_new
    taus_new <- rep(0, n_particles)
    # Now compute the estimate (must add the last particle)
    if(t %in% smoothing_indexes){
      to_change <- which(predictions_table$obs_index == t)
      predictions_table$Prediction_somme[to_change] <- sum((taus_old + particles_set[, t]) * weights_set[, t])
    }
  }
  list(particles_set = particles_set,
       weights_set = weights_set,
       predictions = predictions_table)
}



perform_kalman_smoothing <- function(data_, density_method){
  observations <- pull(data_, y)
  obs_times <- pull(data_, t)
  n_obs <- nrow(data_)
  KF_gain <- KF_mean <- 
    KF_variance <- KF_prediction_var <- rep(NA, n_obs)
  predictions_table <- map_dfr(1:n_obs,
                               function(t){
                                 tibble(t_state = t,
                                        t_obs = t:n_obs,
                                        Prediction = rep(NA, n_obs - t + 1))
                               })
  # Initialization
  # Model X_{t+delta} = A*X(t) + b + N(0, Q)
  # Y(t) = H * X(t) + N(0, sigma_obs^2)
  if(density_method == "ou"){
    dyn_A <- exp(-rho * delta)
    dyn_b <- (1 - exp(-rho * delta)) * mu
    dyn_Q <- .5 * sigma^2 / rho * (1 - exp(-2 * rho *delta))
    obs_H <- 1
  }
  else if(density_method == "euler"){
    dyn_A <- (1 - rho * delta)
    dyn_b <- rho * delta * mu
    dyn_Q <- sigma^2 * delta
    obs_H <- .9
  }
  
  
  KF_gain[1] <- sigma1^2 * obs_H / (obs_H^2 * sigma1^2 + sigma_obs^2)
  KF_mean[1] <- m1  + 
    KF_gain[1] * (observations[1] - obs_H * m1)
  KF_variance[1] <- (1 - KF_gain[1] * obs_H) * sigma1^2
  KF_prediction_var[1] <- dyn_A^2 * KF_variance[1] + dyn_Q
  for(t in 2:n_obs){
    KF_gain[t] <- KF_prediction_var[t - 1] * obs_H / 
      (obs_H^2 * KF_prediction_var[t - 1] + sigma_obs^2)
    KF_mean[t] <- dyn_A * KF_mean[t - 1] + dyn_b +
      KF_gain[t] * (observations[t] - obs_H * (dyn_A * KF_mean[t - 1] + dyn_b))
    KF_variance[t] <- (1 - KF_gain[t] * obs_H) * KF_prediction_var[t - 1]
    KF_prediction_var[t] <- dyn_A^2 * KF_variance[t] + dyn_Q
  }
  predictions_table <- predictions_table %>% 
    mutate(Prediction = ifelse(t_state == t_obs, KF_mean[t_state], NA))
  for(t_final in n_obs:2){
    KS_mean <- KS_variance <- KS_gain <- rep(NA, t_final)
    KS_mean[t_final] <- KF_mean[t_final]
    KS_variance[t_final] <- KF_variance[t_final]
    for(t in (t_final - 1):1){
      KS_gain[t] <- KF_variance[t] * dyn_A / KF_prediction_var[t]
      KS_mean[t] <- KF_mean[t] + 
        KS_gain[t] * (KS_mean[t + 1] - (dyn_A * KF_mean[t] + dyn_b))
      KS_variance[t] <- KF_variance[t] + 
        KS_gain[t]^2 * (KS_variance[t + 1] - KF_prediction_var[t])
      predictions_table$Prediction[predictions_table$t_state == t & 
                                     predictions_table$t_obs == t_final] <- KS_mean[t]
    }
  }
  # Smoothing
  predictions_table
}



online_kalman_smoothing <- function(data_, density_method, 
                                    smoothing_lag = 1,
                                    obs_H_ = 1){
    observations <- pull(data_, y)
    obs_times <- pull(data_, t)
    n_obs <- nrow(data_)
    smoothing_indexes <- seq(n_obs, 2, by = -smoothing_lag)
    KF_gain <- KF_mean <- 
      KF_variance <- KF_prediction_var <- rep(NA, n_obs)
    predictions_table <- tibble(obs_index = c(1, rev(smoothing_indexes))) %>% 
      mutate(t = obs_times[obs_index]) %>% 
      mutate(Prediction_somme = NA)
    # Initialization
    # Model X_{t+delta} = A*X(t) + b + N(0, Q)
    # Y(t) = H * X(t) + N(0, sigma_obs^2)
    if(density_method == "ou"){
      dyn_A <- exp(-rho * delta)
      dyn_b <- (1 - exp(-rho * delta)) * mu
      dyn_Q <- .5 * sigma^2 / rho * (1 - exp(-2 * rho *delta))
      obs_H <- obs_H_
    }
    else if(density_method == "euler"){
      dyn_A <- (1 - rho * delta)
      dyn_b <- rho * delta * mu 
      dyn_Q <- sigma^2 * delta
      obs_H <- obs_H_
    }
    
    
    KF_gain[1] <- sigma1^2 * obs_H / (obs_H^2 * sigma1^2 + sigma_obs^2)
    KF_mean[1] <- m1  + 
      KF_gain[1] * (observations[1] - obs_H * m1)
    KF_variance[1] <- (1 - KF_gain[1] * obs_H) * sigma1^2
    KF_prediction_var[1] <- dyn_A^2 * KF_variance[1] + dyn_Q
    for(t in 2:n_obs){
      KF_gain[t] <- KF_prediction_var[t - 1] * obs_H / 
        (obs_H^2 * KF_prediction_var[t - 1] + sigma_obs^2)
      KF_mean[t] <- dyn_A * KF_mean[t - 1] + dyn_b +
        KF_gain[t] * (observations[t] - obs_H * (dyn_A * KF_mean[t - 1] + dyn_b))
      KF_variance[t] <- (1 - KF_gain[t] * obs_H) * KF_prediction_var[t - 1]
      KF_prediction_var[t] <- dyn_A^2 * KF_variance[t] + dyn_Q
    }
    predictions_table$Prediction_somme[1] <- KF_mean[1]
    for(final_index in smoothing_indexes){
      KS_mean <- KS_variance <- KS_gain <- rep(NA, final_index)
      KS_mean[final_index] <- KF_mean[final_index]
      KS_variance[final_index] <- KF_variance[final_index]
      for(t in (final_index - 1):1){
        KS_gain[t] <- KF_variance[t] * dyn_A / KF_prediction_var[t]
        KS_mean[t] <- KF_mean[t] + 
          KS_gain[t] * (KS_mean[t + 1] - (dyn_A * KF_mean[t] + dyn_b))
        KS_variance[t] <- KF_variance[t] + 
          KS_gain[t]^2 * (KS_variance[t + 1] - KF_prediction_var[t])
      }
      predictions_table$Prediction_somme[predictions_table$obs_index == final_index] <- sum(KS_mean)
      if(final_index == n_obs){
        predictions_table$Prediction_etat_lissage = KS_mean
      }
    }
    # Smoothing
    predictions_table %>% 
      mutate(Prediction_etat_filtrage = KF_mean) %>% 
      mutate(method = density_method,
             H = obs_H_)
}
