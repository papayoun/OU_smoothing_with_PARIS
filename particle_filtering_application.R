rm(list = ls())
source("generating_data.R") # Get true parameters
source("particle_filtering_utils.R") # Get useful functions
ou_data <- read.table("ou_data.txt",  header = TRUE, sep = ";")


n_particles <- 200
n_backward <- 2
my_H <- .95

truth <- online_kalman_smoothing(ou_data, "ou") %>% 
  pull(Prediction_somme)
biased_kalman <- online_kalman_smoothing(ou_data, density_method = "euler",
                                        obs_H_ = my_H) %>% 
  select(obs_index, t, Prediction_somme) %>% 
  mutate(method = "Biased Kalman",
         difference = Prediction_somme - truth,
         Replicate = "0")
library(parallel)
set.seed(123)
paris <- mclapply(1:30, 
                  function(i){
                    online_paris_smoothing(data_ = ou_data, n_particles = n_particles, 
                                           n_backward = n_backward, 
                                           density_method = "ou")$predictions %>% 
                      mutate(method = "PaRIS",
                             difference = Prediction_somme - truth)
                  },
                  mc.cores = parallel::detectCores() - 1) %>% 
  bind_rows(.id = "Replicate")
biased_paris <- mclapply(1:30, 
                         function(i){
                           online_paris_smoothing(data_ = ou_data, n_particles = n_particles, 
                                                  n_backward = n_backward,
                                                  obs_H_ = my_H,
                                                  density_method = "euler")$predictions %>% 
                             mutate(method = "Biased PaRIS",
                                    difference = Prediction_somme - truth)
                         },
                         mc.cores = parallel::detectCores() - 1)  %>% 
  bind_rows(.id = "Replicate")
# Remove all functions to only keep 
results <- bind_rows(biased_kalman, paris, biased_paris) %>% 
  mutate(Smoother = factor(method,
                           levels = c("PaRIS", "Biased PaRIS", "Biased Kalman"))) %>% 
  select(-method)
rm(biased_kalman, biased_paris, paris)
rm(list = lsf.str())
save.image("experiments_results.RData")


  