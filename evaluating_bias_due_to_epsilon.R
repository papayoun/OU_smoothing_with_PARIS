rm(list = ls())
source("generating_data.R") # Get true parameters
source("particle_filtering_utils.R") # Get useful functions
ou_data <- read.table("ou_data.txt",  header = TRUE, sep = ";")


n_particles <- 200
n_backward <- 2
alphas <- seq(0.5, 1.5, length.out = 41)

truth <- online_kalman_smoothing(ou_data, "ou") %>% 
  pull(Prediction_somme)
biased_kalman <- map_dfr(alphas, function(my_H){
  online_kalman_smoothing(ou_data, density_method = "ou",
                          obs_H_ = my_H) %>% 
    mutate(method = "Biased Kalman",
           difference = Prediction_somme - truth,
           epsilon = abs(1 - H),
           signe_epsilon = ifelse(H > 1, "+", "-"),
           Replicate = "0") %>% 
    slice(nrow(.))
})

ggplot(biased_kalman) +
  aes(x = epsilon, y = abs(difference),
      color = signe_epsilon) +
  geom_point()
set.seed(123)
paris <- mclapply(1:30,
                  function(i){
                    map_dfr(alphas, function(my_alpha){
                      online_paris_smoothing(data_ = ou_data, 
                                             n_particles = n_particles, 
                                             n_backward = n_backward, 
                                             obs_H_ = my_alpha,
                                             density_method = "ou")$predictions %>% 
                        mutate(method = "PaRIS",
                               difference = Prediction_somme - truth,
                               H = my_alpha,
                               epsilon = abs(1 - H),
                               signe_epsilon = ifelse(H >=1, "+", "-")) %>% 
                        slice(nrow(.))
                    })
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

results <- bind_rows(biased_kalman, 
                     mutate(paris, signe_epsilon = ifelse(H > 1, "+", "-")) ) %>% 
  mutate(Smoother = factor(method,
                           labels = c("Kalman", "PaRIS"),
                           levels = c("Biased Kalman", "PaRIS"))) %>% 
  select(-method)
results_epsilon_neg <- results %>% 
  dplyr::filter(signe_epsilon == "-")
results_epsilon_neg %>% 
  ggplot() +
  aes(x = epsilon, y = abs(difference), color = Smoother) +
  geom_point(data = dplyr::filter(results_epsilon_neg, Smoother == "PaRIS")) +
  geom_path(data = dplyr::filter(results_epsilon_neg, Smoother == "Kalman")) +
  labs(x = expression(epsilon), y = expression("|"~L^epsilon~"-"~L~"|"))
