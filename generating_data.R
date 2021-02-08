mu <- 5; rho <- 1; sigma <- 1
sigma1 <- 1
m1 <- 0 # Mean of initial state
delta <- 1
n_obs <- 51
require(tidyverse)
set.seed(123)
source("particle_filtering_utils.R")
ou_data <- tibble(t = seq(0, length.out = n_obs, by = delta)) %>% 
  mutate(x = get_ou_samples(times_ = t, x0_ = rnorm(1, m1, sigma1)))
sigma_obs <- ou_data %>% 
  pull(x) %>% 
  diff() %>% 
  abs(.) %>% 
  mean() %>% 
  {. * 2} 
sigma_obs <- 1

ou_data <- ou_data %>% 
  mutate(y = rnorm(n_obs, mean = x, sd = sigma_obs))
ou_data %>% 
  write.table(file = "ou_data.txt", row.names = FALSE, 
              col.names = TRUE, sep = ";")
rm(ou_data)
