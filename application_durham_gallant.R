rm(list = ls())
source("generating_data.R") # Get true parameters
source("particle_filtering_utils.R") # Get useful functions
ou_data <- read.table("ou_data.txt",  header = TRUE, sep = ";")
tail(ou_data)
x0_ <- rnorm(2, mu)
xF_ <- rnorm(2, mu)
n_interval <- 10
durham_params_ <- list(epsilon = delta / n_interval, M = 50)
get_transition_density(xF_, x0_, delta_= delta, method = "ou", durham_params_ = durham_params_)
get_transition_density(xF_, x0_, delta_= delta, method = "euler", durham_params_ = durham_params_)
get_transition_density(xF_, x0_, delta_= delta, method = "durham", durham_params_ = durham_params_)
n_particles <- 10
n_backward <- 2

ex_ou <- online_paris_smoothing(data_ = ou_data, n_particles = n_particles, 
                       n_backward = n_backward, 
                       density_method = "ou")
ex_durham <- online_paris_smoothing(data_ = ou_data[1:5, ], n_particles = n_particles, 
                                    n_backward = n_backward, 
                                    density_method = "durham", durham_params_ = durham_params_)
