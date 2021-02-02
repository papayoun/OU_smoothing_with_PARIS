rm(list = ls())
# Utils -------------------------------------------------------------------

get_q <- function(x, a, b, x0){
  dnorm(x, a * x0, b)
}
get_g <- function(x, c, d, y){
  dnorm(y, c * x, d)
}
get_a <- function(rho, delta, method = "ou"){
  if(method == "ou"){
    return(exp(-rho * delta))
  }
  else if(method == "euler"){
    return(1 - rho * delta)
  }
  else{
    stop("method must be 'ou' or 'euler'")
  }
}
get_b <- function(rho, sigma, delta, method = "ou"){
  if(method == "ou"){
    return(sigma * sqrt(0.5 * (1 - exp(-2 * rho * delta)) / rho))
  }
  else if(method == "euler"){
    return(sigma * sqrt(delta))
  }
  else{
    stop("method must be 'ou' or 'euler'")
  }
}

get_l <- function(x, a, b, c, d, x0, y){
  get_q(x, a, b, x0) * get_g(x, c, d, y)
}

get_gamma2 <- function(b, c, d){
  gamma2_inv <- c^2 / d^2 + 1 / b^2
  1 / gamma2_inv
}

get_integrand <- function(rho, sigma, delta, c, d, x0, y, method = "ou"){
  a <- get_a(rho, delta, method); 
  b <- get_b(rho, sigma, delta, method)
  gamma2 <- get_gamma2(b, c, d)
  initial_integrand <- sqrt(gamma2) / (b * d * sqrt(2 * pi))
  eterm <- exp(-0.5 * ((y / d)^2 + (a / b * x0)^2 - gamma2 * (a / b^2 * x0 + c / d^2 * y)^2))
  if(eterm > 1)
    warning("eterm superieur à 1")
  initial_integrand * eterm
}


# Example -----------------------------------------------------------------

set.seed(123)
foo <- function(my_x0, my_y, method){
  z <- get_integrand(rho = 0.5, sigma = 1, delta = 2 , c = 1, d = 1, my_x0, my_y, method)
  data.frame(x0 = my_x0, y = my_y, z = z, method = method)
}

plan_xp <- expand.grid(my_x0 = seq(-4, 4, length.out = 101),
                       my_y = seq(-4, 4, length.out = 101),
                       method = c("euler", "ou"))
library(tidyverse)
resultats <- purrr::pmap_dfr(plan_xp, foo)
resultats %>% 
  group_by(x0, y) %>% 
  summarise(diff = abs(diff(z))) %>% 
  ggplot() +
  aes(x = x0, y = y, fill = diff) +
  geom_raster() +
  scale_fill_viridis_c()

my_alpha <- runif(1)
my_y <- runif(1, -10, 10)
foo2 <- function(x, y){
  dnorm(y, x, 1) - dnorm(y, my_alpha * x, 1)
}
foo2_grad <- function(x, y){
  (y - x) * dnorm(y, x, 1) - my_alpha * (y - my_alpha * x) *dnorm(y, my_alpha * x, 1)
}
optimize(function(x) abs(foo2(x, my_y)), c(-10, 10), maximum = TRUE)

tibble(x = seq(-10, 10, length.out = 501), y = my_y, g = dnorm(y, x, 1), geps = dnorm(y, my_alpha * x, 1),
       diff = abs(foo2(x, y)), grad = foo2_grad(x, y), bound = abs(x *y)) %>% 
  pivot_longer(-c("x", "y"), names_to = "Fonction", values_to = "Valeur") %>% 
  ggplot(aes(x = x, y = Valeur, color = Fonction)) +
    geom_line() + 
  geom_hline(yintercept = 0)
