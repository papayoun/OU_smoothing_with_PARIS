
# Libairies ---------------------------------------------------------------

rm(list = ls())

library(tidyverse) # Data manipulation and plots


# Theme set ---------------------------------------------------------------

my_theme <- function(base_size = 8){
  theme_bw()  %+replace%
    theme(
      panel.border = element_rect(colour = "black", 
                                  fill = rgb(0, 0, 0, 0)),
      # plot.background = element_rect(fill = "white"),# bg around panel
      legend.background = element_blank(), 
      text = element_text(family = "LM Roman 10", size = base_size),
      axis.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      legend.title = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      strip.background = element_rect(fill = "lightgoldenrod1",
                                      color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1)),
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5))
}
theme_set(my_theme(base_size = 12))


# output_set --------------------------------------------------------------

out_width <- 14 # cm
out_height <- 2 * out_width / 3 # format 3/2


# Varying epsilon, fixed n ------------------------------------------------

load("experiments_results_vary_epsilon_fixed_n.RData")
plot1 <- results %>% 
  ggplot() +
  aes(x = epsilon, y = abs(difference), color = Smoother) +
  geom_point(data = dplyr::filter(results, Smoother == "PaRIS")) +
  geom_path(data = dplyr::filter(results, Smoother == "Kalman")) +
  labs(x = expression(epsilon),
       y = expression("|"~phi[0:n]^epsilon~h[n]~"-"~phi[0:n]~h[n]~"|"))
ggsave(filename = "vary_epsilon_fixed_n.png", 
         plot = plot1,
         width = out_width, height = out_height, units = "cm",
         device = "png")

# Fixed epsilon, varying n ------------------------------------------------

load("experiments_results_fixed_epsilon_vary_n.RData")
plot2 <- results %>% 
  ggplot() +
  aes(x = obs_index, y = difference, color = Smoother) +
  geom_path(aes(group = interaction(Smoother, Replicate)), 
                data = dplyr::filter(results, str_detect(Smoother, "PaRIS"))) +
  geom_path(data = dplyr::filter(results, Smoother == "Biased Kalman")) +
  labs(x = expression(n),
       y = expression(phi[0:n]^epsilon~h[n]~"-"~phi[0:n]~h[n])) +
  scale_color_discrete(labels = c("Kalman,"~epsilon==0.1,
                                  "PaRIS,"~epsilon==0.1,
                                  "PaRIS,"~epsilon==0))
ggsave(filename = "fixed_epsilon_vary_n.png", 
       plot = plot2,
       width = out_width, height = out_height, units = "cm",
       device = "png")
