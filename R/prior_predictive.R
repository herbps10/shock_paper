library(invgamma)
library(tidyverse)
library(truncnorm)

source("R/plot_theme.R")

simulate_regularized_horseshoe <- function(N, scale_global = 1, local_df = 1, slab_df = 1, s = 1) {
  tau <- abs(rt(1, df = 1) * scale_global)
  lambda <- rt(N, df = local_df)
  csq <- rinvgamma(1, 0.5 * slab_df, 0.5 * slab_df * s^2)
  
  lambda_tilde = sqrt(csq * lambda^2 / (csq + tau^2 * lambda^2))
  beta <- rtruncnorm(N, mean = 0, sd = tau * lambda_tilde, a = 0, b = Inf)
  #beta <- rnorm(N, 0, sd = tau * lambda_tilde)
  beta
}

simulate_csq <- function(N, slab_df, s) {
  rinvgamma(N, 0.5 * slab_df, 0.5 * slab_df * s^2)
}

ggdist::pstudent_t(100, df = 5, sigma = 10, lower.tail = FALSE) * 2 * 100
ggdist::pstudent_t(50, df = 5, sigma = 10, lower.tail = FALSE) * 2 * 100
ggdist::pstudent_t(30, df = 5, sigma = 10, lower.tail = FALSE) * 2 * 100
ggdist::pstudent_t(20, df = 5, sigma = 9, lower.tail = FALSE) * 2 * 100
ggdist::pstudent_t(10, df = 5, sigma = 9, lower.tail = FALSE) * 2 * 100

obj <- function(x) {
  abs(ggdist::pstudent_t(20, df = x[1], sigma = x[2], lower.tail = FALSE) * 2 - 0.1) +
    abs(ggdist::pstudent_t(100, df = x[1], sigma = x[2], lower.tail = FALSE) * 2)
}

pars <- optim(c(5, 10), obj)$par

pars <- round(pars)

ggdist::pstudent_t(20, df = pars[1], sigma = pars[2], lower.tail = FALSE) * 2
ggdist::pstudent_t(50, df = pars[1], sigma = pars[2], lower.tail = FALSE) * 2
ggdist::pstudent_t(threshold, df = pars[1], sigma = pars[2], lower.tail = FALSE) * 2

curve(ggdist::dstudent_t(x, df = pars[1], sigma = pars[2]) * 2, 0, 100)

set.seed(2352)
N <- 2632
results <- expand_grid(
  scale_global = 10^seq(from = -1, to = -3, by = -0.1),
  s            = c(10),
  slab_df      = c(6),
  index        = 1:5e3
) %>%
  mutate(simulations = pmap(list(scale_global, s, slab_df), function(scale_global, s, slab_df) {
    simulate_regularized_horseshoe(1e3, scale_global = scale_global, s = s, slab_df = slab_df)
  }))

threshold <- 2 * fit$samples$summary("epsilon_scale")$median

above_threshold <- results %>%
  mutate(above_threshold = map_dbl(simulations, function(x) mean(abs(x) > threshold))) %>%
  group_by(scale_global) %>%
  summarize(above_threshold = mean(above_threshold))

above_threshold %>%
  ggplot(aes(x = scale_global, y = above_threshold)) +
  geom_point() +
  scale_x_log10() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = expression(tau[0]), y = expression(P(delta[ct] > delta^"*")))
ggsave("plots/prior_predictive_above_threshold.pdf", width = 7, height = 3)

results %>%
  unnest(simulations) %>%
  filter(abs(simulations) > threshold) %>%
  group_by(scale_global) %>%
  summarize(above_20 = mean(abs(simulations) > 20),
            above_100 = mean(abs(simulations) > 100))

results %>%
  unnest(simulations) %>%
  filter(simulations > threshold) %>%
  group_by(scale_global) %>%
  ggplot(aes(x = simulations)) +
  geom_histogram(color = "white") +
  facet_wrap(~scale_global, scales = "free_y")
