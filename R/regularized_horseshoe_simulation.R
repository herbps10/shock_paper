library(invgamma)
library(tidyverse)
library(truncnorm)

source("R/plot_theme.R")

simulate_horseshoe <- function(N, scale_global = 1) {
  tau <- rt(N, df = 1) * scale_global
  lambda <- rcauchy(N)
  
  beta <- rnorm(N, mean = 0, sd = abs(tau * lambda))
  
  beta
}

simulate_regularized_horseshoe <- function(N, scale_global = 1, local_df = 1, slab_df = 1, s = 1) {
  tau <- abs(rt(1, df = 1) * scale_global)
  lambda <- rt(N, df = local_df)
  csq <- rinvgamma(1, 0.5 * slab_df, 0.5 * slab_df * s^2)
  
  lambda_tilde = sqrt(csq * lambda^2 / (csq + tau^2 * lambda^2))
  beta <- rtruncnorm(N, mean = 0, sd = tau * lambda_tilde, a = 0, b = Inf)
  beta
}

expand_grid(
  slab_df = c(1, 3, 100),
  s = c(0.1, 0.5, 1),
  x = exp(seq(-5, 5, 0.1))
) %>%
  mutate(d = dinvgamma(x, shape = 0.5 * slab_df, rate = 0.5 * slab_df * s^2)) %>%
  ggplot(aes(x = x, y = d)) +
  geom_line() +
  facet_grid(s ~ slab_df, scale = "free_y") +
  scale_x_log10()

results_horseshoe <- expand_grid(
  scale_global = c(1e-2, 1e-3, 1e-4, 1e-5),
  index = 1:200
) %>%
  mutate(simulations = map(scale_global, function(scale_global) {
    simulate_horseshoe(1e3, scale_global = scale_global)
  }))

threshold_results_horseshoe <- results_horseshoe %>%
  mutate(thresholds = list(tibble(threshold = c(1, 0.5, 0.1, 0.01)))) %>%
  unnest(thresholds) %>%
  mutate(above_threshold = map2_dbl(simulations, threshold, function(simulations, threshold) mean(abs(simulations) > threshold)))

threshold_results_horseshoe %>%
  group_by(threshold, scale_global) %>%
  summarize(above_threshold = mean(above_threshold)) %>%
  ggplot(aes(x = scale_global, y = above_threshold, color = factor(threshold))) +
  geom_line() +
  geom_point() + 
  scale_x_log10()

set.seed(2352)
N <- 2632
results <- expand_grid(
  #scale_global = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
  scale_global = c(1e-2, 1e-4, 1e-6),
  s            = c(0.1, 0.5, 1, 10),
  slab_df      = c(1, 3, 100),
  index        = 1:5e3
) %>%
  mutate(simulations = pmap(list(scale_global, s, slab_df), function(scale_global, s, slab_df) {
    simulate_regularized_horseshoe(1e3, scale_global = scale_global, s = s, slab_df = slab_df)
  }))


results %>%
  unnest(simulations) %>%
  filter(simulations > delta_threshold) %>%
  group_by(scale_global, s, slab_df) %>%
  summarize(q = quantile(simulations, 0.99)) %>%
  ggplot(aes(x = scale_global, y = q)) +
  geom_point(aes(color = factor(s))) +
  facet_wrap(~slab_df, scales = "free_x") +
  scale_x_log10(breaks = c(1e-6, 1e-4, 1e-2), sec.axis = sec_axis(~ ., name = expression(nu), breaks = NULL, labels = NULL)) +
  labs(y = expression(99^th~percentile~of~delta[ct],~delta[ct]>delta[ct]^"*"),
       x = expression(Global~scale~tau[0]), y = expression(Average~delta[ct]>delta^"*")) +
  pub_theme +
  theme(panel.border = element_rect(color = "black"))

q <- seq(0.05, 0.95, 0.05)
y <- rtruncnorm(1e6, a = 0)
qn <- quantile(y[y > 0.1], q)

qq <- results %>%
  unnest(simulations) %>%
  filter(simulations > 0.1) %>%
  group_by(scale_global, s, slab_df) %>%
  summarize(q = list(tibble(qemp = quantile(simulations, q), qn))) %>%
  unnest(q)

ggplot(qq, aes(x = qn, y = qemp, color = factor(s))) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  facet_wrap(~slab_df)

results %>%
  unnest(simulations) %>%
  filter(simulations > delta_threshold) %>%
  ggplot(aes(x = simulations, color = factor(scale_global))) +
  geom_density() +
  scale_x_log10(sec.axis = sec_axis(~ ., name = expression(s), breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = expression(nu), breaks = NULL, labels = NULL)) +
  facet_grid(slab_df ~ s, scales = "free_x")
ggsave("plots/prior_predictive_densities.pdf", width = 8, height = 5)

results %>%
  unnest(simulations) %>%
  filter(simulations > delta_threshold) %>%
  ggplot(aes(x = simulations, color = factor(slab_df))) +
  stat_ecdf(n = 50) +
  scale_x_log10(sec.axis = sec_axis(~ ., name = expression(s), breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = expression(tau[0]), breaks = NULL, labels = NULL)) +
  facet_grid(scale_global ~ s) +
  labs(x = expression(delta), y = "ECDF") 

ggsave("plots/prior_predictive_ecdf.pdf", width = 8, height = 5)

threshold_results <- results %>%
  mutate(above_threshold = map2_dbl(simulations, threshold, function(simulations, threshold) mean(abs(simulations) > threshold)))

threshold_results %>%
  group_by(scale_global, s, slab_df, threshold) %>%
  summarize(above_threshold = mean(above_threshold)) %>%
  ggplot(aes(x = scale_global, y = above_threshold)) +
  geom_line() +
  geom_point() +
  scale_x_log10(sec.axis = sec_axis(~ ., name = expression(s), breaks = NULL, labels = NULL)) +
  scale_y_log10(sec.axis = sec_axis(~ ., name = expression(nu), breaks = NULL, labels = NULL), labels = scales::percent) +
  facet_grid(slab_df ~ s) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = expression(Global~scale~tau[0]), y = expression(Average~delta[ct]>delta^"*")) +
  pub_theme +
  theme(panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/mcpr_prior_predictive.pdf", width = 8, height = 5)

threshold_results %>% 
  filter(scale_global == 1e-2, slab_df == 1, s == 1) %>% 
  group_by(threshold) %>%
  summarize(above_threshold = mean(above_threshold) * 100)