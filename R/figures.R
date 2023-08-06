epsilon_scale <- tidybayes::spread_draws(fit$samples$draws("epsilon_scale"), epsilon_scale) %>%
  mutate(model = "No shocks") %>%
  bind_rows(
    tidybayes::spread_draws(fits$fit[[1]]$samples$draws("epsilon_scale"), epsilon_scale) %>%
      mutate(model = "Shocks")
  )
  
epsilon_scale %>%
  ggplot(aes(x = epsilon_scale)) +
  geom_density(aes(color = model)) +
  geom_function(fun = function(x) invgamma::dinvgamma(x, 0.1, 0.1), aes(color = "Prior")) +
  xlim(c(0, 2))

#
# Long-term projections
#
ci_width_comparison <- fit$posteriors$temporal %>%
  filter(variable == "eta", year == 2095) %>%
  mutate(ci_width_no_shocks = `97.5%` - `2.5%`) %>%
  select(name, year, ci_width_no_shocks) %>%
  left_join(
    fits$fit[[1]]$posteriors$temporal %>%
      filter(variable == "eta", year == 2095) %>%
      mutate(ci_width_shocks = `97.5%` - `2.5%`) %>%
      select(name, year, ci_width_shocks)    
  )

projection_comparison <- fit$posteriors$temporal %>%
  filter(year == 2095) %>%
  mutate(median_no_shocks = `50%`) %>%
  select(name, year, median_no_shocks) %>%
  left_join(
    fits$fit[[1]]$posteriors$temporal %>%
      filter(year == 2095, variable == "eta") %>%
      mutate(median_shocks = `50%`) %>%
      select(name, year, median_shocks)
  )

posterior_median_labels <- projection_comparison %>%
  filter(median_shocks < 70 | median_shocks > 90)
  
p1 <- projection_comparison %>%
  ggplot(aes(x = median_no_shocks, y = median_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  #geom_text(data = posterior_median_labels, aes(label = name), hjust = 0.5, nudge_y = 0.6, nudge_x = 0.8, size = 3) +
  pub_theme +
  labs(x = "no shocks",
       y = "shocks",
       subtitle = "posterior median")

ci_width_labels <- ci_width_comparison %>%
  filter(ci_width_shocks > 19.4)

p2 <- ci_width_comparison %>%
  ggplot(aes(x = ci_width_no_shocks, y = ci_width_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  #geom_text(data = ci_width_labels, aes(label = name), hjust = 1, nudge_x = -0.1, nudge_y = 0, size = 3) +
  pub_theme +
  labs(x = "no shocks",
       y = "shocks",
       subtitle = "95% credible interval width")

(p1 + p2) + plot_annotation("Male period life expectancy by country, 2095-2100", tag_levels = "A")
ggsave("plots/life_projection_comparison.pdf", height = 4, width = 10)

#
# Fit comparisons
#
eta <- fit$posteriors$temporal %>%
  filter(variable == "eta") %>%
  mutate(model = "No shocks") %>%
  bind_rows(
    fits$fit[[2]]$posteriors$temporal %>%
      filter(variable == "eta") %>%
      mutate(model = "Shocks") 
  )

countries <- projection_comparison %>% 
  mutate(abs_diff = abs(median_no_shocks - median_shocks)) %>%
  arrange(abs_diff) %>%
  pull(name)

# Pick top and bottom 3
countries <- c(countries[1:4], countries[(length(countries) - 3):length(countries)])

eta %>%
  filter(name %in% countries) %>%
  mutate(name = factor(name, levels = countries)) %>%
  ggplot(aes(x = year + 2.5, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = model), color = "transparent", alpha = 0.2) +
  geom_line(aes(color = model)) + 
  geom_point(data = fit$data %>% filter(name %in% countries) %>% mutate(name = factor(name, levels = countries)), aes(y = e0), alpha = 0.5) +
  facet_wrap(~name, nrow = 2) +
  labs(x = "Year", y = expression(e[0]))

ggsave("plots/life_fit_examples.pdf", width = 10, height = 5)

#
# Compare transition functions
#

plot_mean_transition(fit)
plot_mean_transition(fits$fit[[1]])
plot_mean_transition(fits$fit[[2]])
plot_mean_transition(fits$fit[[3]])

#
# Prior/posterior plots
#

stan_data <- fits$fit[[1]]$stan_data

# Prior on c
caux <- fits %>%
  mutate(caux = map(fit, function(fit) {
    spread_draws(fit$samples$draws("caux"), caux)
  })) %>%
  select(scale_global, caux) %>%
  unnest(c(caux))

caux %>%
  ggplot(aes(x = caux, color = factor(scale_global))) +
  geom_density() +
  geom_function(aes(color = "Prior"), fun = function(x) invgamma::dinvgamma(x, 0.5 * stan_data$slab_df, 0.5 * stan_data$slab_df))

# Prior on tau0
global_shrinkage <- fits %>%
  mutate(global_shrinkage = map(fit, function(fit) {
    spread_draws(fit$samples$draws("global_shrinkage"), global_shrinkage)
  })) %>%
  select(scale_global, global_shrinkage) %>%
  unnest(c(global_shrinkage))

global_shrinkage %>%
  ggplot(aes(x = global_shrinkage, color = factor(scale_global))) +
  geom_density() +
  labs(x = expression(tau), y = "Posterior distribution\n(kernel density estimate)", color = expression(tau[0])) +
  ggtitle("Global scale posterior distribution")
  #geom_function(aes(color = "Prior"), fun = function(x) 2 * ggdist::dstudent_t(x, 1, mu = 0, sigma = stan_data$scale_global))

ggsave("plots/global_scale_posterior.pdf", width = 7, height = 3)

#
# Largest shocks
#

largest_shocks <- function(fit) {
  threshold_shocks <- 2 * fit$samples$summary("epsilon_scale")$median
  print(threshold_shocks)
  fit$posteriors$temporal %>%
    filter(variable == "shock") %>%
    arrange(-abs(`50%`)) %>%
    filter(abs(`97.5%`) > threshold_shocks)
}

largest_shocks(fits$fit[[1]]) %>% select(name, year, `2.5%`, `50%`, `97.5%`)
largest_shocks(fits$fit[[2]]) %>% select(name, year, `2.5%`, `50%`, `97.5%`)
largest_shocks(fits$fit[[3]]) %>% select(name, year, `2.5%`, `50%`, `97.5%`)

countries <- largest_shocks(fits$fit[[2]]) %>% select(name, year, `2.5%`, `50%`, `97.5%`) %>% pull(name) %>% unique()

plot_shock_corrected(fits$fit[[2]], countries[1:6])
ggsave("plots/largest_shocks.pdf", width = 8, height = 4)

plot_shock_corrected(fits$fit[[2]], countries)
ggsave("plots/all_largest_shocks.pdf", width = 10, height = 8)

#
# Plot all
#

countries <- sort(unique(fit$data$name))

pdf("plots/life_all_countries.pdf", width = 10, height = 4)
for(country in countries) {
  print(country)
  p1 <- plot_indicator(fit, country) +
    pub_theme +
    labs(x = "Year", y = expression(e[0]), subtitle = "no shocks")
  
  p2 <- plot_indicator(fit2, country) +
    pub_theme +
    labs(x = "Year", y = expression(e[0]), subtitle = expression(shocks~(tau[0]==0.01)))
  
  p <- (p1 + p2) + plot_layout(guides = "collect")
  print(p)
}
dev.off()