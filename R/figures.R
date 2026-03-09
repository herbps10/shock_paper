epsilon_scale <- tidybayes::spread_draws(fit_noshock$samples$draws("epsilon_scale"), epsilon_scale) |>
  mutate(model = "No shocks") |>
  bind_rows(
    tidybayes::spread_draws(fit_shock$samples$draws("epsilon_scale"), epsilon_scale) |>
      mutate(model = "Shocks")
  )
  
epsilon_scale |>
  ggplot(aes(x = epsilon_scale)) +
  geom_density(aes(color = model)) +
  geom_function(fun = function(x) invgamma::dinvgamma(x, 0.1, 0.1), aes(color = "Prior")) +
  xlim(c(0, 2))

#
# Long-term projections
#

fit_shock <- fits$fit[[1]]
fit_noshock <- fits$fit[[2]]

ci_width_comparison <- fit_noshock$posteriors$temporal |>
  filter(variable == "eta", year == 2100) |>
  mutate(ci_width_no_shocks = `90%` - `10%`) |>
  select(name, year, ci_width_no_shocks) |>
  left_join(
    fit_shock$posteriors$temporal |>
      filter(variable == "eta", year == 2100) |>
      mutate(ci_width_shocks = `90%` - `10%`) |>
      select(name, year, ci_width_shocks)    
  )

projection_comparison <- fit_noshock$posteriors$temporal |>
  filter(year == 2100, variable == "eta") |>
  mutate(median_no_shocks = `50%`) |>
  select(name, year, median_no_shocks) |>
  left_join(
    fit_shock$posteriors$temporal |>
      filter(year == 2100, variable == "eta") |>
      mutate(median_shocks = `50%`) |>
      select(name, year, median_shocks)
  )

posterior_median_labels <- projection_comparison |>
  filter(median_shocks < 70 | median_shocks > 90)
  
p1 <- projection_comparison |>
  ggplot(aes(x = median_no_shocks, y = median_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  #geom_text(data = posterior_median_labels, aes(label = name), hjust = 0.5, nudge_y = 0.6, nudge_x = 0.8, size = 3) +
  pub_theme +
  labs(x = "no shocks",
       y = "shocks",
       subtitle = "posterior median")

ci_width_labels <- ci_width_comparison |>
  filter(ci_width_shocks > 19.4)

p2 <- ci_width_comparison |>
  ggplot(aes(x = ci_width_no_shocks, y = ci_width_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  #geom_text(data = ci_width_labels, aes(label = name), hjust = 1, nudge_x = -0.1, nudge_y = 0, size = 3) +
  pub_theme +
  labs(x = "no shocks",
       y = "shocks",
       subtitle = "80% credible interval width")

(p1 + p2) + plot_annotation("Male period life expectancy by country, 2100", tag_levels = "A")
ggsave("plots/life_projection_comparison.pdf", height = 4, width = 10)

#
# Fit comparisons
#
eta <- fit_noshock$posteriors$temporal |>
  filter(variable == "eta") |>
  mutate(model = "No shocks") |>
  bind_rows(
    fit_shock$posteriors$temporal |>
      filter(variable == "eta") |>
      mutate(model = "Shocks") 
  )

countries <- projection_comparison |> 
  mutate(abs_diff = abs(median_no_shocks - median_shocks), diff = median_no_shocks - median_shocks) |>
  arrange(abs_diff)


# Pick top and bottom 3
#arranged_countries <- c(countries$name[1:4], arrange(countries, diff)$name[1:4], arrange(countries, -diff)$name[1:4])
arranged_countries <- c(arrange(countries, abs_diff)$name[1:4], arrange(countries, -abs_diff)$name[1:4])

eta |>
  filter(name %in% arranged_countries) |>
  mutate(name = factor(name, levels = arranged_countries)) |>
  ggplot(aes(x = year + 2.5, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = model), color = "transparent", alpha = 0.2) +
  geom_line(aes(color = model)) + 
  geom_point(data = fit_shock$data |> filter(name %in% arranged_countries) |> mutate(name = factor(name, levels = arranged_countries)), aes(y = e0), alpha = 0.5) +
  facet_wrap(~name, nrow = 2, scales = "free_y") +
  labs(x = "Year", y = expression(e[0])) +
  pub_theme +
  theme(legend.position = "bottom")

ggsave("plots/life_fit_examples.pdf", width = 10, height = 4)

lower_shocks_countries <- projection_comparison |>
  mutate(abs_diff = abs(median_no_shocks - median_shocks), diff = median_no_shocks - median_shocks) |>
  filter(diff > 1) |>
  arrange(diff) |>
  pull(name)

eta |>
  filter(name %in% lower_shocks_countries) |>
  mutate(name = factor(name, levels = lower_shocks_countries)) |>
  ggplot(aes(x = year + 2.5, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = model), color = "transparent", alpha = 0.2) +
  geom_line(aes(color = model)) + 
  geom_point(data = fit_shock$data |> filter(name %in% lower_shocks_countries) |> mutate(name = factor(name, levels = lower_shocks_countries)), aes(y = e0), alpha = 0.5) +
  facet_wrap(~name, nrow = 4, scales = "free_y") +
  labs(x = "Year", y = expression(e[0])) +
  pub_theme +
  theme(legend.position = "bottom")

plot_comparison(fit_shock, "Uganda")

plot_transition(fit_shock,   "Lesotho")
plot_transition(fit_noshock, "Lesotho")

areas <- c("Kenya", "Uganda", "Zimbabwe", "Haiti")
data <- fit_shock$data |>
  filter(name %in% areas) |>
  arrange(name, period) |>
  group_by(name) |>
  mutate(diff = c(diff(e0), NA))

fit_shock$posteriors$transition_functions |>
  filter(name %in% areas, .width == 0.5) |>
  ggplot(aes(x = 15 + x * (110 - 15), y = transition_function_pred)) +
  geom_line(aes(color = "Shock")) +
  geom_line(aes(color = "No Shock"), data = filter(fit_noshock$posterior$transition_functions, name %in% areas, .width == 0.5)) +
  scale_fill_brewer() +
  geom_point(data = data, aes(x = e0, y = diff)) +
  facet_wrap(~name)

fit_shock$posteriors$temporal |>
  filter(variable == "shock") |>
  filter(name == "Kenya") |>
  ggplot(aes(x = year, y= `50%`)) +
  geom_line()

fits$fit[[1]]$posteriors$temporal |>
  filter(name %in% areas, variable == "eta") |>
  ggplot(aes(x = year, y = `50%`, color = "Shock")) +
  geom_line() +
  geom_line(aes(color = "No shock"), data = fits$fit[[2]]$posteriors$temporal |> filter(name %in% areas, variable == "eta")) +
  facet_wrap(~name)

fit_shock$posteriors$temporal |>
  filter(name %in% areas, variable == "eta") |>
  ggplot(aes(x = year, y = `50%`, color = "Shock")) +
  geom_line() +
  geom_line(aes(color = "No shock"), data = fit_noshock$posteriors$temporal |> filter(name %in% areas, variable == "eta")) +
  facet_wrap(~name)


#
# Compare transition functions
#

p1 <- plot_mean_transition(fit_noshock) + ylim(c(0, 11)) + labs(x = expression(e0), y = expression(f[b])) +
  theme(legend.position = "none") +
  ggtitle(label = "No shocks")
p2 <- plot_mean_transition(fit_shock) + ylim(c(0, 11)) + labs(x = expression(e0), y = expression(f[b])) +
  theme(legend.position = "none") +
  ggtitle(label = "Shocks")

p1 / p2

ggsave("plots/transition_function_comparisons.pdf", width = 8, height = 5)

fit_noshock$posteriors$transition_function_mean |> mutate(model = "No shock") |>
  bind_rows(fit_shock$posteriors$transition_function_mean |> mutate(model = "Shock")) |>
  ggplot(aes(x = 15 + (x * (110 - 15)), y = transition_function_mean)) +
  geom_point(data = datM_diffs, aes(x = e0, y = diff), size = 0.1) +
  geom_line(aes(color = model)) +
  geom_smooth(se = FALSE)
  

#
# Prior/posterior plots
#

stan_data <- fit_shock$stan_data

# Prior on c
caux <- fits |>
  mutate(caux = map(fit, function(fit) {
    spread_draws(fit$samples$draws("caux"), caux)
  })) |>
  select(scale_global, caux) |>
  unnest(c(caux))

caux |>
  ggplot(aes(x = caux, color = factor(scale_global))) +
  geom_density() +
  geom_function(aes(color = "Prior"), fun = function(x) invgamma::dinvgamma(x, 0.5 * stan_data$slab_df, 0.5 * stan_data$slab_df))

# Prior on tau0
global_shrinkage <- fits |>
  filter(model == "shock2") |>
  mutate(global_shrinkage = map(fit, function(fit) {
    spread_draws(fit$samples$draws("global_shrinkage"), global_shrinkage)
  })) |>
  select(scale_global, global_shrinkage) |>
  unnest(c(global_shrinkage))

global_shrinkage |>
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
  fit$posteriors$temporal |>
    filter(variable == "shock") |>
    arrange(-abs(`50%`)) |>
    filter(abs(`97.5%`) > threshold_shocks)
}

largest_shocks(fit_shock) |> select(name, year, `2.5%`, `50%`, `97.5%`)

countries <- largest_shocks(fit_shock) |> select(name, year, `2.5%`, `50%`, `97.5%`) |> pull(name) |> unique()

plot_shock_corrected(fit_shock, countries[1:6]) +
  theme(legend.position = "bottom")
ggsave("plots/largest_shocks.pdf", width = 8.5, height = 4)

plot_shock_corrected(fit_shock, countries)
ggsave("plots/all_largest_shocks.pdf", width = 10, height = 8)

#
# Plot all
#

countries <- sort(unique(fit_shock$data$name))

pdf("plots/life_all_countries_yearly.pdf", width = 10, height = 4)
for(country in countries) {
  print(country)
  p1 <- plot_indicator(fit_noshock, country) +
    pub_theme +
    labs(x = "Year", y = expression(e[0]), subtitle = "no shocks")
  
  p2 <- plot_indicator(fit_shock, country) +
    pub_theme +
    labs(x = "Year", y = expression(e[0]), subtitle = expression(shocks~(tau[0]==0.01)))
  
  p <- (p1 + p2) + plot_layout(guides = "collect")
  print(p)
}
dev.off()
