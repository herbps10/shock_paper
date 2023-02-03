library(tidyverse)
library(wpp2022)
library(BayesTransitionModels)
library(tidybayes)
library(bayesLife)
library(patchwork)

source("R/plot_theme.R")
source("R/lifeplus.R")
source("R/process_lifeplus.R")

data(e0M)
data(include_2010)

# Number of countries before filtering
e0M %>%
  filter(country_code < 900) %>%
  distinct(name) %>%
  nrow()

included_codes <- include_2010 %>% filter(include_code %in% 1:2) %>% pull(country_code)

datM <- e0M %>%
  pivot_longer(cols = `1950-1955`:`2015-2020`, names_to = "period", values_to = "e0") %>%
  filter(country_code %in% included_codes) %>%
  mutate(year = parse_integer(str_sub(period, 1, 4)),
         source = "WPP2022")

# Number of countries after filtering
datM %>%
  distinct(name) %>%
  nrow()

# Histogram of year-over-year changes in e0
e0_differences <- datM %>%
  group_by(name) %>%
  mutate(diff = c(NA, diff(e0)))

e0_differences %>%
  ggplot(aes(x = diff)) +
  geom_histogram(color = "white", binwidth = 1)

quantile(e0_differences$diff, na.rm = TRUE, c(0, 0.025, 0.05, 0.5, 0.95, 0.975))
mean(e0_differences$diff < -2, na.rm = TRUE)

countries <- c("Bangladesh", "Bosnia and Herzegovina", "Cambodia", "El Salvador", "Lebanon", "Timor-Leste")
  
datM %>%   
  filter(name %in% countries) %>%
  ggplot(aes(x = year + 2.5, y = e0)) +
  geom_point() +
  geom_line(alpha = 0.5) +
  facet_wrap(~name) +
  pub_theme +
  labs(x = "Year", y = expression(e[0]))

ggsave("plots/life_examples.pdf", height = 7, width = 10)

fit <- lifeplus(
  datM,
  y = "e0", 
  year = "year",
  area = "name",
  source = "source",
  start_year = 1950,
  end_year = 2100,
  
  model = "spline",
  
  spline_degree = 2,
  num_knots = 5, 
  hierarchical_splines = c("intercept", "name"),
  
  parallel_chains = 4,
  iter_warmup = 250,
  iter_sampling = 500,
  
  extra_stan_data = list(
    scale_global = 1e-2,
    slab_scale = 1,
    slab_df = 1
  )
)

fit1 <- lifeplus(
  datM,
  y = "e0", 
  year = "year",
  area = "name",
  source = "source",
  start_year = 1950,
  end_year = 2100,
  
  model = "shock",
  
  spline_degree = 2,
  num_knots = 5, 
  hierarchical_splines = c("intercept", "name"),
  
  parallel_chains = 4,
  iter_warmup = 250,
  iter_sampling = 500,
  
  extra_stan_data = list(
    scale_global = 1e-2,
    slab_scale = 1,
    slab_df = 1
  )
)

fit2 <- lifeplus(
  datM,
  y = "e0", 
  year = "year",
  area = "name",
  source = "source",
  start_year = 1950,
  end_year = 2100,
  
  model = "shock2",
  
  spline_degree = 2,
  num_knots = 5, 
  hierarchical_splines = c("intercept", "name"),
  
  adapt_delta = 0.99,
  parallel_chains = 4,
  iter_warmup = 250,
  iter_sampling = 500,
  
  extra_stan_data = list(
    scale_global = 1e-2,
    slab_scale = 1,
    slab_df = 1
  )
)

fit$samples %>% spread_draws(epsilon_scale) %>% mutate(model = "no shocks") %>%
  bind_rows(
    fit2$samples %>% spread_draws(epsilon_scale) %>% mutate(model = "with shocks")
  ) %>%
  ggplot(aes(x = epsilon_scale, color = model)) +
  geom_density() +
  labs(x = "random walk variance")

ggsave("plots/life_random_walk_variance.pdf", height = 5)


plot_shock <- function(fit, areas = c()) {
  shocks <- fit$samples %>% spread_draws(shock[c, t]) %>%
    left_join(fit$time_index) %>%
    left_join(fit$country_index)
  
  if(length(areas) > 0) {
    shocks <- shocks %>% filter(name %in% areas)
  }
  
  shocks %>%
    group_by(year, name) %>%
    median_qi(shock, .width = c(0.5, 0.8, 0.95)) %>%
    ggplot(aes(x = year, y = shock)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
    scale_fill_brewer() +
    facet_wrap(~name)
}

plot_P_tilde <- function(fit) {
  fit$samples %>% spread_draws(P_tilde2[c]) %>%
    left_join(fit$country_index) %>%
    group_by(name) %>%
    mutate(P_tilde2 = 15 + P_tilde2) %>%
    ggplot(aes(x = P_tilde2)) +
    geom_density() +
    facet_wrap(~name)
}

plot_P_tilde(fit)
plot_P_tilde(fit2)

plot_shock(fit1, countries)
plot_shock(fit2, countries)


ci_width_comparison <- fit$posteriors$temporal %>%
  filter(year == 2095) %>%
  mutate(ci_width_no_shocks = `97.5%` - `2.5%`) %>%
  select(name, year, ci_width_no_shocks) %>%
  left_join(
    fit2$posteriors$temporal %>%
      filter(year == 2095) %>%
      mutate(ci_width_shocks = `97.5%` - `2.5%`) %>%
      select(name, year, ci_width_shocks)    
  )

ci_width_comparison %>%
  ggplot(aes(x = ci_width_no_shocks, y = ci_width_shocks)) +
  geom_point() +
  geom_abline(slope = 1) +
  coord_fixed()

ggsave("plots/life_ci_width_comparison_2100.pdf", height = 5)

ci_width_comparison %>%
  summarize_at(vars(ci_width_no_shocks, ci_width_shocks), median)

countries <- c("Somalia", "Cambodia")
p1 <- plot_indicator(fit, areas = countries) + 
  labs(x = "Year", y = expression(e[0]), subtitle = "no shocks") +
  pub_theme

p2 <- plot_indicator(fit2, areas = countries) + 
  labs(x = "Year", y = expression(e[0]), subtitle = "level shocks") +
  pub_theme

(p1 / p2) + plot_layout(guides = "collect")
ggsave("plots/life_fit_examples.pdf", width = 10, height = 8)

projection_comparison <- fit$posteriors$temporal %>%
  filter(year == 2095) %>%
  mutate(median_no_shocks = `50%`) %>%
  select(name, year, median_no_shocks) %>%
  left_join(
    fit2$posteriors$temporal %>%
      filter(year == 2095) %>%
      mutate(median_shocks = `50%`) %>%
      select(name, year, median_shocks)
  )

projection_comparison %>%
  mutate(diff = abs(median_no_shocks - median_shocks)) %>%
  summarize(median(diff))

posterior_median_labels <- projection_comparison %>%
  filter(median_no_shocks < 60)
  
p1 <- projection_comparison %>%
  ggplot(aes(x = median_no_shocks, y = median_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  geom_text(data = posterior_median_labels, aes(label = name), hjust = 0.2, nudge_y = 0.6, nudge_x = 0.8, size = 3) +
  pub_theme +
  labs(x = "no shocks",
       y = "level shocks",
       subtitle = "posterior median")

ci_width_labels <- ci_width_comparison %>%
  filter(ci_width_shocks > 24 | ci_width_no_shocks > 30)

p2 <- ci_width_comparison %>%
  ggplot(aes(x = ci_width_no_shocks, y = ci_width_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  geom_text(data = ci_width_labels, aes(label = name), hjust = 1, nudge_x = -0.5, size = 3) +
  pub_theme +
  labs(x = "no shocks",
       y = "level shocks",
       subtitle = "95% credible interval width")

(p1 + p2) + plot_annotation("Male period life expectancy by country, 2095-2100", tag_levels = "A")
ggsave("plots/life_projection_comparison.pdf", height = 5, width = 10)

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

# Largest shocks
shocks <- fit2$samples %>% spread_draws(shock[c, t]) %>%
  left_join(fit2$country_index) %>%
  left_join(fit2$time_index)

shocks_summarized <- shocks %>%
  group_by(name, year) %>%
  median_qi(shock) %>%
  arrange(shock)

#
# Validations
#

validation_cutoff <- function(model, cutoff_year, scale_global) {
  lifeplus(
    datM,
    y = "e0", 
    year = "year",
    area = "name",
    source = "source",
    start_year = 1950,
    end_year = 2015,
    held_out = datM$year >= cutoff_year,
    
    model = model,
    
    spline_degree = 2,
    num_knots = 5, 
    hierarchical_splines = c("intercept", "name"),
    
    parallel_chains = 4,
    iter_warmup = 250,
    iter_sampling = 500,
    
    extra_stan_data = list(
      scale_global = 1e-2,
      slab_scale = 1,
      slab_df = 1
    )
  )
}

validations <- expand_grid(
  model = c("spline", "shock", "shock2"),
  cutoff_year = c(2005, 2010, 2015),
  scale_global = c(1e-3, 1e-2, 1e-1)
) %>%
  mutate(fit = pmap(list(model, cutoff_year, scale_global), validation_cutoff))

plot_indicator(validations$fit[[1]], countries)

validation_measures <- function(fit) {
  fit$data %>%
    mutate(held_out = fit$held_out) %>%
    filter(held_out == 1) %>%
    group_by(name) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    left_join(fit$posteriors$temporal %>%
      filter(variable == "eta")) %>%
    mutate(below = e0 < `2.5%`,
           above = e0 > `97.5%`,
           covered = `2.5%` <= e0 & `97.5%` >= e0,
           error = e0 - `50%`) %>%
    summarize(below = mean(below),
              above = mean(above),
              coverage = mean(covered),
              ci_width = median(`97.5%` - `2.5%`),
              median_error = median(error),
              mean_squared_error = mean(error^2),
              median_abs_error = median(abs(error)))
}

validation_table <- validations %>%
  mutate(validation = map(fit, validation_measures)) %>%
  select(-fit) %>%
  unnest(cols = c(validation)) %>%
  arrange(cutoff_year, model, scale_global) %>%
  select(cutoff_year, model, scale_global, below, coverage, above, ci_width, median_error, median_abs_error) %>%
  mutate_at(vars(below, coverage, above, ci_width, median_error, median_abs_error), signif, 3) %>%
  mutate_at(vars(below, coverage, above), `*`, 100) %>%
  mutate_at(vars(below, coverage, above), paste0, "%") %>%
  mutate(model = case_when(
    model == "shock2" ~ "level shocks",
    model == "shock" ~ "rate shocks",
    model == "spline" ~ "no shocks"
  )) %>%
  mutate(cutoff_year = cutoff_year - 5,
         cutoff_year = glue::glue("{cutoff_year}-{cutoff_year + 5}"),
         index = 1:n()) %>%
  rename(Cutoff = cutoff_year,
         Model = model,
         `% Below` = below,
         `% Included` = coverage,
         `% Above` = above,
         `CI Width` = ci_width,
         ME = median_error,
         MAE = median_abs_error
         ) %>%
  #mutate(Cutoff = ifelse(index %% 2 == 1, Cutoff, "")) %>%
  select(-index)

validation_table %>%
  filter(scale_global == 0.01) %>%
  select(-scale_global) %>%
  knitr::kable(format = "latex")

validation_table %>%
  knitr::kable(format = "latex")

fit_validation$samples %>% spread_draws(epsilon_scale) %>% ggplot(aes(x = epsilon_scale)) + geom_density()
fit2_validation$samples %>% spread_draws(epsilon_scale) %>% ggplot(aes(x = epsilon_scale)) + geom_density()

shock_validation_2005 <- validations %>% filter(model == "shock2", cutoff_year == 2005) %>%
  pull(fit) %>%
  .[[1]]

fit_validation_2005 <- validations %>% filter(model == "spline", cutoff_year == 2005) %>%
  pull(fit) %>%
  .[[1]]

errors <- shock_validation_2005$data %>%
  mutate(held_out = shock_validation_2005$held_out) %>%
  filter(held_out == 1) %>%
  group_by(name) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  left_join(shock_validation_2005$posteriors$temporal %>%
              filter(variable == "eta")) %>%
  mutate(error = e0 - `50%`)

errors %>%
  ggplot(aes(x = error)) +
  geom_histogram(color = "white")

overpredictions <- errors %>% arrange(error) %>% select(name, error)
underpredictions <- errors %>% arrange(-error) %>% select(name, error)

plot_indicator(shock_validation_2005, underpredictions$name[1:6])
ggsave("plots/life_underpredictions.pdf", height = 6, width = 9)

shock_validation_2005$samples %>%
  spread_draws(shock[c, t]) %>%
  left_join(shock_validation_2005$time_index) %>%
  left_join(shock_validation_2005$country_index) %>%
  filter(name %in% underpredictions$name[1:6]) %>%
  group_by(name, year) %>%
  median_qi(shock, .width = c(0.8, 0.9, 0.95)) %>%
  ggplot(aes(x = year, y = shock)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer() +
  facet_wrap(~name)
  

plot_indicator(fit_validation_2005, underpredictions$name[1:6])

plot_indicator(shock_validation_2005, overpredictions$name[1:6])
ggsave("plots/life_overpredictions.pdf", height = 6, width = 9)
