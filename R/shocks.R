library(BayesTransitionModels)
library(tidyverse)

data <- national_data(
  countries = c("Somalia", "Afghanistan", "India", "Bangladesh", "Zimbabwe",
                "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala", 
                "Rwanda", "Timor-Leste", "Viet Nam", "Ghana", "Nigeria",
                "Ethiopia", "Uganda", "Serbia", "Ukraine", "Libya", "Lesotho", "Swaziland"),
  start_year = 1970
)

my_fit <- function(scale_global = 0.01, slab_scale = 1, slab_df = 1, model = here::here("./stan/fpem_spline_shock.stan")) { 
  fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",
    
    start_year = 1970, 
    end_year = 2030, 
    
    model = model,
    
    smoothing = TRUE,
    
    # Spline setup
    spline_degree = 2,
    num_knots = 7,
    
    # Hierarchical setup
    hierarchical_level     = c("intercept", "name_country"), 
    hierarchical_splines   = c("intercept", "name_country"),
    hierarchical_asymptote = c("intercept", "name_country"),
    
    # Prior settings
    tau_prior = "normal(0, 2)",
    rho_prior = "uniform(0, 1)",
    
    extra_stan_data = list(
      scale_global = scale_global,
      slab_scale = slab_scale,
      slab_df = slab_df
    ),
    
    # Stan sampler settings
    adapt_delta = 0.95,
    max_treedepth = 10,
    iter_warmup = 250,
    iter_sampling = 500,
    seed = 5,
    parallel_chains = 4,
    refresh = 50
  )
}

basic_fit <- my_fit(model = "spline")
fit <- my_fit(0.01, model = here::here("./stan/fpem_spline_shock.stan"))
fit2 <- my_fit(1e-4, model = here::here("./stan/fpem_spline_shock.stan"))

plot_indicator(basic_fit)
plot_indicator(fit)
plot_indicator(fit2)

shock <- tidybayes::spread_draws(fit$samples, shock[c, t])
shock2 <- tidybayes::spread_draws(fit2$samples, shock[c, t])

-shock$shock %>% log %>% hist()
-shock2$shock %>% log %>% hist()
-shock3$shock %>% log %>% hist()

mean(-shock$shock > 0.05)
mean(-shock2$shock > 0.05)

hist(shock$shock[shock$shock < 0.001 & shock$shock > -0.1])

plot_indicator(fit2, areas = "Timor-Leste")
plot_smoother(fit2, areas = "Timor-Leste")

plot_smoothing_hyperparameters(basic_fit)
plot_smoothing_hyperparameters(fit2)
plot_smoothing_hyperparameters(fit2_obs_nonse)

plot_transition(basic_fit, "name_country")
plot_transition(fit2, "name_country")

eta_uncertainty <- function(fit, target_year) fit$posteriors$temporal %>%
  filter(year == target_year, variable == "eta") %>%
  mutate(ci_width = `97.5%` - `2.5%`)
  
eta_uncertainty(basic_fit, 2030) %>%
  mutate(model = "Basic") %>%
  bind_rows(
    eta_uncertainty(fit2, 2030) %>%
      mutate(model = "Horsehoe Shocks")
  ) %>%
  ggplot(aes(x = ci_width)) + geom_histogram(color = "white") +
  facet_wrap(~model)

eta_uncertainty(basic_fit, 2030) %>%
  select(name_country, ci_width) %>%
  rename(ci_width_basic = ci_width) %>%
  left_join(
    eta_uncertainty(fit2, 2030) %>%
      select(name_country, ci_width) %>%
      rename(ci_width_shock = ci_width)
  ) %>%
  ggplot(aes(x = ci_width_basic, y = ci_width_shock)) +
  geom_point() +
  geom_text(aes(label = name_country)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed()

plot_data_hyperparameters(fit1)
plot_data_hyperparameters(fit2)
plot_data_hyperparameters(fit3)

plot_shock <- function(fit, areas = c(), n = 0) {
  shock <- tidybayes::spread_draws(fit$samples, shock[c, t])
  
  dat <- shock %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(fit$time_index) %>%
    left_join(fit$country_index)
  
  if(n > 0) {
    sample_paths <- shock %>%
      group_by(.iteration, .draw, .chain) %>%
      nest() %>%
      ungroup() %>%
      sample_n(n) %>%
      mutate(index = 1:n()) %>%
      unnest() %>%
      left_join(fit$time_index) %>%
      left_join(fit$country_index)
  }
  
  if(length(areas) > 0) {
    dat <- dat %>% filter(name_country %in% areas)
    if(n > 0) {
      sample_paths <- sample_paths %>% filter(name_country %in% areas)
    }
  }
  
  p <- dat %>%
    ggplot(aes(x = year, y = shock)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper), size = 0.5) +
    scale_fill_brewer() +
    facet_wrap(~name_country)
  
  if(n > 0) {
    p <- p + geom_line(data = sample_paths, aes(group = index))
  }
  
  p
}

plot_shock(fit1)
plot_shock(fit2)
plot_shock(fit3)

plot_indicator(fit2, "Timor-Leste")

plot_shock(fit, "Rwanda", n = 5)
plot_shock(fit2, "Rwanda", n = 5)

plot_indicator(fit2, "Indonesia")
plot_shock(fit2, "Indonesia", n = 5)

plot_smoothing_hyperparameters(basic_fit)
plot_smoothing_hyperparameters(fit2)

country <- "Rwanda"
plot_indicator(basic_fit, country)
plot_indicator(fit, country)
plot_indicator(fit2, country)
plot_indicator(fit3, country)

basic_fit$posteriors$nonse %>%
  mutate(setup = "Default") %>%
  bind_rows(fit2$posteriors$nonse %>%
              mutate(setup = "Horseshoe shock")) %>%
  group_by(data_series_type, setup) %>%
  median_qi(nonse) %>%
  ggplot(aes(x = nonse, y = data_series_type, color = setup)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0, position = position_dodge(width = 0.25))

basic_fit$posteriors$ar %>%
  mutate(setup = "Default") %>%
  bind_rows(fit2$posteriors$ar %>%
              mutate(setup = "Horseshoe shock")) %>%
  pivot_longer(c(est_rho, est_tau)) %>%
  ggplot(aes(x = value, color = setup)) +
  geom_density() +
  facet_wrap(~name, scales = "free")


