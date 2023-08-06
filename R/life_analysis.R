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

# Histogram of changes in e0
e0_differences <- datM %>%
  group_by(name) %>%
  mutate(diff = c(NA, diff(e0)))

e0_differences %>%
  ggplot(aes(x = diff)) +
  geom_histogram(color = "white", binwidth = 1) +
  geom_boxplot(aes(y = -20), width = 30, alpha = 0.5) +
  labs(x = expression(paste("Difference in ", e[0], ": ", eta[ct] - eta[ct-1])), y = "Count")

countries <- c("Republic of Korea", "Bosnia and Herzegovina", "Cambodia", "Lebanon", "Timor-Leste", "Syrian Arab Republic")

datM %>%   
  filter(name %in% countries) %>%
  ggplot(aes(x = year + 2.5, y = e0)) +
  geom_point() +
  geom_line(alpha = 0.5) +
  facet_wrap(~name) +
  pub_theme +
  labs(x = "Year", y = expression(e[0]))

ggsave("plots/life_examples.pdf", height = 5, width = 10)

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
  num_knots = 7, 
  hierarchical_splines = c("intercept", "name"),
  
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1e3,
  
  extra_stan_data = list(
    scale_global = 1e-2,
    slab_scale = 10,
    slab_df = 6
  )
)

print(fit$samples$cmdstan_diagnose())

# What is the 2*SD(eps) threshold?
threshold <- 2 * fit$samples$summary("epsilon_scale")$median

# How many of the observed differences fall below this threshold?
mean(e0_differences$diff < -threshold, na.rm = TRUE)

fits <- expand_grid(
  scale_global = c(1e-3, 1e-2, 1e-1),
  num_knots = c(7)
) %>%
  mutate(fit = map(scale_global, function(scale_global) {
    lifeplus(
      datM,
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,
      
      model = "shock2",
      
      spline_degree = 2,
      num_knots = 7, 
      hierarchical_splines = c("intercept", "name"),
      
      adapt_delta = 0.99,
      max_treedepth = 12,
      parallel_chains = 4,
      iter_warmup = 250,
      iter_sampling = 1e3,
      
      extra_stan_data = list(
        scale_global = scale_global,
        slab_scale = 10,
        slab_df = 6
      )
    )
  }))

#BayesTransitionModels:::plot_indicator(fits_with_shocks$fit[[1]], areas = "Somalia")
#BayesTransitionModels:::plot_temporal("eta_crisisfree", fits_with_shocks$fit[[1]], areas = "Somalia")
#BayesTransitionModels:::plot_temporal("shock", fits_with_shocks$fit[[1]], areas = "Cambodia")
#
#name <- "Algeria"
#BayesTransitionModels:::plot_temporal("eta", fit, name)
#BayesTransitionModels:::plot_temporal("eta", fits_with_shocks$fit[[1]], name)
#BayesTransitionModels:::plot_temporal("eta_crisisfree", fits_with_shocks$fit[[1]], name)
#plot_shock(fits_with_shocks$fit[[1]], name)

#
# Validations
#

validation_cutoff <- function(model, cutoff_year, scale_global, num_knots) {
  fit <- lifeplus(
    datM,
    y = "e0", 
    year = "year",
    area = "name",
    source = "source",
    start_year = 1950,
    end_year = 2100,
    held_out = datM$year >= cutoff_year,
    
    model = model,
    
    spline_degree = 2,
    num_knots = num_knots, 
    hierarchical_splines = c("intercept", "name"),
    
    adapt_delta = 0.99,
    max_treedepth = 12,
    parallel_chains = 4,
    iter_warmup = 250,
    iter_sampling = 500,
    
    extra_stan_data = list(
      scale_global = scale_global,
      slab_scale = 10,
      slab_df = 6
    )
  )
  print(fit$samples$cmdstan_diagnose())
  return(fit)
}

validations <- bind_rows(validations, expand_grid(
  model = c("spline", "shock2"),
  cutoff_year = c(2005, 2010, 2015),
  scale_global = c(1e-2),
  num_knots = c(9)
) %>%
  mutate(fit = pmap(list(model, cutoff_year, scale_global, num_knots), validation_cutoff)))

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
           error = e0 - `50%`)
}

validation_results <-  validations %>%
  filter(model != "shock") %>%
  mutate(validation = map(fit, validation_measures))

plot_indicator(validations$fit[[2]], "Iraq")

validation_results %>%
  filter(cutoff_year == 2015) %>%
  select(-fit) %>%
  unnest(validation) %>%
  select(model, name, error) %>%
  pivot_wider(names_from = "model", values_from = "error") %>%
  ggplot(aes(x = spline, y = shock2)) +
  geom_point() +
  geom_abline(lty = 2)

validation_results_summary <- validation_results %>%
  mutate(validation = map(validation, function(validation) {
    validation %>%
      summarize(below = mean(below),
              above = mean(above),
              coverage = mean(covered),
              ci_width = median(`97.5%` - `2.5%`),
              median_error = median(error),
              mean_squared_error = mean(error^2),
              median_abs_error = median(abs(error)))
  })) %>%
  select(-fit) %>%
  unnest(cols = c(validation)) %>%
  arrange(cutoff_year, model, scale_global)

validation_table <- validation_results_summary %>%
  select(cutoff_year, model, scale_global, num_knots, below, coverage, above, ci_width, median_error, median_abs_error) %>%
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

validation_table

validation_table %>%
  select(-scale_global) %>%
  knitr::kable(format = "latex")

validation_table %>%
  knitr::kable(format = "latex")
