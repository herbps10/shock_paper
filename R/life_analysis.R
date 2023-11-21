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
data(pop5)
data(include_2010)


# Number of countries before filtering
e0M %>%
  filter(country_code < 900) %>%
  distinct(name) %>%
  nrow()

included_codes <- include_2010 %>% filter(include_code %in% 1:2) %>% pull(country_code)

large_countries <- pop5 %>%
  filter(`2020` >= 100) %>%
  pull(name)

datM <- e0M %>%
  filter(name %in% large_countries) %>%
  pivot_longer(cols = `1950-1955`:`2015-2020`, names_to = "period", values_to = "e0") %>%
  filter(country_code %in% included_codes) %>%
  mutate(year = parse_integer(str_sub(period, 1, 4)),
         source = "WPP2022")

datM_diffs <- datM %>%
  group_by(name) %>%
  mutate(diff = c(NA, diff(e0)))

datM_diffs %>%
  filter(e0 > 55) %>%
  ggplot(aes(x = e0, y = diff)) +
  geom_point() +
  geom_smooth()

datM %>%
  group_by(name) %>%
  mutate(diff = c(NA, diff(e0)))  %>%
  group_by(period) %>%
  summarize(mean = mean(diff))

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
  scale_global = c(1e-2),
  model = c("spline", "shock2"),
  num_knots = c(7)
) %>%
  mutate(fit = pmap(list(scale_global, model), function(scale_global, model) {
    lifeplus(
      datM,
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,
      
      model = model,
      
      spline_degree = 2,
      num_knots = 7, 
      hierarchical_splines = c("intercept", "name"),
      
      normal_data_model = TRUE,
      data_model_df = 5,
      
      adapt_delta = 0.99,
      max_treedepth = 12,
      parallel_chains = 4,
      iter_warmup = 500,
      iter_sampling = 500,
      
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

validation_cutoff <- function(model, cutoff_year, scale_global, num_knots, normal_data_model, data_model_df) {
  #datM <- datM %>%
  #  #filter(name %in% c("Sri Lanka", "Suriname", "Sweden", "Syrian Arrab Republic", "Thailand", "Trinidad and Tobago", "United Kingdom", "Bosnia and Herzegovina", "Brazil", "Cambodia", "Dem. People's Republic of Korea"))
  #  filter(name %in% c("Dem. People's Republic of Korea"))
  fit <- lifeplus(
    datM,
    y = "e0", 
    year = "year",
    area = "name",
    source = "source",
    start_year = 1950,
    end_year = 2100,
    held_out = datM$year >= cutoff_year,
    
    normal_data_model = normal_data_model,
    data_model_df = data_model_df,
    
    model = model,
    
    spline_degree = 2,
    num_knots = num_knots, 
    hierarchical_splines = c("intercept", "name"),
    
    adapt_delta = 0.999,
    max_treedepth = 12,
    parallel_chains = 4,
    iter_warmup = 500,
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

validations <- expand_grid(
  model = c("spline", "shock2"),
  cutoff_year = c(2005, 2010, 2015),
  #cutoff_year = 2005,
  scale_global = c(1e-2),
  num_knots = c(7),
  #normal_data_model = c(FALSE, TRUE),
  data_model_df = c(3),
  normal_data_model = TRUE,
) %>%
  mutate(fit = pmap(list(model, cutoff_year, scale_global, num_knots, normal_data_model, data_model_df), validation_cutoff))


validation_measures <- function(fit, cutoff_year) {
  fit$data %>%
    mutate(held_out = fit$held_out) %>%
    filter(held_out == 1) %>%
    group_by(name) %>%
    filter(year == max(year)) %>%
    #filter(year == cutoff_year) %>%
    ungroup() %>%
    left_join(fit$posteriors$temporal %>%
      filter(variable %in% c("eta", "eta_crisisfree"))) %>%
    group_by(variable) %>%
    mutate(below = e0 < `2.5%`,
           above = e0 > `97.5%`,
           covered = `2.5%` <= e0 & `97.5%` >= e0,
           below0.1 = e0 < `10%`,
           above0.9 = e0 > `90%`,
           covered0.8 = `10%` <= e0 & `90%` >= e0,
           error = e0 - `50%`)
}

validation_results <-  validations %>%
  filter(model != "shock") %>%
  mutate(validation = map2(fit, cutoff_year, validation_measures))

validation_results_summary <- validation_results %>%
  mutate(validation = map(validation, function(validation) {
    validation %>%
      group_by(variable) %>%
      summarize(n_below = sum(below),
                below = mean(below),
              above = mean(above),
              coverage = mean(covered),
              ci_width = median(`97.5%` - `2.5%`),
              below0.1 = mean(below0.1),
              above0.9 = mean(above0.9),
              coverage0.8 = mean(covered0.8),
              ci_width0.8 = median(`90%` - `10%`),
              median_error = median(error),
              mean_squared_error = mean(error^2),
              median_abs_error = median(abs(error)))
  })) %>%
  select(-fit) %>%
  unnest(cols = c(validation)) %>%
  arrange(cutoff_year, model, scale_global)

validation_table <- validation_results_summary %>%
  #select(normal_data_model, data_model_df, scale_global, num_knots, cutoff_year, variable, model, scale_global, below, coverage, above, ci_width, below0.1, coverage0.8, above0.9, ci_width0.8, median_error, median_abs_error) %>%
  select(scale_global, num_knots, cutoff_year, variable, model, scale_global, below, coverage, above, ci_width, below0.1, coverage0.8, above0.9, ci_width0.8, median_error, median_abs_error) %>%
  mutate_at(vars(below, coverage, above, below0.1, above0.9, coverage0.8, ci_width, ci_width0.8, median_error, median_abs_error), signif, 3) %>%
  mutate_at(vars(below, coverage, above, below0.1, above0.9, coverage0.8), `*`, 100) %>%
  mutate_at(vars(below, coverage, above, below0.1, above0.9, coverage0.8), paste0, "%") %>%
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
         `% Below (0.8)` = below0.1,
         `% Included (0.8)` = coverage0.8,
         `% Above (0.8)` = above0.9,
         `CI Width (0.8)` = ci_width0.8,
         ME = median_error,
         MAE = median_abs_error
         ) %>%
  #mutate(Cutoff = ifelse(index %% 2 == 1, Cutoff, "")) %>%
  select(-index)

validation_table %>% View()

validation_table %>%
  mutate(variable = ifelse(variable == "eta", "crisis", "crisis-free")) %>%
  select(-scale_global, -num_knots, -ME, -MAE) %>%
  #select(scale_global, Cutoff, variable, Model, ME, MAE) %>%
  knitr::kable(format = "latex") 

validation_table %>%
  mutate(variable = ifelse(variable == "eta", "crisis", "crisis-free")) %>%
  select(Cutoff, variable, Model, ME, MAE) %>%
  #select(scale_global, Cutoff, variable, Model, ME, MAE) %>%
  knitr::kable(format = "latex") 


plot_life_transition <- function(fit) {
  fit$posteriors$transition_function_mean %>% 
    ggplot(aes(x = x, y = transition_function_mean)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
    geom_line() +
    scale_fill_brewer()
}

validations$fit[[2]]$posteriors$transition_function_mean %>%
  mutate(model = "Shock") %>%
  bind_rows(validations$fit[[1]]$posteriors$transition_function_mean %>%
              mutate(model = "No shocks")) %>%
  ggplot(aes(x = x, y = transition_function_mean, color = model)) +
  geom_line()

validations$fit[[1]] %>% plot_indicator("Suriname")

plot_life_transition(validations$fit[[1]]) + ylim(c(0, 10))
plot_life_transition(validations$fit[[3]]) + ylim(c(0, 10))
plot_life_transition(validations$fit[[6]])
plot_life_transition(validations$fit[[4]]) +
  geom_point(data = filter(e0_differences, year <= 2000), aes(x = e0 / 85, y = diff), alpha = 0.5, size = 0.1)

bind_rows(
  mutate(validations$fit[[1]]$posteriors$transition_function_mean, model = "No shocks"),
  mutate(validations$fit[[3]]$posteriors$transition_function_mean, model = "Shocks")
) %>%
  ggplot(aes(x = x, y = transition_function_mean, color = model)) +
  #geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  geom_line() +
  scale_fill_brewer()

all_fs <- bind_rows(
  mutate(validations$fit[[1]]$posteriors$transition_functions, model = "No shocks"),
  mutate(validations$fit[[3]]$posteriors$transition_functions, model = "Shocks")
) %>%
  filter(.width == 0.95)

all_fs %>%
  ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_pred, group = name)) +
  geom_line() +
  facet_wrap(~model)

all_fs %>%
  ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_pred, group = name)) +
  geom_line(aes(y = .upper)) +
  facet_wrap(~model)

all_fs %>%
  ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_pred, group = name)) +
  geom_line(aes(y = .lower)) +
  facet_wrap(~model)


country <- ""
bind_rows(
  mutate(validations$fit[[1]]$posteriors$transition_functions, model = "No shocks"),
  mutate(validations$fit[[4]]$posteriors$transition_functions, model = "With shocks")
) %>%
  filter(name == country, .width == 0.95) %>%
  ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_pred, color = model)) +
  #geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = filter(datM_diffs, name == country, year <= 2000), aes(x = e0, y = diff), color = "black") +
  geom_line() +
  geom_line(aes(y = .lower), lty = 2) +
  geom_line(aes(y = .upper), lty = 2) +
  #geom_vline(xintercept = c(45, 50, 55, 60, 65, 70), lty = 2) +
  scale_fill_brewer()

pit <- validations$fit[[1]]$samples$draws("pit") %>%
  spread_draws(pit[n])

gamma <- validations$fit[[1]]$samples$draws("gamma") %>%
  spread_draws(gamma[c, t]) %>%
  left_join(validations$fit[[1]]$country_index) %>%
  left_join(validations$fit[[1]]$time_index) %>%
  group_by(year, name, c, t) %>%
  median_qi()

gamma %>%
  ggplot(aes(x = year, y = gamma)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer() +
  geom_hline(yintercept = -8.66)

countries <- sort(validations$fit[[1]]$country_index$name)

plot_comparison <- function(fit1, fit2, name1, name2, country) {
  bind_rows(
    fit1$posteriors$temporal %>%
      filter(variable == "eta") %>%
      filter(name == country) %>%
      mutate(model = name1) %>%
      mutate_at(vars(`2.5%`, `50%`, `97.5%`), as.numeric) %>%
      select(year, model, `2.5%`, `50%`, `97.5%`),
    fit2$posteriors$temporal %>%
      filter(name == country) %>%
      filter(variable == "eta") %>%
      mutate(model = name2) %>%
      mutate_at(vars(`2.5%`, `50%`, `97.5%`), as.numeric) %>%
      select(year, model, `2.5%`, `50%`, `97.5%`),
  ) %>%
    ggplot(aes(x = year, y = `50%`, color = model)) +
    geom_line(aes(y = `2.5%`), lty = 2) +
    geom_line(aes(y = `97.5%`), lty = 2) +
    geom_line() +
    geom_point(aes(x = year, y = e0), color = "black", data = filter(fit1$data, name == country))
}

country <- "Myanmar"
pdf("plots/shocks_no_shocks_comparison.pdf", width = 6, height = 5)
for(country in countries) {
  print(country)
  p1 <- plot_comparison(validations$fit[[4]], validations$fit[[1]], "Shocks", "No shocks", country = country) + ggtitle("Cutoff = 2005")# + theme(legend.position = "none")
  #p2 <- plot_comparison(validations$fit[[2]], validations$fit[[5]], "Shocks", "No shocks", country = country) + ggtitle("Cutoff = 2010") + theme(legend.position = "none")
  #p3 <- plot_comparison(validations$fit[[3]], validations$fit[[6]], "Shocks", "No shocks", country = country) + ggtitle("Cutoff = 2015")
  #p <- (p1 + p2 + p3) & plot_annotation(title = country)
  p <- p1 + plot_annotation(title = country)
  print(p)
}
dev.off()


shock <- validations$fit[[4]]$samples %>% spread_draws(shock[c, t]) %>% group_by(c, t) %>% median_qi(.width = c(0.8, 0.9, 0.95))

shock %>%
  filter(.width == 0.95) %>%
  left_join(validations$fit[[3]]$country_index) %>%
  left_join(validations$fit[[3]]$time_index) %>%
  arrange(shock)

shock %>%
  left_join(validations$fit[[4]]$country_index) %>%
  left_join(validations$fit[[4]]$time_index) %>%
  filter(name == country) %>%
  ggplot(aes(x = year, y = shock)) + 
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer()

validation_results[c(1,2),] %>%
  unnest(validation) %>%
  group_by(model, name) %>% 
  filter(year == max(year), variable == "eta") %>%
  filter(covered == FALSE) %>%
  pull(name)

validation_results[c(1,2),] %>%
  unnest(validation) %>%
  group_by(model, name) %>% 
  filter(year == max(year), variable == "eta", model == "shock2") %>%
  arrange(error) %>%
  select(name, model, error)

validation_results[c(1,4),] %>%
  unnest(validation) %>%
  group_by(model, name) %>% 
  filter(year == max(year), variable == "eta") %>%
  select(model, name, `50%`) %>%
  pivot_wider(names_from = "model", values_from = "50%") %>%
  mutate(diff = spline - shock2) %>% 
  arrange(-diff)
