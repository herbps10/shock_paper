library(tidyverse)
library(wpp2024)
library(tidybayes)
library(bayesLife)
library(patchwork)

source("R/plot_theme.R")
source("R/lifeplus.R")
source("R/process_lifeplus.R")

data(UNlocations, package = "wpp2024")
data(e0M1, package = "wpp2024")
data(pop1, package = "wpp2024")
data(include_2010, package = "bayesLife")


# Number of countries before filtering
e0M1 |>
  filter(country_code < 900) |>
  distinct(name) |>
  nrow()

included_codes <- include_2010 |> filter(include_code %in% 1:2) |> pull(country_code)

large_countries <- pop1 |> 
  filter(`2023` >= 1e3) |>
  pull(country_code)

datM <- e0M1 |>
  filter(country_code %in% large_countries) |>
  pivot_longer(cols = `1950`:`2023`, names_to = "period", values_to = "e0") |>
  filter(country_code %in% included_codes) |>
  mutate(year = parse_integer(str_sub(period, 1, 4)),
         source = "WPP2024")

datM_diffs <- datM |>
  group_by(name) |>
  mutate(diff = c(NA, diff(e0)))

datM_diffs |>
  filter(e0 > 55) |>
  ggplot(aes(x = e0, y = diff)) +
  geom_point() +
  geom_smooth()

datM |>
  group_by(name) |>
  mutate(diff = c(NA, diff(e0)))  |>
  group_by(period) |>
  summarize(mean = mean(diff))

# Number of countries after filtering
datM |>
  distinct(name) |>
  nrow()

# Histogram of changes in e0
e0_differences <- datM |>
  group_by(name) |>
  mutate(diff = c(NA, diff(e0)))

e0_differences |>
  ggplot(aes(x = diff)) +
  geom_histogram(color = "white", binwidth = 1) +
  geom_boxplot(aes(y = -20), width = 30, alpha = 0.5) +
  labs(x = expression(paste("Difference in ", e[0], ": ", eta[ct] - eta[ct-1])), y = "Count")

countries <- c("Republic of Korea", "Bosnia and Herzegovina", "Cambodia", "Lebanon", "Timor-Leste", "Syrian Arab Republic", "Switzerland", "Norway")

datM |>   
  filter(name %in% countries) |>
  ggplot(aes(x = year + 2.5, y = e0)) +
  geom_point(size = 0.25) +
  geom_line(alpha = 0.5) +
  facet_wrap(~name) +
  #pub_theme +
  labs(x = "Year", y = expression(e[0]))

ggsave("plots/life_examples.pdf", height = 5, width = 10)

# What is the 2*SD(eps) threshold?
threshold <- 2 * fits$fit[[2]]$samples$summary("epsilon_scale")$median

# How many of the observed differences fall below this threshold?
mean(e0_differences$diff < -threshold, na.rm = TRUE)

set.seed(4)
random_countries <- sample(unique(datM$name), 20)
random_countries <- unique(c(random_countries, c("Republic of Korea", "Dem. People's Republic of Korea", "Bangladesh", "Lebanon", "Somalia")))
random_countries <- c("Republic of Korea", "Dem. People's Republic of Korea", "Bangladesh", "Lebanon", "Somalia")

fits <- expand_grid(
  #scale_global = c(1e-3, 1e-2, 1e-1),
  scale_global = 1e-2,
  #model = c("shock2")
  #model = "logistic_shock",
  #outlier_threshold = 1e3
  model = "logistic",
  outlier_threshold = 5
) |>
  bind_rows(
    #tibble(scale_global = 1e-2, model = "logistic", outlier_threshold = 1e3),
    #tibble(scale_global = 1e-2, model = "logistic_shock", outlier_threshold = 1e3)
    tibble(scale_global = 1e-2, model = "logistic", outlier_threshold = 5)
  ) |>
  mutate(fit = pmap(list(scale_global, model, outlier_threshold), function(scale_global, model, outlier_threshold) {
    lifeplus(
      datM |> filter(name %in% random_countries),
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,
      
      hierarchical = TRUE,
      centered = TRUE,
      
      outlier_threshold = outlier_threshold,
      
      model = model,
      
      adapt_delta = 0.95,
      max_treedepth = 12,
      parallel_chains = 4,
      iter_warmup = 250,
      #iter_sampling = 1e3,
      iter_sampling = 250,
      
      extra_stan_data = list(
        scale_global = scale_global,
        slab_scale = 10,
        slab_df = 6
      )
    )
  }))

fit_shock <- fits$fit[[1]]

fit <- fits$fit[[1]]
tidybayes::spread_draws(fit$samples, Delta1[c]) |> 
  median_qi(.width =c(0.5, 0.9, 0.95)) |> 
  left_join(fit$country_index) |> 
  ggplot(aes(x = Delta4, y = reorder(name, sort(name)))) +
  tidybayes::geom_interval(aes(xmin = .lower, xmax = .upper)) + 
  geom_point() + 
  scale_color_brewer()

bayesplot::mcmc_dens(fits$fit[[1]]$samples$draws(c("sigma_Delta1", "sigma_Delta2", "sigma_Delta3", "sigma_Delta4", "sigma_z", "sigma_k", "epsilon_scale")))
bayesplot::mcmc_dens(fits$fit[[2]]$samples$draws(c("mu_Delta1", "mu_Delta2", "mu_Delta3", "mu_Delta4", "mu_z", "mu_k", "epsilon_scale")))

fit_shock   <- fits$fit[[1]]
fit_noshock <- fits$fit[[2]]

left_join(
  fit_shock$posteriors$transition_functions |> filter(.width == 0.5) |> select(name, x, transition_function_pred),
  fit_noshock$posteriors$transition_functions |> filter(.width == 0.5) |> select(name, x, transition_function_pred),
  by = c("name", "x")
) |>
  group_by(name) |>
  summarize(diff = sum(abs(transition_function_pred.x - transition_function_pred.y))) |>
  arrange(diff)
  

Delta_shock <- gather_draws(fit_shock$samples$draws(c("Delta1", "Delta2", "Delta3", "Delta4", "k", "z")), Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c])
Delta_noshock <- gather_draws(fit_noshock$samples$draws(c("Delta1", "Delta2", "Delta3", "Delta4", "k", "z")), Delta1[c], Delta2[c], Delta3[c], Delta4[c], k[c], z[c])

left_join(
  Delta_shock  |>
    group_by(c, .variable) |>
    median_qi() |>
    mutate(model = "shock") |>
    select(c, .variable, .value, .lower, .upper, model),
  Delta_noshock  |>
    group_by(c, .variable) |>
    median_qi() |>
    mutate(model = "no shock") |>
    select(c, .variable, .value, .lower, .upper, model)
  , by = c(".variable", "c")
) |>
  left_join(fit_shock$country_index) |>
  ggplot(aes(y = name)) +
  geom_pointinterval(aes(x = .value.x, xmin = .lower.x, xmax = .upper.x, color = "Shock")) +
  geom_pointinterval(aes(x = .value.y, xmin = .lower.y, xmax = .upper.y, color = "No Shock")) +
  scale_fill_brewer() +
  facet_wrap(~.variable, scales = "free_x")

(Delta_shock |> median_qi()) |> left_join(Delta_noshock |> median_qi(), by = "c") |> left_join(fit_shock$country_index)  |>
  select(name, Delta1.x, Delta1.y, Delta2.x, Delta2.y, Delta3.x, Delta3.y, Delta4.x, Delta4.y, k.x, k.y, z.x, z.y)

name <- "Afghanistan"
name <- random_countries
plot_transition(fit_noshock_naive, name)
plot_transition(fit_noshock) + ylim(c(0, 10))
plot_transition(fit_shock) + ylim(c(0, 10))
plot_shock(fit_shock)

plot_temporal("eta", fit_noshock, plot_data = TRUE) + ylim(c(15, 150))
plot_temporal("eta_crisisfree", fit_shock, plot_data = TRUE) + ylim(c(15, 150))
plot_temporal("eta", fit_shock, plot_data = TRUE) + ylim(c(15, 150))

fit_shock$posteriors$temporal |> filter(year == 2100) |> mutate(ci_width = `99.9%` - `0.1%`)

comp <- left_join(
  fit_shock$posteriors$temporal |> filter(year == 2100, variable == "eta") |> mutate(ci_width = `90%` - `10%`) |> select(name, `50%`, ci_width),
  fit_noshock$posteriors$temporal |> filter(year == 2100, variable == "eta") |> mutate(ci_width = `90%` - `10%`) |> select(name, `50%`, ci_width)
 , by = c("name"))

comp |> filter(ci_width.x > ci_width.y)

comp |>
  ggplot(aes(x = `ci_width.y`, y = `ci_width.x`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

comp |>
  ggplot(aes(x = `50%.x`, y = `50%.y`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

bind_rows(
  fit_noshock$posteriors$temporal |> filter(year == 2100) |> mutate(ci_width = `90%` - `10%`) |>
    mutate(model = "no shocks"),
  fits$fit[[2]]$posteriors$temporal |> filter(year == 2100) |> mutate(ci_width = `90%` - `10%`) |>
    mutate(model = "shocks")
) 
