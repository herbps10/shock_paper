library(tidyverse)
library(wpp2024)
library(tidybayes)
library(bayesLife)
library(patchwork)

source("R/plot.R")
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

#ggsave("plots/life_examples.pdf", height = 5, width = 10)

# What is the 2*SD(eps) threshold?
threshold <- 2 * fits$fit[[2]]$samples$summary("epsilon_scale")$median

# How many of the observed differences fall below this threshold?
mean(e0_differences$diff < -threshold, na.rm = TRUE)

set.seed(7)
random_countries <- unique(c(sample(unique(datM$name), 50)))

fits <- expand_grid(
  scale_global = c(1e-1),
  transition = c("logistic"),
  shock = c(TRUE),
  data_model = c("normal"),
  hierarchical = c(0),
  include_prior = c(0)
) |>
  #mutate(include_prior = ifelse(transition == "gp", 1, 0)) |>
  filter(!(data_model == "mixture" & shock == TRUE)) |>
  mutate(fit = pmap(list(transition, shock, data_model, scale_global, hierarchical, include_prior), \(transition, shock, data_model, scale_global, hierarchical, include_prior) {
    lifeplus(
      datM,
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2050,
      
      transition = transition,
      shock = shock,
      data_model = data_model,
      
      spline_degree = 2,
      num_knots = 7, 
      
      hierarchical = hierarchical,
      centered = FALSE,
      
      outlier_threshold = 5,
      
      adapt_delta = 0.90,
      max_treedepth = 12,
      parallel_chains = 4,
      iter_warmup = 250,
      #iter_sampling = 1e3,
      iter_sampling = 500,
      
      epsilon_prior = c(0, 1),
      
      extra_stan_data = list(
        scale_global = scale_global,
        slab_scale = 10,
        slab_df = 6,
        L = 1.5,
        M = 25,
        heteroskedastic = 0,
        include_prior = include_prior
      )
    )
  }))

fits$fit[[1]]$samples$summary("epsilon_sigma")

plot_shock(fits$fit[[1]], "Bosnia and Herzegovina")

plot_temporal("eta", fits$fit[[1]], "Bosnia and Herzegovina")
plot_transition(fits$fit[[1]], "Somalia") + coord_cartesian(ylim = c(0, 2.5))

np <- nuts_params(fits$fit[[1]]$samples)

bayesplot::mcmc_pairs(
  fits$fit[[1]]$samples$draws(c("epsilon_sigma", "c_slab", "lambda_tilde_sd")), np = np
)

bayesplot::mcmc_pairs(
  fits$fit[[1]]$samples$draws(c("Delta[1,1]", "Delta[1,2]", "Delta[1,3]", "Delta[1,4]", "Delta[1,5]", "Delta[1,6]", "epsilon_sigma")), np = np
)

fits$fit[[1]]$samples$draws(c("lambda_tilde")) |>
  spread_draws(lambda_tilde[k]) |>
  median_qi() |>
  ggplot(aes(x = k, y = lambda_tilde)) +
  geom_point()


fits$fit[[1]]$posteriors$transition_params_corr |>
  ggplot(aes(x = k1, y = k2)) +
  geom_tile(aes(fill = median)) +
  geom_text(aes(label = scales::number(median, accuracy = 0.01, scale = 1)), color = "white") +
  coord_fixed() +
  labs(x = "Delta", y = "Delta")

fits$fit[[1]]$posteriors$transition_params |>
  ggplot(aes(x = `50%`, y = name)) +
  geom_point() +
  facet_wrap(~variable, scales = "free_x")

plot_transition(fits$fit[[1]], "Somalia")
plot_mean_transition(fits$fit[[1]])

fits$fit[[1]]$samples$draws("shock2")

fits$fit[[1]]$posteriors$temporal |>
  filter(variable == "shock2") |>
  ggplot(aes(x = year, y = `50%`)) +
  geom_line(aes(y = `0.1%`)) +
  geom_line(aes(y = `99.9%`))

fits$fit[[1]]$posteriors$temporal |>
  filter(variable %in% c("eta", "eta_shockfree")) |>
  ggplot(aes(x = year, y = `50%`, color = variable)) +
  geom_line(aes(y = `0.1%`)) +
  geom_line(aes(y = `99.9%`))
