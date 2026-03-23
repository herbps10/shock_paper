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

set.seed(4)
random_countries <- unique(c(sample(unique(datM$name), 25))) #, c("Republic of Korea", "Bosnia and Herzegovina", "Cambodia", "Lebanon", "Timor-Leste", "Syrian Arab Republic", "Switzerland", "Norway")))

random_countries <- c("Oman")

fits <- expand_grid(
  scale_global = 1e-1,
  transition = c("logistic"),
  shock = c(FALSE),
  data_model = c("outlier"),
  hierarchical = c(0, 1),
) |>
  mutate(include_prior = ifelse(transition == "gp", 1, 0)) |>
  filter(!(data_model == "mixture" & shock == TRUE)) |>
  mutate(fit = pmap(list(transition, shock, data_model, scale_global, hierarchical, include_prior), \(transition, shock, data_model, scale_global, hierarchical, include_prior) {
    lifeplus(
      datM |> filter(name %in% random_countries),
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,
      
      transition = transition,
      shock = shock,
      data_model = data_model,
      
      spline_degree = 2,
      num_knots = 7, 
      
      hierarchical = hierarchical,
      centered = FALSE,
      
      outlier_threshold = 5,
      
      adapt_delta = 0.95,
      max_treedepth = 12,
      parallel_chains = 4,
      iter_warmup = 5e2,
      #iter_sampling = 1e3,
      iter_sampling = 5e2,
      
      extra_stan_data = list(
        scale_global = scale_global,
        slab_scale = 10,
        slab_df = 6,
        L = 2,
        M = 10,
        heteroskedastic = 0,
        include_prior = include_prior
      )
    )
  }))

fits$fit[[1]]$samples$summary(c("gamma", "inner_epsilon_sigma", "outer_epsilon_sigma"))
fits$fit[[2]]$samples$summary(c("gamma", "inner_epsilon_sigma", "outer_epsilon_sigma"))

ggplot(fits$fit[[1]]$posteriors$transition_params, aes(y = name, x = `50%`)) +
  geom_point(aes(color = "independent")) +
  geom_point(aes(color = "hierarchical"), data =fits$fit[[2]]$posteriors$transition_params) +
  facet_wrap(~variable, scales = "free_x")

likelihood_ratios <- spread_draws(fits$fit[[1]]$samples$draws("likelihood_ratio"), likelihood_ratio[c, t]) |>
  left_join(fits$fit[[1]]$country_index) |>
  left_join(fits$fit[[1]]$time_index) |>
  group_by(name, year) |>
  median_qi(likelihood_ratio, .width = c(0.5))

likelihood_ratios |>
  arrange(-likelihood_ratio)

likelihood_ratios |>
  filter(name %in% c("Oman")) |>
  ggplot(aes(x = year, y = likelihood_ratio)) +
  geom_point() +
  facet_wrap(~name, scales = "free_y")

plot_temporal("eta", fits$fit[[1]], "Somalia")
plot_temporal("eta", fits$fit[[2]], "Somalia")

plot_mean_transition(fits$fit[[1]])
plot_mean_transition(fits$fit[[2]])

plot_transition(fits$fit[[1]], "Somalia")
plot_transition(fits$fit[[2]], "Somalia")

comp <- left_join(
  fits$fit[[1]]$posteriors$temporal |> filter(year == 2100) |> select(name, `10%`, `50%`, `90%`),
  fits$fit[[2]]$posteriors$temporal |> filter(year == 2100) |> select(name, `10%`, `50%`, `90%`),
  by = "name"
)  |>
  mutate(ci_width.x = `90%.x` - `10%.x`, ci_width.y = `90%.y` - `10%.y`)

comp |>
  filter((ci_width.x / ci_width.y) > 3)

comp |> #ggplot(aes(`50%.x`, `50%.y`)) +
  ggplot(aes(ci_width.x, ci_width.y)) +
  geom_point() +
  labs(x = "independent", y = "hierarchical") +
  geom_abline()

fits$fit[[1]]$posteriors$temporal

plot_temporal("eta", fits$fit[[1]])
plot_temporal("eta_shockfree", fits$fit[[1]])
plot_temporal("eta", fits$fit[[2]])
plot_temporal("eta_shockfree", fits$fit[[1]], plot_data = TRUE)

plot_transition(fits$fit[[1]]) + coord_cartesian(xlim = c(80, 110), ylim = c(0, 1))

plot_mean_transition(fits$fit[[1]])

fits$p <- map(arrange(fits, transition, data_model) |> pull(fit), \(x) plot_transition(x) + theme(legend.position = "none") + ggtitle(paste(x$transition, x$data_model, ifelse(x$shock, "shock term", ""))))
gridExtra::grid.arrange(grobs = fits$p, ncol = 3)

fits$p <- map(arrange(fits, transition, data_model) |> pull(fit), \(x) plot_temporal("eta", x) + theme(legend.position = "none") + ggtitle(paste(x$transition, x$data_model, ifelse(x$shock, "shock term", ""))))
gridExtra::grid.arrange(grobs = fits$p, ncol = 3)

plot_temporal("eta", fits$fit[[1]])
plot_temporal("eta_shockfree", fits$fit[[2]])


plot_shock(fits$fit[[2]])

fits$fit[[1]]$samples$summary("epsilon_sigma")
fits$fit[[2]]$samples$summary("epsilon_sigma")

shinystan::launch_shinystan(fits$fit[[1]]$samples)

spread_draws(fits$fit[[1]]$samples$draws("epsilon_sigma_pred"), epsilon_sigma_pred[i]) |>
  left_join(tibble(i = 1:length(fits$fit[[1]]$stan_data$grid), x = fits$fit[[1]]$stan_data$grid)) |>
  group_by(x) |>
  median_qi(epsilon_sigma_pred, .width = c(0.8, 0.9, 0.95)) |>
  ggplot(aes(x = x * 110, y = epsilon_sigma_pred)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer()
 
fits$fit[[1]]$samples$summary("epsilon_sd")
fits$fit[[2]]$samples$summary("epsilon_sd")
fits$fit[[3]]$samples$summary("epsilon_sd")

plot_shock(fits$fit[[1]])
plot_shock(fits$fit[[2]])

plot_shock(fits$fit[[3]])

plot_transition(fits$fit[[1]], "intercept")

plot_transition(fits$fit[[2]])
plot_transition(fits$fit[[3]])

plot_temporal("eta", fits$fit[[1]])
plot_temporal("eta", fits$fit[[2]])
plot_temporal("eta", fits$fit[[3]])

plot_temporal("eta_crisisfree", fits$fit[[1]])
plot_temporal("eta_crisisfree", fits$fit[[2]])
plot_temporal("eta_crisisfree", fits$fit[[3]])

fits$fit[[2]]$posteriors$temporal |>
  filter(year == 2100) |>
  mutate(ci_width = `90%` - `10%`) |>
  select(variable, name, ci_width) |>
  pivot_wider(names_from = "variable", values_from = "ci_width")
 
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
fit_noshock_naive <- fits$fit[[3]]

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

diagnostics <- fits |>
  filter(model == "logistic_shock") |>
  mutate(diagnostics = map(fit, \(fit) fit$samples$diagnostic_summary())) |>
  select(country, diagnostics)

diagnostics |>
  mutate(treedepth = map_int(diagnostics, \(x) max(x$num_max_treedepth))) |>
  arrange(-treedepth)

fit_shock <- fits |>
  filter(model == "logistic_shock", country == "Kazakhstan") |>
  pull(fit)

fit_noshock <- fits |>
  filter(model == "logistic", country == "Kazakhstan") |>
  pull(fit)

fits |>
  select(model, country, params) |>
  unnest(params) |>
  select(model, country, variable, q95) |>
  pivot_wider(names_from = c("model"), values_from = "q95") |>
  mutate(diff = abs(logistic_shock - logistic)) |>
  arrange(-diff) |>
  filter(country == "Kazakhstan")

plot_transition(fit_shock[[1]]) + geom_hline(yintercept = 1.15) + ylim(c(-2.5, 7.5))
plot_transition(fit_noshock[[1]]) + geom_hline(yintercept = 1.15) + ylim(c(-2.5, 7.5))
plot_shock(fit_shock[[1]])
plot_temporal("eta", fit_shock[[1]])
plot_temporal("eta", fit_noshock[[1]])
plot_temporal("eta_crisisfree", fit_shock[[1]])

fit_shock[[1]]$samples$summary("epsilon_variance")
fit_noshock[[1]]$samples$summary("epsilon_variance")

fits$params <- map(fits$fit, \(x) x$samples$summary(c("Delta1", "Delta2", "Delta3", "Delta4", "k", "z", "epsilon_variance")))
fits |>
  select(model, country, params) |>
  unnest(params) |>
  select(model, country, variable, mean) |>
  pivot_wider(names_from = c("model"), values_from = "mean") |>
  mutate(diff = abs(logistic_shock - logistic)) |>
  arrange(-diff)

fits$pred <- map(fits$fit, \(x) x$posteriors$temporal |> filter(variable %in% c("eta", "eta_crisisfree"), year == 2100))

fits |>
  select(model, country, pred) |>
  unnest(pred) |>
  mutate(ci_width = `90%` - `10%`) |>
  filter((model == "logistic_shock" & variable == "eta_crisisfree") | (model == "logistic" & variable == "eta")) |>
  select(model, country, ci_width) |>
  pivot_wider(names_from = "model", values_from = "ci_width") |>
  ggplot(aes(x = logistic, y = logistic_shock)) +
  geom_point() +
  geom_abline()
 
name <- "Niger"
plot_transition(fit_noshock_naive, name)
plot_transition(fit_noshock, name)
plot_transition(fit_shock, name)
plot_shock(fit_shock, "Denmark")

plot_temporal("eta", fit_noshock, name, plot_data = TRUE) + ylim(c(15, 150))
plot_temporal("eta_crisisfree", fit_shock, name, plot_data = TRUE) + ylim(c(15, 150))
plot_temporal("eta", fit_shock, name, plot_data = TRUE) + ylim(c(15, 150))

plot_shock(fits$fit[[1]], "Niger")

plot_transition(fits$fit[[1]], "Lebanon") + ylim(c(0, 5))
plot_transition(fits$fit[[2]], "Lebanon") + ylim(c(0, 5))

fit_shock$posteriors$temporal |> filter(year == 2100) |> mutate(ci_width = `99.9%` - `0.1%`)

comp <- left_join(
  fit_shock$posteriors$temporal |> filter(year == 2100, variable == "eta_crisisfree") |> mutate(ci_width = `90%` - `10%`) |> select(name, `50%`, ci_width),
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
