library(tidyverse)
library(patchwork)
library(invgamma)
devtools::load_all("../BayesTransitionModels/")

source("R/plot_theme.R")

fit <- read_rds("./fit")

fit1 <- read_rds("./fit_shock_1e2")
fit1$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "*2c9b3e.csv")))

fit2 <- read_rds("./fit_shock_1e4")
fit2$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "*41f83f.csv")))

fit3 <- read_rds("./fit_shock_1e6")
fit3$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "*44a5ba.csv")))

fit$data %>%
  filter(name_country %in% c("Timor-Leste")) %>%
  rename(`Data Source` = data_series_type) %>%
  ggplot(aes(x = year, y = contraceptive_use_modern, shape = `Data Source`)) +
  geom_point() +
  facet_wrap(~name_country) +
  pub_theme +
  labs(x = "Year", y = "mCPR")

ggsave("plots/mcpr_example.pdf", height = 4, width = 6)


#
# Prior/Posterior plots for shock hyperparameters
#

caux_priors <- expand_grid(
  x = 10^seq(-4, 2, 0.01),
  s = c(0.1, 0.5, 1),
  nu = c(1, 3, 100)
) %>%
  mutate(
    y = dinvgamma(x, shape = 0.5 * nu, rate = 0.5 * s^2 * nu)
  ) 

ggplot(caux_priors, aes(x = x, y = y, color = factor(nu))) +
  geom_line() +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = expression(s), breaks = NULL, labels = NULL)) +
  facet_wrap(~s, scales = "free") +
  labs(x = expression(c^2), y = "prior density") +
  scale_color_discrete(name = expression(nu)) +
  pub_theme

ggsave("plots/prior_density_csq.pdf", width = 8, height = 3)

caux_prior <- tibble(
  x = 10^seq(-1, 6, 0.1),
  y = dinvgamma(x, shape = 0.5, rate = 0.5 * 1)
)

ggplot(caux_prior, aes(x, y)) +
  geom_line() +
  scale_x_log10()

caux_posts <- map2(
  list(fit1, fit2, fit3), c("1e-2", "1e-4", "1e-6"), 
  function(fit, model) {
    fit$samples %>% spread_draws(caux) %>%
      mutate(model = model)
  }) %>%
  bind_rows()

caux_posts %>%
  ggplot(aes(x = caux)) +
  geom_density(aes(color = "Posterior KDE")) +
  geom_line(data = caux_prior, aes(x = x, y, color = "Prior")) +
  #scale_x_sqrt(sec.axis = sec_axis(~ ., name = expression(tau[0]), breaks = NULL, labels = NULL)) +
  scale_x_log10() +
  scale_color_discrete(name = "") +
  labs(x = expression(c), y = "density") +
  facet_wrap(~model) +
  pub_theme

ggsave("plots/mcpr_csq_prior_posterior.pdf", width = 8, height = 3)

global_shrinkage_posts <- map2(list(fit1, fit2, fit3), c("1e-2", "1e-4", "1e-6"),
  function(fit, model) {
    fit$samples %>% 
      spread_draws(global_shrinkage) %>%
        mutate(model = model)
  }
) %>%
  bind_rows()

global_shrinkage_prior <- expand_grid(
  x = 10^seq(-8, 1, 0.1),
  tau0 = c(1e-2, 1e-4, 1e-6)
) %>%
  mutate(
    y = LaplacesDemon::dhalft(x, scale = tau0, nu = 1),
    model = case_when(
      tau0 == 1e-2 ~ "tau0 = 1e-2",
      tau0 == 1e-4 ~ "tau0 = 1e-4",
      tau0 == 1e-6 ~ "tau0 = 1e-6"
    )
  )

global_shrinkage_posts %>%
  ggplot(aes(x = global_shrinkage)) +
  geom_density(aes(color = model), ) +
  #geom_line(data = global_shrinkage_prior, aes(x, y, color = "Prior")) +
  scale_x_log10() +
  scale_color_discrete(name = expression(tau[0])) +
  labs(y = "Posterior kernel density\nestimator", x = expression(tau)) +
  pub_theme

ggsave("plots/mcpr_tau_posteriors.pdf", width = 8, height = 3)

global_shrinkage_prior %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(~model, scales = "free_y") +
  scale_x_log10()


shocks_raw1 <- spread_draws(fit1$samples, shock[c,t])
shocks_raw2 <- spread_draws(fit2$samples, shock[c,t]) 
shocks_raw3 <- spread_draws(fit3$samples, shock[c,t]) 

mean(shocks_raw1$shock < -0.401) * 100
mean(shocks_raw2$shock < -0.401) * 100
mean(shocks_raw3$shock < -0.401) * 100

# Find most different estimates
fit1$posteriors$temporal %>%
  filter(variable == "eta") %>%
  select(year, name_country, `50%`) %>%
  left_join(
    fit2$posteriors$temporal %>%
      filter(variable == "eta") %>%
      select(year, name_country, `50%`),
    by = c("year", "name_country"), 
    suffix = c(".1e2", ".1e4")
  ) %>%
  mutate(diff = abs(`50%.1e2` - `50%.1e4`)) %>%
  arrange(-diff)

plot_

country <- "Rwanda"
p1 <- plot_indicator(fit,  country) + labs(subtitle = "No shocks") + pub_theme
p2 <- plot_indicator(fit1, country) + labs(subtitle = expression(Shocks~(tau[0]==1e-2))) + pub_theme

country <- "Timor-Leste"
p3 <- plot_indicator(fit,  country) + pub_theme
p4 <- plot_indicator(fit1, country) + pub_theme

(p1 + p2) / (p3 + p4) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") +
  pub_theme
ggsave("plots/mcpr_shock_example.pdf", width = 8, height = 5)



country <- "Rwanda"
p1 <- plot_indicator(fit,  country) + ggtitle("No shocks")
p2 <- plot_indicator(fit1, country) + ggtitle(expression(tau[0]==1e-2))
p3 <- plot_indicator(fit2, country) + ggtitle(expression(tau[0]==1e-4))
p4 <- plot_indicator(fit3, country) + ggtitle(expression(tau[0]==1e-6))
(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")
ggsave("plots/rwanda_example.pdf")

country <- "Timor-Leste"
p1 <- plot_indicator(fit,  country) + ggtitle("No shocks")
p2 <- plot_indicator(fit1, country) + ggtitle(expression(tau[0]==1e-2))
p3 <- plot_indicator(fit2, country) + ggtitle(expression(tau[0]==1e-4))
p4 <- plot_indicator(fit3, country) + ggtitle(expression(tau[0]==1e-6))
(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")
ggsave("plots/timor_leste_example.pdf")

country <- "Somalia"
p1 <- plot_indicator(fit,  country) + ggtitle("No shocks")
p2 <- plot_indicator(fit1, country) + ggtitle(expression(tau[0]==1e-2))
p3 <- plot_indicator(fit2, country) + ggtitle(expression(tau[0]==1e-4))
p4 <- plot_indicator(fit3, country) + ggtitle(expression(tau[0]==1e-6))
(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")
ggsave("plots/somalia_example.pdf")


projection_comparison <- fit$posteriors$temporal %>%
  filter(variable == "eta", year == 2030) %>%
  mutate(median_no_shocks = `50%`, ci_width_no_shocks = `97.5%` - `2.5%`) %>%
  select(name_country, median_no_shocks, ci_width_no_shocks) %>%
  left_join(
    fit1$posteriors$temporal %>%
      filter(variable == "eta", year == 2030) %>%
      mutate(median_shocks = `50%`, ci_width_shocks = `97.5%` - `2.5%`) %>%
      select(name_country, median_shocks, ci_width_shocks)
  ) 


pm10 <- tribble(
  ~x, ~y,
  -1, -1 - 0.05,
  2, 2 - 0.05,
  2, 2 + 0.05,
  -1, -1 + 0.05
)

pm5 <- tribble(
  ~x, ~y,
  -1, -1 - 0.025,
  2, 2 - 0.025,
  2, 2 + 0.025,
  -1, -1 + 0.025
)

p1 <- projection_comparison %>%
  ggplot(aes(x = median_no_shocks, median_shocks)) +
  geom_point() +
  geom_polygon(data = pm10, aes(x, y), alpha = 0.2) +
  geom_polygon(data = pm5, aes(x, y), alpha = 0.2) +
  coord_cartesian(xlim = range(projection_comparison$median_no_shocks),
                  ylim = range(projection_comparison$median_shocks)) +
  labs(subtitle = "posterior median", x = "no shocks", y = expression(shocks~(tau[0]==0.01))) +
  geom_abline(slope = 1, lty = 2) +
  pub_theme


p2 <- projection_comparison %>%
  ggplot(aes(x = ci_width_no_shocks, y = ci_width_shocks)) +
  geom_point() +
  geom_abline(slope = 1, lty = 2) +
  geom_polygon(data = pm10, aes(x, y), alpha = 0.2) +
  geom_polygon(data = pm5, aes(x, y), alpha = 0.2) +
  coord_cartesian(xlim = range(projection_comparison$ci_width_no_shocks),
                  ylim = range(projection_comparison$ci_width_shocks)) +
  labs(subtitle = "95% credible interval width", x = "no shocks", y = expression(shocks~(tau[0]==0.01)), lty = 2) +
  pub_theme

(p1 + p2) + plot_annotation(title = "mCPR by country, 2030")
ggsave("plots/ci_width_comparison.pdf", width = 10, height = 5)

shocks1 <- shocks_raw1 %>%
  median_qi(shock) %>%
  left_join(fit2$country_index) %>%
  left_join(fit2$time_index)

shocks2 <- shocks_raw2 %>%
  median_qi(shock) %>%
  left_join(fit1$country_index) %>%
  left_join(fit1$time_index)

shocks3 <- shocks_raw3 %>%
  median_qi(shock) %>%
  left_join(fit1$country_index) %>%
  left_join(fit1$time_index)

largest_shocks1 <- shocks1 %>%
  arrange(.lower) %>% filter(.lower < -0.401) %>% 
  distinct(name_country)

largest_shocks2 <- shocks2 %>%
  arrange(.lower) %>% filter(.lower < -0.401) %>% 
  distinct(name_country)

largest_shocks3 <- shocks3 %>%
  arrange(.lower) %>% filter(.lower < -0.401) %>% 
  distinct(name_country)

largest_shocks1
largest_shocks2
largest_shocks3

plot_indicator(fit2, largest_shocks$name_country)

plot_indicator(fit1, "Rwanda")
ggsave("plots/rwanda_0.01.pdf", width = 8, height = 5)

plot_indicator(fit1, "Timor-Leste")
ggsave("plots/timor_leste_0.01.pdf", width = 8, height = 5)

plot_indicator(fit0.001, "Rwanda")

plot_indicator(fit1, "Somalia")
plot_indicator(fit2, "Somalia")
plot_indicator(fit3, "Somalia")

rwanda1 <- shocks_raw1 %>%
  left_join(fit1$time_index) %>%
  left_join(fit1$country_index) %>%
  filter(name_country == "Rwanda", year %in% 1993:1999)

rwanda2 <- shocks_raw2 %>%
  left_join(fit1$time_index) %>%
  left_join(fit1$country_index) %>%
  filter(name_country == "Rwanda", year %in% 1993:1999)

rwanda3 <- shocks_raw3 %>%
  left_join(fit1$time_index) %>%
  left_join(fit1$country_index) %>%
  filter(name_country == "Rwanda", year %in% 1993:1999)

rwanda1 %>% mutate(model = "1e2") %>%
  bind_rows(rwanda2 %>% mutate(model = "1e4")) %>%
  bind_rows(rwanda3 %>% mutate(model = "1e6")) %>%
  ggplot(aes(x = -shock)) +
  geom_density(aes(color = model)) +
  scale_x_log10() +
  facet_wrap(~year)

shocks1 %>%
  filter(name_country == "Somalia") %>%
  ggplot(aes(x = year, y = shock)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer()

shocks2 %>%
  filter(name_country == "Somalia") %>%
  ggplot(aes(x = year, y = shock)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer()

shocks3 %>%
  filter(name_country == "Somalia") %>%
  ggplot(aes(x = year, y = shock)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer()

fit$posteriors$ar %>%
  pivot_longer(c(est_rho, est_tau)) %>%
  mutate(model = "No shocks") %>%
  bind_rows(fit1$posteriors$ar %>% pivot_longer(c(est_rho, est_tau)) %>% mutate(model = "tau_0 = 1e-2")) %>%
  #bind_rows(fit2$posteriors$ar %>% pivot_longer(c(est_rho, est_tau)) %>% mutate(model = "tau_0 = 1e-4")) %>%
  #bind_rows(fit3$posteriors$ar %>% pivot_longer(c(est_rho, est_tau)) %>% mutate(model = "tau_0 = 1e-6")) %>%
  ggplot(aes(x = value, color = model)) +
  geom_density() +
  facet_wrap(~name, scale = "free")

fit$posteriors$ar %>%
  mutate(model = "no shocks") %>%
  bind_rows(fit1$posteriors$ar %>% mutate(model = "shocks")) %>%
  mutate(unconditional_sd = est_tau / sqrt(1 - est_rho^2)) %>%
  ggplot(aes(x = unconditional_sd, color = model)) +
  geom_density() +
  scale_color_manual(values = scales::hue_pal()(2), labels = c("no shocks", expression(shocks~(tau[0]==0.01)))) +
  labs(x = "Unconditional standard deviation", y = "Kernel density estimate", subtitle = "mCPR AR(1) stochastic smoothing component ") +
  pub_theme +
  theme(legend.position = "top", legend.title = element_blank()) 

ggsave("plots/mcpr_ar_comparisons.pdf", width = 8, height = 6)



countries <- sort(unique(fit$data$name_country))
pdf("plots/mcpr_all_countries.pdf", width = 10, height = 4)
for(country in countries) {
  print(country)
  p1 <- plot_indicator(fit,  country) + labs(subtitle = "No shocks") + pub_theme
  p2 <- plot_indicator(fit1, country) + labs(subtitle = expression(Shocks~(tau[0]==1e-2))) + pub_theme
 
  p <- (p1 + p2) + plot_layout(guides = "collect") + pub_theme
  print(p)
}
dev.off()


compare_fits <- function(fits, areas, fit_names) {
  post <- map2(fits, fit_names, function(fit, fit_name) fit$posteriors$temporal %>%
    filter(name_country %in% areas, variable == "eta") %>%
    mutate(model = fit_name)) %>%
    bind_rows()
  
  data <- fits[[1]]$data
  
  filtered_data <- data %>%
    dplyr::filter(!!sym(fits[[1]]$area) %in% areas)
  
  if(!is.null(fits[[1]]$se)) {
    filtered_data <- filtered_data %>% mutate(
      lower = truncnorm::qtruncnorm(0.025, mean = !!sym(fits[[1]]$y), sd = !!sym(fits[[1]]$se), a = 0, b = 1),
      upper = truncnorm::qtruncnorm(0.975, mean = !!sym(fits[[1]]$y), sd = !!sym(fits[[1]]$se), a = 0, b = 1)
    )
  }
  
  ggplot(post, aes_string(x = fits[[1]]$year, y = "`50%`")) +
    geom_line(aes(y = .data$`2.5%`, color = model), lty = 2) +
    geom_line(aes(y = .data$`97.5%`, color = model), lty = 2) +
    geom_line(aes(color = model)) +
    ggplot2::geom_errorbar(aes_string(y = fits[[1]]$y, ymin = "lower", ymax = "upper"), alpha = 0.3, width = 0, data = filtered_data) +
    ggplot2::geom_point(aes_string(y = fits[[1]]$y), data = filtered_data, alpha = 0.7) +
    scale_fill_brewer(direction = -1) +
    facet_wrap(vars(!!sym(fits[[1]]$area))) +
    labs(x = fits[[1]]$year)
}

countries <- sort(unique(fit$data$name_country))
pdf("plots/mcpr_all_countries_comparison.pdf", width = 7, height = 4)
for(country in countries) {
  print(country)
  p <- compare_fits(list(fit, fit1, fit2, fit3), country, c("No shocks", "Shocks tau=1e-2", "Shocks tau=1e-4", "Shocks tau=1e-6"))
  print(p)
}
dev.off()

fit_cv <- read_rds("./fit_shock_cv_2010_spline_0")
fit_cv_1e2 <- read_rds("./fit_shock_cv_2010_shock_0.01")
fit_cv_1e4 <- read_rds("./fit_shock_cv_2010_shock_1e.04")
fit_cv_1e6 <- read_rds("./fit_shock_cv_2010_shock_1e.06")
  
countries <- sort(unique(fit_cv$data$name_country))
pdf("plots/mcpr_cv_2010_comparison.pdf", width = 7, height = 4)
for(country in countries) {
  print(country)
  p <- compare_fits(list(fit_cv, fit_cv_1e2, fit_cv_1e4, fit_cv_1e6), country, c("No shocks", "Shocks tau=1e-2", "Shocks tau=1e-4", "Shocks tau=1e-6")) +
    geom_vline(xintercept = 2010, lty = 2, alpha = 0.5)
  print(p)
}
dev.off()
