library(tidyverse)
devtools::load_all("../fpemplus")

fit0.01 <- read_rds("./fit_shock0.01")
fit0.01$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "model_3c539d0c08595f54f657c753cfa6e15b-202211101859")))

fit0.001 <- read_rds("./fit_shock0.001")
fit0.001$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "model_d9f3aceb315b96509a213ead9b144cb9-202211152048")))

shocks <- spread_draws(fit0.001$samples, shock[c,t]) %>%
  median_qi(shock) %>%
  left_join(fit0.001$country_index) %>%
  left_join(fit0.001$time_index)

largest_shocks <- shocks %>%
  arrange(.lower) %>% filter(.lower < -0.8) %>% 
  distinct(name_country)

plot_indicator(fit0.001, largest_shocks$name_country)
