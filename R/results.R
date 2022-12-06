library(tidyverse)
devtools::load_all("../fpemplus")

fit0.01 <- read_rds("./fit_shock0.01")
fit0.01$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "model_3c539d0c08595f54f657c753cfa6e15b-202211101859")))

fit0.001 <- read_rds("./fit_shock0.001")
fit0.001$samples <- cmdstanr::as_cmdstan_fit(files = str_c("server-output/", list.files(path = "server-output/", pattern = "model_d9f3aceb315b96509a213ead9b144cb9-202211152048")))

shocks_raw <- spread_draws(fit0.01$samples, shock[c,t]) 

mean(shocks_raw$shock < -0.401) * 100

shocks <- shocks_raw %>%
  median_qi(shock) %>%
  left_join(fit0.01$country_index) %>%
  left_join(fit0.01$time_index)

largest_shocks <- shocks %>%
  arrange(.lower) %>% filter(.lower < -0.8) %>% 
  distinct(name_country)

largest_shocks

plot_indicator(fit0.01, largest_shocks$name_country)

plot_indicator(fit0.01, "Rwanda")
ggsave("plots/rwanda_0.01.pdf", width = 8, height = 5)

plot_indicator(fit0.01, "Timor-Leste")
ggsave("plots/timor_leste_0.01.pdf", width = 8, height = 5)

plot_indicator(fit0.001, "Rwanda")

plot_indicator(fit0.01, "Somalia")
plot_indicator(fit0.001, "Somalia")


shocks0.01 <- spread_draws(fit0.01$samples, shock[c,t]) %>%
  median_qi(shock) %>%
  left_join(fit0.01$country_index) %>%
  left_join(fit0.01$time_index)

shocks %>%
  filter(name_country == "Rwanda") %>%
  ggplot(aes(x = year, y = shock)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
  scale_fill_brewer()
