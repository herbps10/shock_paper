library(patchwork)

plot_comparison <- function(fit, country) {
  fit$posteriors$temporal %>%
    filter(name == country) %>%
    ggplot(aes(x = year, y = `50%`, color = variable)) +
    geom_line() +
    geom_line(aes(y = `2.5%`), lty = 2) +
    geom_line(aes(y = `97.5%`), lty = 2) +
    ggtitle(country)
}


plot_shock <- function(fit, areas = c()) {
  fit$posteriors$temporal %>%
    filter(variable == "shock", name %in% areas) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
    geom_point() +
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

plot_mean_transition <- function(fit) {
  fit$posteriors$transition_function_mean %>%
    ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_mean)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
    scale_fill_brewer()
}


plot_transition <- function(fit, areas = c()) {
  if(length(areas) == 0) {
    areas <- unique(fit$data$name)
  }
  
  fit$posteriors$transition_functions %>%
    filter(name %in% areas) %>%
    ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_pred)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
    scale_fill_brewer() +
    facet_wrap(~name)
}

plot_with_shocks <- function(fit, area) {
  p1 <- fit$data %>%
    filter(name == area) %>%
    ggplot(aes(x = year, y = e0)) +
    geom_line() +
    geom_point() +
    xlim(c(1950, 2015)) +
    labs(x = "", y = expression(e[0])) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(label = "", subtitle = area)
  
  p2 <- fit$posteriors$temporal %>%
    filter(variable == "shock", name == area) %>%
    ggplot(aes(x = year, y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
    geom_point() +
    xlim(c(1950, 2015)) +
    labs(x = "Year", y = expression(delta[ct]))
  
  p1 / p2 + plot_layout(heights = c(5, 2))
}

plot_shock_corrected <- function(fit, areas) {
  threshold <- 2 * fit$samples$summary("epsilon_scale")$median
  fit$data %>% 
    filter(name %in% areas) %>%
    left_join(
      fit$posteriors$temporal %>% 
      filter(variable == "shock", name %in% areas, `97.5%` < -threshold)) %>% 
    ggplot(aes(x = year, y = e0)) + 
    geom_point(aes(shape = "Observations", color = "Observations")) + 
    geom_point(aes(shape = "Shock-corrected", y = e0 - `50%`, color = "Shock-corrected")) +
    geom_segment(aes(x = year, xend = year, y = e0, yend = e0 - `50%`), lty = 3, alpha = 0.5) +
    geom_errorbar(aes(color = "Shock-corrected", ymin = e0 - `97.5%`, ymax = e0 - `2.5%`, width = 0)) +
    scale_color_manual(values = c("black", "blue")) +
    guides(shape = FALSE) +
    facet_wrap(~name) +
    labs(color = "", x = "Year", y = expression(e[0]))
}

