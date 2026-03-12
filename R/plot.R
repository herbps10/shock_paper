library(patchwork)

plot_comparison <- function(fit, country) {
  fit$posteriors$temporal |>
    filter(name == country) |>
    ggplot(aes(x = year, y = `50%`, color = variable)) +
    geom_line() +
    geom_line(aes(y = `2.5%`), lty = 2) +
    geom_line(aes(y = `97.5%`), lty = 2) +
    ggtitle(country)
}


plot_shock <- function(fit, areas = fit$country_index$name) {
  fit$posteriors$temporal |>
    filter(variable == "shock2", name %in% areas) |>
    ggplot(aes(x = year, y = `50%`)) +
    geom_errorbar(aes(ymin = `1%`, ymax = `99%`), width = 0) +
    geom_point() +
    facet_wrap(~name)
}

plot_P_tilde <- function(fit) {
  fit$samples |> spread_draws(P_tilde2[c]) |>
    left_join(fit$country_index) |>
    group_by(name) |>
    mutate(P_tilde2 = 15 + P_tilde2) |>
    ggplot(aes(x = P_tilde2)) +
    geom_density() +
    facet_wrap(~name)
}

plot_mean_transition <- function(fit) {
  fit$posteriors$transition_function_mean |>
    ggplot(aes(x = 15 + x * (85 - 15), y = transition_function_mean)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
    scale_fill_brewer()
}


plot_transition <- function(fit, areas = c()) {
  if(length(areas) == 0) {
    areas <- unique(fit$data$name)
  }
  
  
  data <- fit$data |>
    filter(name %in% areas) |>
    arrange(name, period) |>
    group_by(name) |>
    mutate(diff = c(diff(e0), NA))
  
  if(fit$model == "logistic" || fit$model == "logistic_shock") {
    fit$posteriors$transition_functions |>
      filter(name %in% areas) |>
      ggplot(aes(x = x, y = transition_function_pred)) +
      geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
      scale_fill_brewer() +
      geom_point(data = data, aes(x = e0, y = diff)) +
      facet_wrap(~name)
  }
  else {
    fit$posteriors$transition_functions |>
      filter(name %in% areas) |>
      ggplot(aes(x = 15 + x * (110 - 15), y = transition_function_pred)) +
      geom_lineribbon(aes(ymin = .lower, ymax = .upper)) +
      scale_fill_brewer() +
      geom_point(data = data, aes(x = e0, y = diff)) +
      facet_wrap(~name)
  }
}

plot_with_shocks <- function(fit, area) {
  p1 <- fit$data |>
    filter(name == area) |>
    ggplot(aes(x = year, y = e0)) +
    geom_line() +
    geom_point() +
    xlim(c(1950, 2015)) +
    labs(x = "", y = expression(e[0])) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(label = "", subtitle = area)
  
  p2 <- fit$posteriors$temporal |>
    filter(variable == "shock", name == area) |>
    ggplot(aes(x = year, y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
    geom_point() +
    xlim(c(1950, 2015)) +
    labs(x = "Year", y = expression(delta[ct]))
  
  p1 / p2 + plot_layout(heights = c(5, 2))
}

plot_shock_corrected <- function(fit, areas = fit$country_index$name) {
  threshold <- 2 * fit$samples$summary("epsilon_scale")$median
  fit$data |> 
    filter(name %in% areas) |>
    left_join(
      fit$posteriors$temporal |> 
      filter(variable == "shock", name %in% areas, `97.5%` < -threshold)) |> 
    ggplot(aes(x = year, y = e0)) + 
    geom_point(aes(shape = "Observations", color = "Observations"), size = 0.5) + 
    geom_point(aes(shape = "Shock-corrected", y = e0 - `50%`, color = "Shock-corrected"), size = 0.5) +
    geom_segment(aes(x = year, xend = year, y = e0, yend = e0 - `50%`), lty = 3, alpha = 0.5) +
    geom_errorbar(aes(color = "Shock-corrected", ymin = e0 - `97.5%`, ymax = e0 - `2.5%`, width = 0)) +
    scale_color_manual(values = c("black", "blue")) +
    guides(shape = FALSE) +
    facet_wrap(~name) +
    labs(color = "", x = "Year", y = expression(e[0])) +
    pub_theme
}


plot_temporal <- function(x, fit, areas = c(), plot_data = FALSE, color_sources = FALSE) {
  if(length(areas) == 0) areas <- fit$country_index[[fit$area]]
  
  post <- fit$posteriors$temporal %>%
    filter(!!sym(fit$area) %in% areas) %>%
    filter(.data$variable == x)
  
  
  p <- ggplot2::ggplot(post, aes_string(x = fit$year, y = "`50%`")) +
    ggplot2::geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = "95%")) +
    ggplot2::geom_ribbon(aes(ymin = .data$`10%`,  ymax = .data$`90%`, fill = "80%")) +
    ggplot2::geom_ribbon(aes(ymin = .data$`25%`,  ymax = .data$`75%`, fill = "50%")) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(direction = -1) +
    ggplot2::facet_wrap(vars(!!sym(fit$area))) +
    labs(fill = "Posterior\nQuantile", x = fit$year)
  
  if(plot_data == TRUE) {
    data_tibble <- tibble::tibble(i = 1:length(fit$held_out), held_out = as.logical(fit$held_out))
    data <- fit$data %>%
      dplyr::mutate(i = 1:n()) %>%
      dplyr::left_join(data_tibble, by = "i")
    
    data[[fit$source]] <- factor(data[[fit$source]])
    
    filtered_data <- data %>%
      dplyr::filter(!!sym(fit$area) %in% areas)
    
    if(!is.null(fit$se)) {
      filtered_data <- filtered_data %>% mutate(
        lower = truncnorm::qtruncnorm(0.025, mean = !!sym(fit$y), sd = !!sym(fit$se), a = 0, b = 1),
        upper = truncnorm::qtruncnorm(0.975, mean = !!sym(fit$y), sd = !!sym(fit$se), a = 0, b = 1)
      )
    }
    
    some_held_out <- any(fit$held_out == 1)
    
    if(color_sources == TRUE) {
      if(some_held_out == TRUE) {
        point_aes <- aes_string(y = fit$y, color = fit$source, shape = "held_out")
      }
      else {
        point_aes <- aes_string(y = fit$y, color = fit$source)
      }
      if(!is.null(fit$se)) {
        p <- p + ggplot2::geom_errorbar(aes_string(y = fit$y, ymin = "lower", ymax = "upper", color = fit$source), alpha = 0.3, width = 0, data = filtered_data)
      }
      p <- p + ggplot2::geom_point(point_aes, data = filtered_data, alpha = 0.7)
    }
    else {
      if(some_held_out == TRUE) {
        point_aes <- aes_string(y = fit$y, shape = "held_out")
      }
      else {
        point_aes <- aes_string(y = fit$y)
      }
      if(!is.null(fit$se)) {
        p <- p + ggplot2::geom_errorbar(aes_string(y = fit$y, ymin = "lower", ymax = "upper"), alpha = 0.3, width = 0, data = filtered_data)
      }
      p <- p + ggplot2::geom_point(point_aes, data = filtered_data, alpha = 0.7)
    }
  }
  
  p
}
