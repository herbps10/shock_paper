library(tidyverse)
library(wpp2024)
library(tidybayes)

cmdstanr::set_cmdstan_path("/gpfs/data/diazi07lab/cmdstan-2.38.0/")

source("R/lifeplus.R")
source("R/process_lifeplus.R")

data(UNlocations, package = "wpp2024")
data(e0M1, package = "wpp2024")
data(pop1, package = "wpp2024")
data(include_2010, package = "bayesLife")

included_codes <- include_2010 |> filter(include_code %in% 1:2) |> pull(country_code)

large_countries <- pop1 |>
  filter(`2023` >= 1e3) |>
  pull(name)

datM <- e0M1 |>
  filter(name %in% large_countries) |>
  pivot_longer(cols = `1950`:`2023`, names_to = "period", values_to = "e0") |>
  filter(country_code %in% included_codes) |>
  mutate(year = parse_integer(str_sub(period, 1, 4)),
         source = "WPP2024")

index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

fits <- expand_grid(
  scale_global = c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8),
  model = c("logistic_shock", "logistic"),
  outlier_threshold = c(5, 1e3),
  centered = c(TRUE, FALSE),
  hierarchical = c(TRUE),
  config = c("low")
) |>
  filter(!(model == "logistic_shock" & outlier_threshold == 5)) |>
  filter(!(model == "logistic" & scale_global != 1e-2)) |>
  filter(!(hierarchical == FALSE & centered == TRUE))

fits <- fits[(1:nrow(fits) %% 9 + 1) == index, ]

#random_countries <- c(sample(unique(datM$name), 25), "Niger")
#print(random_countries)
#datM <- datM |> filter(name %in% random_countries)

print(glue::glue("Starting index {index} with {nrow(datM)} rows for {length(unique(datM$name))} countries"))
print(glue::glue("Scale: {fits$scale_global} model: {fits$model} outlier: {fits$outlier_threshold}"))

fits <- fits |>
  mutate(fit = pmap(list(scale_global, model, outlier_threshold, centered, config), function(scale_global, model, outlier_threshold, centered, config) {
    if(config == "high") {
      adapt_delta <- 0.999
      max_treedepth <- 15
      iter_warmup <- 500
      iter_sampling <- 500
    }
    else {
      adapt_delta <- 0.999
      max_treedepth <- 14
      iter_warmup <- 250
      iter_sampling <- 250
    }
    lifeplus(
      datM,
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,

      outlier_threshold = outlier_threshold,
      
      hierarchical = hierarchical,
      centered = centered,
      
      model = model,
      
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,

      parallel_chains = 8,
      chains = 8,

      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      
      output_dir = "/gpfs/scratch/susmah01/shock_paper/draws/",
      
      extra_stan_data = list(
        scale_global = scale_global,
        slab_scale = 10,
        slab_df = 6
      )
    )
  }))

print(glue::glue("Finished index {index} with {nrow(datM)} rows for {length(unique(datM$name))} countries"))

write_rds(fits, glue::glue("/gpfs/scratch/susmah01/shock_paper/fits-{index}.rds"))
