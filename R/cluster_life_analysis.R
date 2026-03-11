library(tidyverse)
library(wpp2024)
#library(BayesTransitionModels)
library(tidybayes)
#library(bayesLife)

cmdstanr::set_cmdstan_path("/gpfs/data/diazi07lab/cmdstan-2.38.0/")

source("R/lifeplus.R")
source("R/process_lifeplus.R")

data(UNlocations, package = "wpp2024")
data(e0M1, package = "wpp2024")
data(pop1, package = "wpp2024")
data(include_2010, package = "bayesLife")

included_codes <- include_2010 |> filter(include_code %in% 1:2) |> pull(country_code)

large_countries <- pop1 |>
  filter(`2020` >= 1e3) |>
  pull(name)

datM <- e0M1 |>
  filter(name %in% large_countries) |>
  pivot_longer(cols = `1950`:`2023`, names_to = "period", values_to = "e0") |>
  filter(country_code %in% included_codes) |>
  mutate(year = parse_integer(str_sub(period, 1, 4)),
         source = "WPP2024")

index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

fits <- expand_grid(
  scale_global = 1e-2,
  model = "logistic_shock",
  outlier_threshold = 1e3
) |>
  bind_rows(
    tibble(scale_global = 1e-2, model = "logistic", outlier_threshold = 5),
    tibble(scale_global = 1e-2, model = "logistic", outlier_threshold = 1e3)
  )

fits <- fits[index, ]

random_countries <- sample(unique(datM$name), 50)
print(random_countries)
datM <- datM |> filter(name %in% random_countries)

print(glue::glue("Starting index {index} with {nrow(datM)} rows for {length(unique(datM$name))} countries"))


fits <- fits |>
  mutate(fit = pmap(list(scale_global, model, outlier_threshold), function(scale_global, model, outlier_threshold) {
    lifeplus(
      datM,
      y = "e0", 
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,

      outlier_threshold = outlier_threshold,
      
      model = model,
      
      adapt_delta = 0.95,
      max_treedepth = 12,

      parallel_chains = 8,
      chains = 8,

      iter_warmup = 100,
      iter_sampling = 200,
      
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
