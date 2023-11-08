
# la using code from life_analysis.R
# wpp2022 from github
# library(devtools)
# options(timeout = 600)
# install_github("PPgp/wpp2022")
#
# # just adding a new install of bayestransition too
# install_github("AlkemaLab/BayesTransitionModels")

# I got errors becasue of int[N] declarations, changed to array[N]
# Error in '/var/folders/pw/nmbdn6h568bdb6hlj9r6rzhcvkpjdp/T/RtmpGiKKLh/model-1248952071ee0.stan', line 17, column 2: Declaration
# of arrays by placing brackets after a variable name was removed in Stan
# 2.33.0. Instead use the array keyword before the type. This can be
# changed automatically using the auto-format flag to stanc

# # automating updates unsuccessful so far
# # because the include files are not in same folder as model file
# # leaving this aside for now
# # see
# # https://mc-stan.org/docs/stan-users-guide/stanc-pretty-printing.html
#
# mod <- cmdstanr::cmdstan_model("stan/life_spline.stan", compile = FALSE)
# mod$format(canonicalize = TRUE)


# just running one setting
fits <- expand_grid(
  scale_global = 1e-2, #c(1e-3, 1e-2, 1e-1),
  # actually paper has 0.1 but posteriors are similar
  num_knots = c(7)
) %>%
  mutate(fit = map(scale_global, function(scale_global) {
    lifeplus(
      datM,
      y = "e0",
      year = "year",
      area = "name",
      source = "source",
      start_year = 1950,
      end_year = 2100,

      model = "shock2",

      spline_degree = 2,
      num_knots = 7,
      hierarchical_splines = c("intercept", "name"),

      adapt_delta = 0.99,
      max_treedepth = 12,
      parallel_chains = 4,
      iter_warmup = 250,
      iter_sampling = 1e3,

      extra_stan_data = list(
        scale_global = scale_global,
        slab_scale = 10,
        slab_df = 6
      )
    )
  }))

saveRDS(fits, "outputs/fits_wshock_la20231019.rds")
#fits <- readRDS("outputs/fits_wshock_la20231019.rds")

fit <- lifeplus(
  datM,
  y = "e0",
  year = "year",
  area = "name",
  source = "source",
  start_year = 1950,
  end_year = 2100,

  model = "spline",

  spline_degree = 2,
  num_knots = 7,
  hierarchical_splines = c("intercept", "name"),

  parallel_chains = 4,
  iter_warmup = 250,
  iter_sampling = 500,

  extra_stan_data = list(
    scale_global = 1e-2,
    slab_scale = 10,
    slab_df = 6
  )
)

saveRDS(fit, "outputs/fit_la20231019.rds")
#fit <- readRDS("outputs/fit_la20231019.rds")

# check prior and post for the hyperpar related to splines coefficients
# in rmd



