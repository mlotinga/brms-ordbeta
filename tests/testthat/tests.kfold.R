context("Tests for kfold and kfold helper functions")

test_that("kfold.brmsfit keeps fixed internal k_threshold", {
  skip_on_cran()

  set.seed(1234)
  dat <- data.frame(
    x = rnorm(8),
    y = rnorm(8)
  )
  fit <- suppressWarnings(brm(
    y ~ x, data = dat, family = gaussian(),
    chains = 1, iter = 400, warmup = 200, refresh = 0
  ))

  kfold1 <- kfold(
    fit, K = 2, chains = 1, iter = 300, warmup = 150
  )
  kfold3 <- kfold(
    fit, K = 2, chains = 1, iter = 300, warmup = 150,
    r_eff = NULL
  )

  # 'dims' and 'k_threshold' are in attributes of kfold output
  expect_all_true(c("dims", "k_threshold") %in% names(attributes(kfold1)))

  # r_eff is 1 by default
  expect_all_true(kfold1$diagnostics$r_eff == 1)
  # if r_eff is NULL use loo::relative_eff()
  expect_all_false(kfold3$diagnostics$r_eff == 1)

  # k_threshold follows loo sample-size-dependent rule
  expected_threshold <- min(1 - 1 / log10(dim(kfold1)[1]), 0.7)
  expect_equal(attr(kfold1, "k_threshold"), expected_threshold)

  # expected length of pareto-k, n_eff, r_eff equal to pointwise
  expected_len <- nrow(kfold1$pointwise)
  expect_equal(length(kfold1$diagnostics$pareto_k), expected_len)
  expect_equal(length(kfold1$diagnostics$n_eff), expected_len)
  expect_equal(length(kfold1$diagnostics$r_eff), expected_len)

  # diagnostics follow loo-style names
  expect_equal(
    sort(names(kfold1$diagnostics)),
    sort(c("pareto_k", "n_eff", "r_eff"))
  )
})

test_that(".kfold_r_eff repeats scalar inputs", {
  log_weights <- matrix(log(c(1, 2, 3, 4, 5, 6, 7, 8)), nrow = 2, ncol = 4)
  r_eff <- 0.75

  out <- brms:::.kfold_r_eff(log_weights, r_eff = r_eff, chains = 2)
  expect_equal(out, rep(r_eff, NCOL(log_weights)))
})

test_that(".kfold_r_eff errors on invalid vector length", {
  log_weights <- matrix(log(c(1, 2, 3, 4, 5, 6, 7, 8)), nrow = 2, ncol = 4)

  expect_error(
    brms:::.kfold_r_eff(log_weights, r_eff = rep(0.5, 3), chains = 2),
    "'r_eff' must have one value or one value per observation.",
    fixed = TRUE
  )
})

test_that(".kfold_r_eff replaces NA values", {
  log_weights <- matrix(log(c(1, 2, 3, 4, 5, 6, 7, 8)), nrow = 2, ncol = 4)
  r_eff <- c(0.25, NA_real_, 0.25, 0.25)

  expect_message(
    out <- brms:::.kfold_r_eff(log_weights, r_eff = r_eff, chains = 2),
    "Replacing NAs in `r_eff` with 1s",
    fixed = TRUE
  )
  expect_equal(out[2], 1)
})

test_that(".kfold_r_eff matches loo::relative_eff when r_eff is NULL", {
  log_weights <- matrix(
    log(c(0.8, 1.2, 0.7, 1.1, 1.5, 0.9, 1.3, 0.6)),
    nrow = 4, ncol = 2
  )
  chains <- 2
  chain_id <- rep(1:chains, each = NROW(log_weights) / chains)

  out <- brms:::.kfold_r_eff(log_weights, r_eff = NULL, chains = chains)
  ref <- loo::relative_eff(exp(log_weights), chain_id = chain_id)

  expect_equal(out, ref, tolerance = 1e-12)
})

test_that(".kfold_n_eff matches loo ESS definition", {
  log_weights <- matrix(
    log(c(0.2, 1.1, 0.7, 1.5, 0.4, 0.8, 1.9, 0.6)),
    nrow = 4, ncol = 2
  )
  r_eff <- c(0.5, 0.8)

  out <- brms:::.kfold_n_eff(log_weights = log_weights, r_eff = r_eff)

  norm_const_log <- matrixStats::colLogSumExps(log_weights)
  log_weights_norm <- sweep(log_weights, 2, norm_const_log, check.margin = FALSE)
  weights_norm <- exp(log_weights_norm)
  ref <- 1 / colSums(weights_norm^2) * r_eff

  expect_equal(out, ref, tolerance = 1e-12)
})

test_that(".kfold_k_threshold matches loo threshold formula", {
  for (S in c(100, 320, 1000, 2200, 10000)) {
    expect_equal(
      brms:::.kfold_k_threshold(S),
      min(1 - 1 / log10(S), 0.7),
      tolerance = 1e-12
    )
  }
})
