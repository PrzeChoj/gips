# Tests for internal plot helpers and plot.gips() convergence output


# ── deduplicate_step_runs() ───────────────────────────────────────────────────

test_that("deduplicate_step_runs() keeps only first of each run plus last row", {
  df <- data.frame(step = 1:6, value = c(1, 1, 2, 2, 2, 3), series = "x")

  result <- gips:::deduplicate_step_runs(df)

  # First of run-1 (step 1), first of run-2 (step 3), first+last of run-3
  # (step 6 is both first-of-run and last, so just one row)
  expect_equal(result$step,  c(1, 3, 6))
  expect_equal(result$value, c(1, 2, 3))
})

test_that("deduplicate_step_runs() always keeps the last row, even when it repeats", {
  df <- data.frame(step = 1:4, value = c(1, 2, 2, 2), series = "x")

  result <- gips:::deduplicate_step_runs(df)

  expect_equal(result$step,  c(1, 2, 4))
  expect_equal(result$value, c(1, 2, 2))
})

test_that("deduplicate_step_runs() keeps all rows when there are no repeats", {
  df <- data.frame(step = 1:5, value = c(10, 20, 30, 40, 50), series = "x")

  result <- gips:::deduplicate_step_runs(df)

  expect_equal(nrow(result), 5L)
  expect_equal(result$step,  1:5)
  expect_equal(result$value, c(10, 20, 30, 40, 50))
})

test_that("deduplicate_step_runs() keeps only first and last row when all values are equal", {
  df <- data.frame(step = 1:5, value = rep(7, 5), series = "x")

  result <- gips:::deduplicate_step_runs(df)

  expect_equal(nrow(result), 2L)
  expect_equal(result$step,  c(1, 5))
  expect_equal(result$value, c(7, 7))
})

test_that("deduplicate_step_runs() returns a single-row data frame unchanged", {
  df <- data.frame(step = 1L, value = 42, series = "x")

  result <- gips:::deduplicate_step_runs(df)

  expect_equal(nrow(result), 1L)
  expect_equal(result$step,  1L)
  expect_equal(result$value, 42)
})

test_that("deduplicate_step_runs() preserves all columns", {
  df <- data.frame(step = 1:3, value = c(5, 5, 6), series = "my_series", extra = letters[1:3])

  result <- gips:::deduplicate_step_runs(df)

  expect_named(result, c("step", "value", "series", "extra"))
  # Row 1 (value=5) and row 3 (value=6, the first of the new run, also the last)
  expect_equal(result$extra, c("a", "c"))
})


# ── Step-function correctness (dedup is visually lossless) ───────────────────

# Helper: given steps and values defining a step function, evaluate it at x.
# geom_step(direction = "hv") maps x to the value of the *last* data point
# with step <= x.
eval_step_fn <- function(steps, values, x) {
  idx <- max(which(steps <= x))
  values[idx]
}

test_that("deduplicate_step_runs() is visually lossless for geom_step()", {
  # Build a realistic cummax-like sequence (many repeated values)
  raw <- c(-10, -10, -10, -8, -8, -8, -8, -5, -5, -3, -3, -3, -3, -3, -3, -1)
  n   <- length(raw)

  df_full <- data.frame(step = seq_len(n), value = raw)
  df_dedup <- gips:::deduplicate_step_runs(df_full)

  # The deduped data must be a strict subset of rows
  expect_lt(nrow(df_dedup), nrow(df_full))

  # For every integer x in range, the step-function value must be the same
  for (x in seq_len(n)) {
    v_full  <- eval_step_fn(df_full$step,  df_full$value,  x)
    v_dedup <- eval_step_fn(df_dedup$step, df_dedup$value, x)
    expect_equal(v_full, v_dedup,
      label = paste0("step function value at x = ", x)
    )
  }
})


# ── plot_gips_convergence(): x-extent and row-count after deduplication ───────

# Shared setup: a small gips object optimized with MH so that the stored
# log_posteriori_values vector has many repeated entries (MH rejects most steps).
local({
  g_base <- gips(diag(1.1, 5), 10)
  g_mh <- suppressMessages(find_MAP(g_base, optimizer = "MH", max_iter = 50,
    show_progress_bar = FALSE))

  num_of_steps <- length(
    attr(g_mh, "optimization_info")[["log_posteriori_values"]]
  )

  test_that("convergence plot data x-range spans full iteration count (type = 'all')", {
    skip_if_not_installed("ggplot2")
    gg <- plot(g_mh, type = "all")
    expect_true(inherits(gg, "ggplot"))
    expect_equal(max(gg$data$step), num_of_steps)
    expect_equal(min(gg$data$step), 1L)
  })

  test_that("convergence plot data x-range spans full iteration count (type = 'best')", {
    skip_if_not_installed("ggplot2")
    gg <- plot(g_mh, type = "best")
    expect_equal(max(gg$data$step), num_of_steps)
  })

  test_that("convergence plot data x-range spans full iteration count (type = 'both')", {
    skip_if_not_installed("ggplot2")
    gg <- plot(g_mh, type = "both")
    expect_equal(max(gg$data$step), num_of_steps)
  })

  test_that("convergence plot data has fewer rows than num_of_steps (deduplication works)", {
    skip_if_not_installed("ggplot2")
    # With MH and 50 steps at p=5, there will be many repeated log-posterior
    # values (rejected proposals). After deduplication the data frame must be
    # smaller (if all 50 values happened to be distinct this would fail, but
    # that is astronomically unlikely for MH).
    gg_all <- plot(g_mh, type = "all")
    expect_lt(nrow(gg_all$data), num_of_steps)
  })
})


# ── Continuing MH runs and plotting n0 trajectory ────────────────────────────

test_that("MH -> continue -> continue -> plot(type = 'n0') works end-to-end", {
  skip_if_not_installed("ggplot2")
  # Start with MH optimization
  g_base <- gips(diag(1.1, 5), 10)
  g_mh <- suppressMessages(find_MAP(g_base, optimizer = "MH", max_iter = 10,
    show_progress_bar = FALSE))
  
  expect_true(inherits(g_mh, "gips"))
  
  # First continue call
  g_continue1 <- suppressMessages(find_MAP(g_mh, optimizer = "continue", max_iter = 5,
    show_progress_bar = FALSE))
  
  expect_true(inherits(g_continue1, "gips"))
  
  # Second continue call
  g_continue2 <- suppressMessages(find_MAP(g_continue1, optimizer = "continue", max_iter = 5,
    show_progress_bar = FALSE))
  
  expect_true(inherits(g_continue2, "gips"))
  
  # Plot n0 trajectory
  gg_n0 <- plot(g_continue2, type = "n0", logarithmic_y = FALSE)
  
  expect_true(inherits(gg_n0, "ggplot"))
})
