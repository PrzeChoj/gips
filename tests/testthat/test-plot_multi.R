test_that("plot() with type = 'heatmap' works on multi-sample gips (non-optimized)", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  
  # Should produce a ggplot with multiple panels (facets) for each group
  p <- plot(g, type = "heatmap")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with type = 'MLE' works on multi-sample gips (non-optimized)", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  
  # Should produce a ggplot with multiple panels (facets) for each group
  p <- plot(g, type = "MLE")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with type = 'block_heatmap' works on multi-sample gips (non-optimized)", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  
  # Should produce a ggplot with multiple panels (facets) for each group
  p <- plot(g, type = "block_heatmap")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with type = 'best' works on multi-sample gips (optimized with BF)", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  # Should produce a convergence plot
  p <- plot(g_map, type = "best")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with type = 'all' works on multi-sample gips (optimized with BF)", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  # Should produce a convergence plot
  p <- plot(g_map, type = "all")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with type = 'both' works on multi-sample gips (optimized with BF)", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  # Should produce a convergence plot with both trajectories
  p <- plot(g_map, type = "both")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with type = 'n0' works on multi-sample gips (optimized with MH)", {
  g <- gips(list(diag(5), 2*diag(5)), c(15L, 18L))
  g_map <- find_MAP(g, optimizer = "MH", max_iter = 2, show_progress_bar = FALSE)
  
  # Should produce an n0 trajectory plot
  p <- plot(g_map, type = "n0")
  expect_s3_class(p, "ggplot")
})


test_that("plot() default behavior works on non-optimized multi-sample gips (type = 'heatmap')", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  
  # Should default to 'heatmap' and produce a message
  expect_message(p <- plot(g), regexp = "automatically set")
  expect_s3_class(p, "ggplot")
})


test_that("plot() default behavior works on optimized multi-sample gips (type = 'both')", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  # Should default to 'both' and produce a message
  expect_message(p <- plot(g_map), regexp = "automatically set")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with custom convergence parameters works on multi-sample optimization", {
  g <- gips(list(diag(2), 2*diag(2)), c(10L, 12L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  # Custom parameters
  p <- plot(g_map, type = "both", logarithmic_x = TRUE, logarithmic_y = FALSE,
            title_text = "Custom Title", xlabel = "Custom X", ylabel = "Custom Y",
            color = c("purple", "orange"), show_legend = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot() with larger multi-sample works (p=3, G=3)", {
  g <- gips(list(diag(2), 2*diag(2), 4*diag(2)), c(20L, 25L, 30L))
  
  # Heatmap should work
  p <- plot(g, type = "heatmap")
  expect_s3_class(p, "ggplot")
})


test_that("plot() with larger optimized multi-sample works (p=3, G=3, BF)", {
  g <- gips(list(diag(2), 2*diag(2), 4*diag(2)), c(20L, 25L, 30L))
  g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  
  # Convergence plots should work
  p1 <- plot(g_map, type = "all")
  p2 <- plot(g_map, type = "best")
  p3 <- plot(g_map, type = "both")
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})


test_that("plot() with >10 groups issues a warning and returns a ggplot in non-interactive mode", {
  # Build 11 identity groups (non-interactive session: skips prompt, just warns)
  S_list <- replicate(11, diag(2), simplify = FALSE)
  n_vec <- rep(10L, 11)

  g <- gips(S_list, n_vec)

  expect_warning(
    p <- plot(g, type = "heatmap"),
    regexp = "cluttered"
  )
  expect_s3_class(p, "ggplot")
  
  expect_warning(
    p <- plot(g, type = "block_heatmap"),
    regexp = "cluttered"
  )
  expect_s3_class(p, "ggplot")
})
