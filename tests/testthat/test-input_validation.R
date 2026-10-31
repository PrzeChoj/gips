test_that("gips rejects nonfinite S with a clear error", {
  for (multi in c(FALSE, TRUE)) {
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      S <- diag(c(1, value))
      if (multi) S <- list(diag(2), S)

      expect_error(
        gips(S, if (multi) c(10, 12) else 10),
        "`S` matrix must contain only finite values.",
        fixed = TRUE, class = "rlang_error"
      )
    }
  }
})

test_that("validate_gips rejects nonfinite S after an attribute change", {
  for (multi in c(FALSE, TRUE)) {
    g <- if (multi) {
      gips(list(diag(2), diag(2)), c(10, 12))
    } else {
      gips(diag(2), 10)
    }

    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      S <- diag(c(1, value))
      attr(g, "S") <- if (multi) list(diag(2), S) else S

      expect_error(
        validate_gips(g),
        "`S` matrix must contain only finite values.",
        fixed = TRUE, class = "rlang_error"
      )
    }
  }
})

test_that("gips rejects NA in D_matrix with a clear error", {
  for (multi in c(FALSE, TRUE)) {
    S <- if (multi) list(diag(2), diag(2)) else diag(2)
    D <- diag(c(1, NA_real_))
    if (multi) D <- list(diag(2), D)

    expect_error(
      gips(S, if (multi) c(10, 12) else 10, D_matrix = D),
      "`D_matrix` must not contain any `NA`s.",
      fixed = TRUE, class = "rlang_error"
    )
  }
})

test_that("validate_gips rejects NA in D_matrix after an attribute change", {
  for (multi in c(FALSE, TRUE)) {
    g <- if (multi) {
      gips(list(diag(2), diag(2)), c(10, 12))
    } else {
      gips(diag(2), 10)
    }
    D <- diag(c(1, NA_real_))
    attr(g, "D_matrix") <- if (multi) list(diag(2), D) else D

    expect_error(
      validate_gips(g),
      "`D_matrix` must not contain any `NA`s.",
      fixed = TRUE, class = "rlang_error"
    )
  }
})

test_that("gips and validate_gips reject infinite delta", {
  for (multi in c(FALSE, TRUE)) {
    S <- if (multi) list(diag(2), diag(2)) else diag(2)
    n <- if (multi) c(10, 12) else 10
    g <- gips(S, n)
    
    for (value in c(Inf, -Inf)) {
      delta <- if (multi) c(3, value) else value
      
      expect_error(
        gips(S, n, delta = delta),
        "`delta` must be finite.",
        fixed = TRUE, class = "rlang_error"
      )
      
      attr(g, "delta") <- delta
      expect_error(
        validate_gips(g),
        "`delta` must be finite.",
        fixed = TRUE, class = "rlang_error"
      )
    }
  }
})

test_that("gips and validate_gips reject nonnumeric D_matrix", {
  invalid_D <- list(
    matrix("x", 2, 2),
    matrix(list(2, 1, 1, 2), 2)
  )
  
  for (multi in c(FALSE, TRUE)) {
    S <- if (multi) list(diag(2), diag(2)) else diag(2)
    n <- if (multi) c(10, 12) else 10
    g <- gips(S, n, D_matrix = NULL)
    
    for (D in invalid_D) {
      if (multi) D <- list(diag(2), D)
      
      expect_error(
        gips(S, n, D_matrix = D),
        "`D_matrix` must be a numeric matrix.",
        fixed = TRUE, class = "rlang_error"
      )
      
      attr(g, "D_matrix") <- D
      expect_error(
        validate_gips(g),
        "`D_matrix` must be a numeric matrix.",
        fixed = TRUE, class = "rlang_error"
      )
    }
  }
})
