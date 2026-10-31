test_that("gips and validate_gips require at least two variables", {
  for (multi in c(FALSE, TRUE)) {
    S <- if (multi) list(diag(2), diag(2)) else diag(2)
    n <- if (multi) c(10, 12) else 10
    expect_silent(g <- gips(S, n))
    expect_silent(validate_gips(g))

    for (p in 0:1) {
      S <- diag(p)
      if (multi) S <- list(S, S)

      expect_error(
        gips(S, n),
        "at least two variables (p >= 2)",
        fixed = TRUE, class = "rlang_error"
      )

      attr(g, "S") <- S
      expect_error(
        validate_gips(g),
        "at least two variables (p >= 2)",
        fixed = TRUE, class = "rlang_error"
      )
    }
  }
})

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

test_that("multi-sample D_matrix defaults and complete lists remain valid", {
  S <- list(diag(c(1, 3)), diag(c(4, 8)))
  n <- c(10, 12)
  expect_silent(g <- gips(S, n, D_matrix = NULL))
  expect_equal(attr(g, "D_matrix"), list(2 * diag(2), 6 * diag(2)))
  expect_silent(validate_gips(g))

  D <- list(diag(2), 3 * diag(2))
  expect_silent(g <- gips(S, n, D_matrix = D))
  expect_equal(attr(g, "D_matrix"), D)
  expect_silent(validate_gips(g))
})

test_that("gips and validate_gips reject NULL elements in a D_matrix list", {
  S <- list(diag(2), 4 * diag(2))
  n <- c(10, 12)
  g <- gips(S, n)

  for (i in seq_along(S)) {
    D <- list(diag(2), 4 * diag(2))
    D[i] <- list(NULL)
    expected_message <- paste0("You provided `D_matrix[[", i, "]]` as `NULL`.")
    expect_error(
      gips(S, n, D_matrix = D), expected_message,
      fixed = TRUE, class = "rlang_error"
    )

    attr(g, "D_matrix") <- D
    expect_error(
      validate_gips(g), expected_message,
      fixed = TRUE, class = "rlang_error"
    )
  }
})

test_that("multi-sample D_matrix errors identify the invalid matrix", {
  S <- list(diag(2), diag(2))

  for (i in seq_along(S)) {
    D <- list(diag(2), diag(2))
    D[[i]] <- diag(c(1, 0))
    expect_error(
      gips(S, c(10, 12), D_matrix = D),
      paste0(
        "In `D_matrix[[", i, "]]`: ",
        "`D_matrix` must either be `NULL` or a positive-definite matrix."
      ),
      fixed = TRUE, class = "rlang_error"
    )
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

test_that("gips and validate_gips require exactly symmetric S", {
  for (multi in c(FALSE, TRUE)) {
    S <- diag(2)
    input <- if (multi) list(diag(2), S) else S
    n <- if (multi) c(10, 12) else 10
    expect_silent(g <- gips(input, n))
    expect_silent(validate_gips(g))

    S[1, 2] <- .Machine$double.eps
    input <- if (multi) list(diag(2), S) else S
    expect_error(gips(input, n), "symmetric", fixed = TRUE, class = "rlang_error")
    attr(g, "S") <- input
    expect_error(validate_gips(g), "symmetric", fixed = TRUE, class = "rlang_error")

    S <- (S + t(S))/2
    input <- if (multi) list(diag(2), S) else S
    expect_silent(gips(input, n))
    attr(g, "S") <- input
    expect_silent(validate_gips(g))
  }
})

test_that("the symmetrization hint uses an absolute threshold of 1e-2", {
  differences <- c(0.005, 0.01, 0.02)
  expected_hints <- c(TRUE, TRUE, FALSE)
  for (scale in c(1e-12, 1, 1e12)) {
    for (i in seq_along(differences)) {
      S <- diag(scale, 2)
      S[1, 2] <- differences[i]
      checked <- check_S_matrix(S)
      expect_match(paste(checked$abort_text, collapse = " "), "symmetric", fixed = TRUE)
      has_hint <- any(grepl("S <- (S + t(S))/2", checked$abort_text, fixed = TRUE))
      expect_identical(has_hint, expected_hints[i])
      expect_equal(checked$additional_info, as.integer(expected_hints[i]))
    }
  }
})

test_that("gips and validate_gips require exactly symmetric D_matrix", {
  for (multi in c(FALSE, TRUE)) {
    S <- if (multi) list(diag(2), diag(2)) else diag(2)
    n <- if (multi) c(10, 12) else 10
    expect_silent(gips(S, n, D_matrix = NULL))
    D <- diag(c(2, 3))
    input <- if (multi) list(diag(2), D) else D
    expect_silent(g <- gips(S, n, D_matrix = input))
    expect_silent(validate_gips(g))

    D[1, 2] <- .Machine$double.eps
    input <- if (multi) list(diag(2), D) else D
    err <- tryCatch(gips(S, n, D_matrix = input), error = identity)
    expect_s3_class(err, "rlang_error")
    expect_match(conditionMessage(err), "symmetric", fixed = TRUE)
    expect_match(conditionMessage(err), "1 problems identified", fixed = TRUE)
    attr(g, "D_matrix") <- input
    expect_error(validate_gips(g), "symmetric", fixed = TRUE, class = "rlang_error")

    D <- (D + t(D))/2
    input <- if (multi) list(diag(2), D) else D
    expect_silent(gips(S, n, D_matrix = input))
    attr(g, "D_matrix") <- input
    expect_silent(validate_gips(g))
  }
})

test_that("the D_matrix symmetrization hint uses an absolute threshold of 1e-2", {
  differences <- c(0.005, 0.01, 0.02)
  expected_hints <- c(TRUE, TRUE, FALSE)
  for (scale in c(1e-12, 1, 1e12)) {
    for (i in seq_along(differences)) {
      D <- diag(scale, 2)
      D[1, 2] <- differences[i]
      checked <- check_D_matrix(D, 2)
      expect_length(checked, 2L)
      expect_match(paste(checked, collapse = " "), "symmetric", fixed = TRUE)
      has_hint <- any(grepl("D_matrix <- (D_matrix + t(D_matrix))/2", checked, fixed = TRUE))
      expect_identical(has_hint, expected_hints[i])
    }
  }
})
