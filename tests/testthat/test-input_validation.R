test_that("gips silently accepts missing or matching axis names", {
  variable_names <- c("height", "weight")
  for (axis_names in list(
    list(NULL, NULL), list(variable_names, NULL),
    list(NULL, variable_names), list(variable_names, variable_names)
  )) {
    S <- diag(c(2, 1))
    dimnames(S) <- axis_names
    for (multi in c(FALSE, TRUE)) {
      input <- if (multi) list(diag(2), S) else S
      n <- if (multi) c(10, 12) else 10
      expect_silent(g <- gips(input, n))
      expect_identical(attr(g, "S"), input)
    }
  }
})

test_that("gips warns about conflicting axis names before optimization without repeating in methods", {
  S <- diag(c(2, 3, 4))
  rownames(S) <- c("height", "weight", "age")
  colnames(S) <- c("weight", "height", "age")

  for (multi in c(FALSE, TRUE)) {
    input <- if (multi) list(diag(3), S) else S
    n <- if (multi) c(10, 12) else 10
    matrix_name <- if (multi) "S[[2]]" else "S"
    expected_warning <- paste0("Row and column names differ in `", matrix_name, "`.")
    expect_warning(
      g <- gips(input, n), expected_warning, fixed = TRUE,
      class = "gips_mismatched_dimnames"
    )
    expect_identical(attr(g, "S"), input)
    expect_silent(validate_gips(g))
    expect_silent(log_posteriori_of_gips(g))
    expect_silent(capture.output(print(g)))
    expect_silent(invisible(summary(g)))

    expect_silent(g_map <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE))
    expect_silent(capture.output(print(g_map)))
    expect_silent(invisible(summary(g_map)))
    expect_identical(attr(g_map, "S"), input)

    # Constructing another object must still report its conflicting names.
    expect_warning(
      gips(input, n), expected_warning, fixed = TRUE,
      class = "gips_mismatched_dimnames"
    )
  }
})

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

test_that("list matrices with scalar numbers get a conversion hint", {
  A <- matrix(list(2, 1, 1, 2), nrow = 2)

  for (multi in c(FALSE, TRUE)) {
    S <- if (multi) list(diag(2), A) else A
    n <- if (multi) c(10, 12) else 10
    conversion <- if (multi) {
      "S[[2]] <- matrix(unlist(S[[2]]), nrow = nrow(S[[2]]), dimnames = dimnames(S[[2]]))"
    } else {
      "S <- matrix(unlist(S), nrow = nrow(S), dimnames = dimnames(S))"
    }
    g <- gips(if (multi) list(diag(2), diag(2)) else diag(2), n)
    attr(g, "S") <- S
    errors <- list(
      tryCatch(gips(S, n), error = identity),
      tryCatch(validate_gips(g), error = identity),
      tryCatch(find_MAP(g, optimizer = "MH", max_iter = 2, show_progress_bar = FALSE),
               error = identity),
      tryCatch(compare_log_posteriories_of_perms(
        "(1,2)", "()", S = S, number_of_observations = n,
        print_output = FALSE
      ), error = identity)
    )

    for (err in errors) {
      expect_s3_class(err, "rlang_error")
      message <- conditionMessage(err)
      expect_match(message, "storage type `list`", fixed = TRUE)
      expect_match(message, conversion, fixed = TRUE)
      expect_match(message, "1 problems identified", fixed = TRUE)
      expect_false(grepl("non-numeric values", message, fixed = TRUE))
    }
  }
})

test_that("optimizer argument validation recognizes a list matrix", {
  expect_error(
    check_find_MAP_arguments(
      S = matrix(list(2, 1, 1, 2), nrow = 2), number_of_observations = 10,
      max_iter = 2, start_perm = gips_perm("()", 2), delta = 3,
      D_matrix = diag(2), was_mean_estimated = FALSE,
      return_probabilities = FALSE, save_all_perms = FALSE, show_progress_bar = FALSE
    ),
    "unlist(S)", fixed = TRUE, class = "rlang_error"
  )
})

test_that("list matrix conversion hints require finite real scalar numbers", {
  invalid_entries <- list(c(1, 2), NULL, "1", list(1), TRUE, NA_real_, NaN, Inf, -Inf, 1i)
  for (entry in invalid_entries) {
    A <- matrix(list(2, 1, 1, entry), nrow = 2)
    err <- tryCatch(gips(A, 14), error = identity)
    expect_s3_class(err, "rlang_error")
    expect_match(conditionMessage(err), "storage type `list`", fixed = TRUE)
    expect_false(grepl("unlist", conditionMessage(err), fixed = TRUE))
  }

  empty <- matrix(list(), nrow = 0, ncol = 0)
  err <- tryCatch(gips(empty, 14), error = identity)
  expect_s3_class(err, "rlang_error")
  expect_false(grepl("unlist", conditionMessage(err), fixed = TRUE))
})

test_that("multi-sample S errors identify the invalid matrix", {
  n <- c(10, 12)
  g <- gips(list(diag(2), diag(2)), n)

  for (i in 1:2) {
    S <- list(diag(2), diag(2))
    S[[i]][1, 1] <- NA_real_
    attr(g, "S") <- S
    errors <- list(
      tryCatch(gips(S, n), error = identity),
      tryCatch(validate_gips(g), error = identity),
      tryCatch(compare_log_posteriories_of_perms(
        "(1,2)", "()", S = S, number_of_observations = n,
        print_output = FALSE
      ), error = identity)
    )

    for (err in errors) {
      expect_s3_class(err, "rlang_error")
      expect_match(
        conditionMessage(err),
        paste0("In `S[[", i, "]]`: `S` matrix must contain only finite values."),
        fixed = TRUE
      )
    }
  }
})

test_that("multi-sample S symmetrization hints refer to the affected matrix", {
  S <- list(diag(2), diag(2))
  n <- c(10, 12)
  g <- gips(S, n)
  S[[2]][1, 2] <- 0.005
  attr(g, "S") <- S
  errors <- list(
    tryCatch(gips(S, n), error = identity),
    tryCatch(validate_gips(g), error = identity),
    tryCatch(compare_log_posteriories_of_perms(
      "(1,2)", "()", S = S, number_of_observations = n,
      print_output = FALSE
    ), error = identity)
  )

  for (err in errors) {
    expect_s3_class(err, "rlang_error")
    expect_match(conditionMessage(err), "In `S[[2]]`:", fixed = TRUE)
    expect_match(conditionMessage(err), "S[[2]] <- (S[[2]] + t(S[[2]]))/2", fixed = TRUE)
    expect_match(conditionMessage(err), "1 problems identified", fixed = TRUE)
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
