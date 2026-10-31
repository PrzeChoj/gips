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
