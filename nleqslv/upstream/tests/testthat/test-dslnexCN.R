test_that("numeric and character delta specifications give identical output", {
  # Dennis & Schnabel example 6.5.1
  dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1] - 1) + x[2]^3 - 2
    y
  }

  xstart <- c(2, 0.5)

  print_result <- function(z) {
    c(
      paste("x:", paste(format(z$x), collapse = " ")),
      paste("fvec:", paste(format(z$fvec), collapse = " ")),
      paste("message:", z$message),
      paste("ok:", all(abs(z$fvec) <= 1e-8))
    )
  }

  # Collect output for numeric deltas
  temp_numeric <- capture.output({
    for (g in c("dbldog", "pwldog")) {
      for (delta in c(-1.0, -2.0)) {
        z <- nleqslv(
          xstart, dslnex,
          global = g,
          control = list(btol = 0.01, delta = delta, trace = 1)
        )
        print_result(z)
      }
    }
  })

  # Collect output for character deltas
  temp_character <- capture.output({
    for (g in c("dbldog", "pwldog")) {
      for (delta in c("cauchy", "newton")) {
        z <- nleqslv(
          xstart, dslnex,
          global = g,
          control = list(btol = 0.01, delta = delta, trace = 1)
        )
        print_result(z)
      }
    }
  })

  # Compare
  expect_equal(temp_numeric, temp_character)
})

