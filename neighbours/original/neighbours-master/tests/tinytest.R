if (requireNamespace("tinytest", quietly = TRUE))
    tinytest.results <-
        tinytest::test_package("neighbours",
                               color = interactive(),
                               verbose = 1)
