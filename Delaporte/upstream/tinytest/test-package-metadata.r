# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# Only test at home. rhub valgrind complains and it should not affect covr.

if (Sys.info()["nodename"] == "HOMEDESKTOP") {

  # Setup
  myPkgs <- c("Delaporte",
              "lamW",
              "Pade",
              "revss",
              "minimaxApprox",
              "MBBEFDLite")
  thisPkg <- "Delaporte"
  otrPkgs <- setdiff(myPkgs, thisPkg)

  # Universal code
  pV <- packageVersion(thisPkg)
  pD <- packageDate(thisPkg)
  cit <- toBibtex(citation(thisPkg))
  lct <- length(cit)
  nws <- news(package = thisPkg)
  lnw <- length(nws)

  # Test CITATION has most recent package version
  expect_true(any(grepl(pV, cit, fixed = TRUE)))

  # Test NEWS has most recent package version
  expect_true(any(grepl(pV, nws, fixed = TRUE)))

  # Test that NEWS has an entry with DESCRIPTION's Date
  expect_true(any(grepl(pD, nws, fixed = TRUE)))

  # Test that CITATION doesn't contain the name of any other of my packages
  expect_false(any(vapply(otrPkgs, grepl, logical(lct), x = cit, fixed = TRUE)))

  # Test that NEWS doesn't contain the name of any other of my packages
  expect_false(any(vapply(otrPkgs, grepl, logical(lnw), x = nws, fixed = TRUE)))
}
