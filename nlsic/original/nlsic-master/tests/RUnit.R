Sys.unsetenv("R_TESTS")
if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("nlsic", quietly=TRUE)) {
  library(RUnit)
  #library(nlsic)
  testSuite <- RUnit::defineTestSuite(
    name = "nlsic unit tests",
    dirs = system.file("unitTests", package = "nlsic"),
    testFuncRegexp = "^[Tt]est.+",
    rngKind = "default",
    rngNormalKind = "default"
  )
  tests <- RUnit::runTestSuite(testSuite)

  RUnit::printTextProtocol(tests)

  if (RUnit::getErrors(tests)$nFail > 0) stop("RUnit test failure")
  if (RUnit::getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
