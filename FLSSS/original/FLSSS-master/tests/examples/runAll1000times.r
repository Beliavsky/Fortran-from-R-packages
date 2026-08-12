

j = 0L
for(i in 1L : 300L)
{
  cat(i, "")
  source("tests/examples/FLSSS.r")
  j = 1L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/FLSSSmultiset.r")
  j = 2L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/GAP.r")
  j = 3L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/GAPintegerized.r")
  j = 4L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/mFLSSSpar.r")
  j = 5L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/mFLSSSparImposeBounds.r")
  j = 6L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/mFLSSSparImposeBoundsIntegerized.r")
  j = 7L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/mFLSSSparIntegerized.r")
  j = 8L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/mmKnapsack.r")
  j = 9L; save(j, file = "tests/examples/j.Rdata")

  source("tests/examples/mmKnapsackIntegerized.r")
  j = 10L; save(j, file = "tests/examples/j.Rdata")
}


# for(i in 1:10)
# {
#   rm(list = ls())
#   cat(1, "")
# }


