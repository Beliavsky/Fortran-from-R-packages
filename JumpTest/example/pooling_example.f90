! SPDX-License-Identifier: MIT
program pooling_example
  use jumptest, only : dp, adjp_result, ppool
  implicit none

  real(dp), parameter :: pvalues(5, 3) = reshape([ &
    0.01_dp, 0.20_dp, 0.50_dp, 0.05_dp, 0.30_dp, &
    0.02_dp, 0.15_dp, 0.40_dp, 0.10_dp, 0.35_dp, &
    0.03_dp, 0.25_dp, 0.60_dp, 0.08_dp, 0.25_dp], [5, 3])
  type(adjp_result) :: pooled

  call ppool(pvalues, pooled, 'SD')
  print '(a)', 'Dependent Stouffer pooled p-values:'
  print '(5(f10.6,1x))', pooled%pvalue
  print '(a)', 'BH-adjusted values:'
  print '(5(f10.6,1x))', pooled%adjp
end program pooling_example
