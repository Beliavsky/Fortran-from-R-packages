! SPDX-License-Identifier: MIT
program ewens_demo
  use ewens, only : dp, i8, ewens_seed, rewens, dewens, number_of_classes, ewens_mle, ewens_k_exact
  implicit none
  integer, allocatable :: x(:)

  call ewens_seed(2923_i8)
  x = rewens(24, 1.0_dp)

  print '(a,i0)', 'sample size: ', size(x)
  print '(a,i0)', 'number of classes: ', number_of_classes(x)
  print '(a,es14.6)', 'Ewens PMF of frequency profile: ', dewens(x, 1.0_dp)
  print '(a,f12.6)', 'MLE theta: ', ewens_mle(x)
  print '(a,f12.6)', 'E[K] at theta=1: ', ewens_k_exact(size(x), 1.0_dp)
end program ewens_demo
