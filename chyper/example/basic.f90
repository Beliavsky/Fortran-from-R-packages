! SPDX-License-Identifier: MIT
program basic
  use chyper, only : dp, dchyper, pchyper, qchyper, mle_s
  implicit none
  integer, parameter :: n(3) = [12, 13, 14]
  integer, parameter :: m(3) = [7, 8, 9]
  integer, parameter :: observed(7) = [0, 0, 1, 1, 0, 2, 0]

  print '(a,f20.15)', 'P(X=3) = ', dchyper(3, 10, n, m)
  print '(a,f20.15)', 'P(X<=3) = ', pchyper(3, 10, n, m)
  print '(a,i0)', 'q(0.9) = ', qchyper(0.9_dp, 10, n, m)
  print '(a,i0)', 'MLE(s) = ', mle_s(observed, n, m)
end program basic
