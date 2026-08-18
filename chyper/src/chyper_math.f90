! SPDX-License-Identifier: MIT
module chyper_math
  use chyper_kinds, only : dp
  implicit none
  private
  public :: log_choose, hypergeom_pmf
contains
  pure real(dp) function log_choose(n, k) result(ans)
    integer, intent(in) :: n, k
    if (k < 0 .or. k > n .or. n < 0) then
      ans = -huge(1.0_dp)
    else
      ans = log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) &
          - log_gamma(real(n - k + 1, dp))
    end if
  end function log_choose

  pure real(dp) function hypergeom_pmf(x, good, bad, draws) result(ans)
    integer, intent(in) :: x, good, bad, draws
    real(dp) :: lp
    if (good < 0 .or. bad < 0 .or. draws < 0 .or. draws > good + bad) then
      ans = 0.0_dp
      return
    end if
    if (x < 0 .or. x > good .or. draws - x < 0 .or. draws - x > bad) then
      ans = 0.0_dp
      return
    end if
    lp = log_choose(good, x) + log_choose(bad, draws - x) &
       - log_choose(good + bad, draws)
    if (lp < log(tiny(1.0_dp))) then
      ans = 0.0_dp
    else
      ans = exp(lp)
    end if
  end function hypergeom_pmf
end module chyper_math
