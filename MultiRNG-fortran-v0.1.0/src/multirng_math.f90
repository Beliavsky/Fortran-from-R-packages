! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multirng_math
  use multirng_kinds, only : dp
  implicit none
  private
  public :: normal_cdf
contains
  elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf
end module multirng_math
