! SPDX-License-Identifier: GPL-3.0-only
module spantest_probability
  use spantest_kinds, only : dp
  use r_distributions, only : core_pnorm => r_pnorm
  use r_distributions, only : core_qnorm => r_qnorm
  use r_distributions, only : core_pt => r_pt
  use r_distributions, only : core_pf => r_pf
  implicit none
  private
  public :: normal_cdf, normal_quantile, student_t_cdf, f_upper_tail

contains

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = core_pnorm(x)
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else
      x = core_qnorm(p)
    end if
  end function normal_quantile

  pure real(dp) function student_t_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (df <= 0.0_dp) then
      p = 0.0_dp
    else
      p = core_pt(x, df)
    end if
  end function student_t_cdf

  pure real(dp) function f_upper_tail(x, df1, df2) result(p)
    real(dp), intent(in) :: x, df1, df2
    if (x < 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
      p = 1.0_dp
    else
      p = core_pf(x, df1, df2, lower_tail=.false.)
    end if
  end function f_upper_tail

end module spantest_probability
