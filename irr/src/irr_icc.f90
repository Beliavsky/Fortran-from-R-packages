! SPDX-License-Identifier: GPL-2.0-or-later
module irr_icc
   use irr_common
   use r_distributions, only: r_pf, r_qf
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   implicit none
   private
   public :: icc
contains
   pure function icc(ratings, model, kind, unit, r0, conf_level) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater measurements; incomplete subjects are removed listwise.
      character(*), intent(in), optional :: model !! 'oneway' (default) or 'twoway' ANOVA model.
      character(*), intent(in), optional :: kind !! 'consistency' (default) or 'agreement'; ignored by the one-way model.
      character(*), intent(in), optional :: unit !! 'single' (default) or 'average' score reliability.
      real(dp), intent(in), optional :: r0 !! Null reliability for the one-sided F test, default zero; must be in (-1,1).
      real(dp), intent(in), optional :: conf_level !! Two-sided confidence level strictly between zero and one, default 0.95.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :)
      real(dp) :: n, r, msr, msw, msc, mse, error, rho, level, a, b, v, fl, fu, ratio, factor
      logical :: average, agreement
      character(:), allocatable :: model_used, type_used, unit_used
      call complete_ratings(ratings, x, res, 2)
      if (res%status /= irr_ok) return
      model_used = 'oneway'
      type_used = 'consistency'
      unit_used = 'single'
      if (present(model)) model_used = model
      if (present(kind)) type_used = kind
      if (present(unit)) unit_used = unit
      rho = 0.0_dp
      level = 0.95_dp
      if (present(r0)) rho = r0
      if (present(conf_level)) level = conf_level
      res%status = irr_invalid
      if (model_used /= 'oneway' .and. model_used /= 'twoway') return
      if (type_used /= 'consistency' .and. type_used /= 'agreement') return
      if (unit_used /= 'single' .and. unit_used /= 'average') return
      if (.not. ieee_is_finite(rho) .or. .not. ieee_is_finite(level)) return
      if (rho <= -1.0_dp .or. rho >= 1.0_dp .or. level <= 0.0_dp .or. level >= 1.0_dp) return
      n = real(res%subjects, dp)
      r = real(res%raters, dp)
      average = unit_used == 'average'
      agreement = model_used == 'twoway' .and. type_used == 'agreement'
      if (.not. average .and. 1.0_dp + (r - 1.0_dp) * rho <= 0.0_dp) return
      call mean_squares(x, msr, msw, msc, mse)
      res%status = irr_ok
      res%df1 = n - 1.0_dp
      if (.not. agreement) then
         error = mse
         res%df2 = (n - 1.0_dp) * (r - 1.0_dp)
         if (model_used == 'oneway') then
            error = msw
            res%df2 = n * (r - 1.0_dp)
         end if
         factor = 1.0_dp - rho
         if (average) then
            res%value = divide(msr - error, msr)
         else
            res%value = divide(msr - error, msr + (r - 1.0_dp) * error)
            factor = factor / (1.0_dp + (r - 1.0_dp) * rho)
         end if
         ratio = divide(msr, error)
         res%statistic = ratio * factor
         fl = ratio / r_qf((1.0_dp + level) / 2.0_dp, res%df1, res%df2)
         fu = ratio * r_qf((1.0_dp + level) / 2.0_dp, res%df2, res%df1)
         if (average) then
            res%lbound = 1.0_dp - divide(1.0_dp, fl)
            res%ubound = 1.0_dp - divide(1.0_dp, fu)
         else
            res%lbound = divide(fl - 1.0_dp, fl + r - 1.0_dp)
            res%ubound = divide(fu - 1.0_dp, fu + r - 1.0_dp)
         end if
         if (error == 0.0_dp .and. msr > 0.0_dp) then
            res%lbound = 1.0_dp
            res%ubound = 1.0_dp
         end if
      else
         if (average) then
            res%value = divide(msr - mse, msr + (msc - mse) / n)
            a = rho / (n * (1.0_dp - rho))
            b = 1.0_dp + rho * (n - 1.0_dp) / (n * (1.0_dp - rho))
         else
            res%value = divide(msr - mse, msr + (r - 1.0_dp) * mse + r * (msc - mse) / n)
            a = r * rho / (n * (1.0_dp - rho))
            b = 1.0_dp + r * rho * (n - 1.0_dp) / (n * (1.0_dp - rho))
         end if
         res%df2 = degrees(a, b, msc, mse, n, r)
         res%statistic = divide(msr, a * msc + b * mse)
         if (res%value == 1.0_dp .and. mse == 0.0_dp .and. msc == 0.0_dp) then
            res%lbound = 1.0_dp
            res%ubound = 1.0_dp
         else
            a = divide(r * res%value, n * (1.0_dp - res%value))
            b = 1.0_dp + divide(r * res%value * (n - 1.0_dp), n * (1.0_dp - res%value))
            v = degrees(a, b, msc, mse, n, r)
            fl = r_qf((1.0_dp + level) / 2.0_dp, n - 1.0_dp, v)
            fu = r_qf((1.0_dp + level) / 2.0_dp, v, n - 1.0_dp)
            if (average) then
               res%lbound = divide(n * (msr - fl * mse), fl * (msc - mse) + n * msr)
               res%ubound = divide(n * (fu * msr - mse), msc - mse + n * fu * msr)
            else
               res%lbound = divide(n * (msr - fl * mse), &
                  fl * (r * msc + (r * n - r - n) * mse) + n * msr)
               res%ubound = divide(n * (fu * msr - mse), &
                  r * msc + (r * n - r - n) * mse + n * fu * msr)
            end if
         end if
      end if
      res%p_value = r_pf(res%statistic, res%df1, res%df2, lower_tail = .false.)
      if (ieee_is_nan(res%value)) res%status = irr_degenerate
   end function

   pure function degrees(a, b, msc, mse, n, r) result(v)
      real(dp), intent(in) :: a !! Coefficient of the rater mean square in the Satterthwaite approximation.
      real(dp), intent(in) :: b !! Coefficient of the residual mean square in the Satterthwaite approximation.
      real(dp), intent(in) :: msc !! Rater mean square from the two-way ANOVA.
      real(dp), intent(in) :: mse !! Residual mean square from the two-way ANOVA.
      real(dp), intent(in) :: n !! Number of subjects, at least two.
      real(dp), intent(in) :: r !! Number of raters, at least two.
      real(dp) :: v
      v = divide((a * msc + b * mse) ** 2, (a * msc) ** 2 / (r - 1.0_dp) + &
         (b * mse) ** 2 / ((n - 1.0_dp) * (r - 1.0_dp)))
   end function
end module
