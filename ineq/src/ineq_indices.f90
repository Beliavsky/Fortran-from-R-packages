! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module ineq_indices
   !! Inequality, concentration, and poverty measures from ineq 0.2-13.
   use ineq_common, only : dp, prepare_sample, normalize_sample, nan_value
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   use r_stability, only : r_log_mean_exp, r_log_sum_exp
   use r_transforms, only : r_expm1, r_log1p
   implicit none
   private
   public :: dp, ineq, gini, rs, atkinson, theil, kolm, var_coeff, entropy
   public :: conc, herfindahl, rosenbluth, pov, watts, sen, sst, foster
contains
   pure real(dp) function ineq(x, parameter, type, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Observations; NaNs are removed by default.
      real(dp), intent(in), optional :: parameter !! Measure parameter; defaults depend on the selected measure.
      character(*), intent(in), optional :: type !! Exact R measure name, default Gini; unknown names return NaN.
      logical, intent(in), optional :: na_rm !! Remove missing observations; defaults to true.
      character(:), allocatable :: method
      method = 'Gini'
      if (present(type)) method = type
      select case (method)
      case ('Gini')
         value = gini(x, na_rm=na_rm)
      case ('RS')
         value = rs(x, na_rm)
      case ('Atkinson')
         value = atkinson(x, parameter, na_rm)
      case ('Theil')
         value = theil(x, parameter, na_rm)
      case ('Kolm')
         value = kolm(x, parameter, na_rm)
      case ('var')
         value = var_coeff(x, na_rm=na_rm)
      case ('square.var')
         value = var_coeff(x, square=.true., na_rm=na_rm)
      case ('entropy')
         value = entropy(x, parameter, na_rm)
      case default
         value = nan_value()
      end select
   end function ineq

   pure real(dp) function gini(x, corr, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Observations, normally nonnegative; signed finite values follow the R formula.
      logical, intent(in), optional :: corr !! Apply the n/(n-1) finite-sample correction; defaults to false.
      logical, intent(in), optional :: na_rm !! Remove NaNs before calculation; defaults to true.
      real(dp), allocatable :: y(:)
      real(dp) :: scale, total, n
      integer :: i
      logical :: valid, correction
      call prepare_sample(x, y, valid, na_rm, .true.)
      value = nan_value()
      if (.not. valid .or. size(y) == 0) return
      correction = .false.
      if (present(corr)) correction = corr
      if (correction .and. size(y) == 1) return
      call normalize_sample(y, scale)
      total = sum(y)
      n = real(size(y), dp)
      value = sum([(y(i) * (2.0_dp * real(i, dp) - n - 1.0_dp), i = 1, size(y))]) / total
      if (correction) then
         value = value / (n - 1.0_dp)
      else
         value = value / n
      end if
   end function gini

   pure real(dp) function rs(x, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Observations used in the relative mean absolute deviation.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:)
      real(dp) :: scale, mu
      logical :: valid
      call prepare_sample(x, y, valid, na_rm)
      value = nan_value()
      if (.not. valid .or. size(y) == 0) return
      call normalize_sample(y, scale)
      mu = sum(y) / real(size(y), dp)
      value = sum(abs(y - mu)) / (2.0_dp * sum(y))
   end function rs

   pure real(dp) function atkinson(x, parameter, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Income observations; logarithmic branches require nonnegative values.
      real(dp), intent(in), optional :: parameter !! Aversion parameter, default 0.5; one selects the geometric mean.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), logs(:)
      real(dp) :: a, scale, mu, power, center
      logical :: valid
      value = nan_value()
      a = 0.5_dp
      if (present(parameter)) a = parameter
      if (.not. ieee_is_finite(a)) return
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. size(y) == 0) return
      if (a == 1.0_dp .and. any(y < 0.0_dp)) return
      if (all(y >= 0.0_dp) .or. all(y <= 0.0_dp)) logs = log(abs(y))
      call normalize_sample(y, scale)
      mu = sum(y) / real(size(y), dp)
      if (mu == 0.0_dp) return
      if (allocated(logs)) logs = (logs - log(scale)) - log(abs(mu))
      if (a == 1.0_dp) then
         value = -r_expm1(sum(logs) / real(size(y), dp))
      else
         y = y / mu
         power = 1.0_dp - a
         if (allocated(logs)) then
            if (abs(power) * maxval(abs(logs)) < 0.5_dp) then
               center = sum(logs) / real(size(y), dp)
               value = -r_expm1(center + r_log1p(sum(r_expm1(power * (logs - center))) / &
                  real(size(y), dp)) / power)
            else
               value = -r_expm1(r_log_mean_exp(power * logs) / power)
            end if
         else
            value = 1.0_dp - (sum(y ** power) / real(size(y), dp)) ** (1.0_dp / power)
         end if
      end if
   end function atkinson

   pure real(dp) function var_coeff(x, square, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Observations; fewer than two retained values yield NaN as in R var.
      logical, intent(in), optional :: square !! Return the squared coefficient when true; defaults to false.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:)
      real(dp) :: scale, mu
      logical :: valid
      value = nan_value()
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. size(y) < 2) return
      call normalize_sample(y, scale)
      mu = sum(y) / real(size(y), dp)
      value = sqrt(sum((y - mu) ** 2) / real(size(y), dp)) / mu
      if (present(square)) then
         if (square) value = value ** 2
      end if
   end function var_coeff

   pure real(dp) function theil(x, parameter, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Observations; the parameter-zero branch removes zero incomes before its mean.
      real(dp), intent(in), optional :: parameter !! Zero selects Theil T; every other finite value selects mean log deviation.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), logs(:)
      real(dp) :: a, scale, mu
      logical :: valid
      value = nan_value()
      a = 0.0_dp
      if (present(parameter)) a = parameter
      if (.not. ieee_is_finite(a)) return
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid) return
      if (a == 0.0_dp) y = pack(y, y /= 0.0_dp)
      if (size(y) == 0) return
      if (.not. (all(y >= 0.0_dp) .or. all(y <= 0.0_dp))) return
      if (a /= 0.0_dp .and. any(y < 0.0_dp)) return
      logs = log(abs(y))
      call normalize_sample(y, scale)
      mu = sum(y) / real(size(y), dp)
      if (mu == 0.0_dp) return
      logs = (logs - log(scale)) - log(abs(mu))
      if (a == 0.0_dp) then
         value = sum(y * logs) / sum(y)
      else
         value = -sum(logs) / real(size(y), dp)
      end if
   end function theil

   pure real(dp) function kolm(x, parameter, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Finite observations; the Kolm index has the same units as x.
      real(dp), intent(in), optional :: parameter !! Nonzero absolute-aversion coefficient, default one; zero yields NaN.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:)
      real(dp) :: a, scale, mu, anchor
      logical :: valid
      value = nan_value()
      a = 1.0_dp
      if (present(parameter)) a = parameter
      if (.not. ieee_is_finite(a) .or. a == 0.0_dp) return
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. size(y) == 0) return
      call normalize_sample(y, scale)
      mu = sum(y) / real(size(y), dp)
      if (a > 0.0_dp) then
         anchor = minval(y)
      else
         anchor = maxval(y)
      end if
      if (abs(a) * maxval(abs(mu - y)) < 0.5_dp / scale) then
         value = r_log1p(sum(r_expm1((a * (mu - y)) * scale)) / real(size(y), dp)) / a
         return
      end if
      ! Center before exponentiation to avoid the upstream exponential overflow.
      value = (mu - anchor) * scale + r_log_mean_exp((-a * (y - anchor)) * scale) / a
   end function kolm

   pure real(dp) function entropy(x, parameter, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Observations for generalized entropy; negative values may make powers undefined.
      real(dp), intent(in), optional :: parameter !! Entropy order, default 0.5; zero and one use the upstream Theil branches.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), logs(:)
      real(dp) :: a, scale, mu
      logical :: valid
      a = 0.5_dp
      if (present(parameter)) a = parameter
      value = nan_value()
      if (.not. ieee_is_finite(a)) return
      if (a == 0.0_dp) then
         value = theil(x, 1.0_dp, na_rm)
         return
      else if (a == 1.0_dp) then
         value = theil(x, 0.0_dp, na_rm)
         return
      end if
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. size(y) == 0) return
      if (all(y >= 0.0_dp) .or. all(y <= 0.0_dp)) logs = log(abs(y))
      call normalize_sample(y, scale)
      mu = sum(y) / real(size(y), dp)
      if (mu == 0.0_dp) return
      y = y / mu
      if (allocated(logs)) then
         logs = (logs - log(scale)) - log(abs(mu))
         if (abs(a - 1.0_dp) < 0.01_dp) then
            ! The normalized mean is analytically one; remove its rounding residual near order one.
            logs = pack(logs, y > 0.0_dp)
            y = pack(y, y > 0.0_dp)
            value = sum(y * r_expm1((a - 1.0_dp) * logs)) / real(size(x) - count(ieee_is_nan(x)), dp) &
               / a / (a - 1.0_dp)
         else
            value = sum(r_expm1(a * logs)) / real(size(y), dp) / a / (a - 1.0_dp)
         end if
      else
         value = sum(y ** a - 1.0_dp) / real(size(y), dp) / a / (a - 1.0_dp)
      end if
   end function entropy

   pure real(dp) function conc(x, parameter, type, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Firm sizes or other concentration observations.
      real(dp), intent(in), optional :: parameter !! Herfindahl order parameter; defaults to one.
      character(*), intent(in), optional :: type !! Exact Herfindahl or Rosenbluth name; default Herfindahl.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      character(:), allocatable :: method
      method = 'Herfindahl'
      if (present(type)) method = type
      select case (method)
      case ('Herfindahl')
         value = herfindahl(x, parameter, na_rm)
      case ('Rosenbluth')
         value = rosenbluth(x, na_rm)
      case default
         value = nan_value()
      end select
   end function conc

   pure real(dp) function herfindahl(x, parameter, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Sizes normalized to shares of their sum.
      real(dp), intent(in), optional :: parameter !! Order minus one, default one; zero reproduces the R infinite-power limit.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), logs(:)
      real(dp) :: a, scale, total
      logical :: valid
      value = nan_value()
      a = 1.0_dp
      if (present(parameter)) a = parameter
      if (.not. ieee_is_finite(a)) return
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid) return
      if (size(y) == 0) then
         value = 0.0_dp
         if (a < 0.0_dp) value = value ** (1.0_dp / a)
         return
      end if
      if (all(y >= 0.0_dp) .or. all(y <= 0.0_dp)) logs = log(abs(y))
      call normalize_sample(y, scale)
      total = sum(y)
      if (total == 0.0_dp) return
      if (a == 0.0_dp) then
         value = 1.0_dp
         return
      end if
      y = y / total
      if (allocated(logs)) then
         logs = (logs - log(scale)) - log(abs(total))
         if (a == -1.0_dp) then
            value = 1.0_dp / real(size(y), dp)
         else if (abs(a) < 0.01_dp) then
            logs = pack(logs, y > 0.0_dp)
            y = pack(y, y > 0.0_dp)
            value = exp(r_log1p(sum(y * r_expm1(a * logs))) / a)
         else
            value = exp(r_log_sum_exp((a + 1.0_dp) * logs) / a)
         end if
      else
         value = sum(y ** (a + 1.0_dp)) ** (1.0_dp / a)
      end if
   end function herfindahl

   pure real(dp) function rosenbluth(x, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Sizes ranked from smallest to largest in the upstream concentration formula.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:)
      real(dp) :: scale, total
      integer :: i
      logical :: valid
      value = nan_value()
      call prepare_sample(x, y, valid, na_rm, .true.)
      if (.not. valid) return
      if (size(y) == 0) then
         value = -1.0_dp
         return
      end if
      call normalize_sample(y, scale)
      total = sum(y)
      if (total == 0.0_dp) return
      value = 1.0_dp / (2.0_dp * sum([(real(size(y) - i + 1, dp) * y(i), i = 1, size(y))]) / total - 1.0_dp)
   end function rosenbluth

   pure real(dp) function pov(x, k, parameter, type, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Income observations.
      real(dp), intent(in) :: k !! Finite poverty threshold; observations strictly below k are poor.
      real(dp), intent(in), optional :: parameter !! Foster exponent plus one; default one.
      character(*), intent(in), optional :: type !! Exact Watts, Sen, SST, or Foster name; default Watts.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      character(:), allocatable :: method
      method = 'Watts'
      if (present(type)) method = type
      select case (method)
      case ('Watts')
         value = watts(x, k, na_rm)
      case ('Sen')
         value = sen(x, k, na_rm)
      case ('SST')
         value = sst(x, k, na_rm)
      case ('Foster')
         value = foster(x, k, parameter, na_rm)
      case default
         value = nan_value()
      end select
   end function pov

   pure real(dp) function watts(x, k, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Income observations; zero poor incomes produce positive infinity for positive k.
      real(dp), intent(in) :: k !! Poverty line; poverty uses the strict comparison x < k.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), poor(:)
      logical :: valid
      value = nan_value()
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. .not. ieee_is_finite(k)) return
      poor = pack(y, y < k)
      value = 0.0_dp
      if (size(poor) == 0) return
      if (k > 0.0_dp .and. all(poor >= 0.0_dp)) then
         value = sum(log(k) - log(poor)) / real(size(y), dp)
      else
         value = sum(log(k / poor)) / real(size(y), dp)
      end if
   end function watts

   pure real(dp) function sen(x, k, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Income observations; an all-zero poor subsample retains upstream NaN behavior.
      real(dp), intent(in) :: k !! Poverty threshold; values equal to k are excluded.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), poor(:)
      real(dp) :: headcount, gap
      logical :: valid
      value = nan_value()
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. .not. ieee_is_finite(k)) return
      poor = pack(y, y < k)
      value = 0.0_dp
      if (size(poor) == 0) return
      headcount = real(size(poor), dp) / real(size(y), dp)
      gap = sum(1.0_dp - poor / k) / real(size(poor), dp)
      value = headcount * (gap + (1.0_dp - gap) * gini(poor))
   end function sen

   pure real(dp) function sst(x, k, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Income observations, sorted internally for Sen-Shorrocks-Thon poverty.
      real(dp), intent(in) :: k !! Poverty line, with a strict less-than test.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), poor(:)
      real(dp) :: n
      integer :: i
      logical :: valid
      value = nan_value()
      call prepare_sample(x, y, valid, na_rm, .true.)
      if (.not. valid .or. .not. ieee_is_finite(k)) return
      poor = pack(y, y < k)
      value = 0.0_dp
      if (size(poor) == 0) return
      n = real(size(y), dp)
      value = sum([((2.0_dp * n - 2.0_dp * real(i, dp) + 1.0_dp) * &
         (1.0_dp - poor(i) / k), i = 1, size(poor))]) / n / n
   end function sst

   pure real(dp) function foster(x, k, parameter, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Income observations.
      real(dp), intent(in) :: k !! Poverty line; an empty poor subsample yields zero.
      real(dp), intent(in), optional :: parameter !! Upstream exponent plus one, default one (the headcount ratio).
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true.
      real(dp), allocatable :: y(:), poor(:)
      real(dp) :: a
      logical :: valid
      value = nan_value()
      a = 1.0_dp
      if (present(parameter)) a = parameter
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. .not. ieee_is_finite(k) .or. .not. ieee_is_finite(a)) return
      poor = pack(y, y < k)
      value = 0.0_dp
      if (size(poor) == 0) return
      value = sum((1.0_dp - poor / k) ** (a - 1.0_dp)) / real(size(y), dp)
   end function foster
end module ineq_indices
