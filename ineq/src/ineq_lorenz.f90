! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module ineq_lorenz
   !! Empirical and theoretical Lorenz curves, asymmetry, and majorization.
   use ineq_common, only : dp, nan_value, prepare_sample, normalize_sample
   use r_ordering, only : r_order
   use r_special, only : r_regularized_beta
   use r_distributions, only : r_pnorm, r_qnorm
   use r_transforms, only : r_log1p, r_expm1
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   implicit none
   private
   integer, parameter, public :: ineq_ok = 0, ineq_invalid = 1, ineq_incomparable = 2
   type, public :: lorenz_curve
      !! Numeric counterpart of the R Lc object; no plotting or S3 metadata.
      real(dp), allocatable :: p(:) !! Population shares, including the origin.
      real(dp), allocatable :: l(:) !! Cumulative income shares.
      real(dp), allocatable :: l_general(:) !! Cumulative income per population; absent for Mehran bounds.
      integer :: status = ineq_invalid !! Zero on success; one for unsupported or undefined inputs.
   end type lorenz_curve
   public :: dp, lc, lc_mehran, lasym, major
   public :: theor_lc, lc_exp, lc_lognorm, lc_pareto, lc_singh, lc_dagum
contains
   pure function lc(x, n) result(curve)
      real(dp), intent(in) :: x(:) !! Class means or individual incomes; NaNs and corresponding frequencies are removed.
      real(dp), intent(in), optional :: n(:) !! Finite frequencies of length size(x), default one per observation.
      type(lorenz_curve) :: curve
      real(dp), allocatable :: y(:), w(:), income(:)
      integer, allocatable :: order(:)
      real(dp) :: sx, sw, total_w, total_x, running_w, running_x
      integer :: i, m
      logical :: valid
      allocate(curve%p(0), curve%l(0), curve%l_general(0))
      call paired_sample(x, n, y, w, valid)
      if (.not. valid .or. size(y) == 0) return
      call r_order(y, order)
      y = y(order)
      w = w(order)
      sx = maxval(abs(y))
      sw = maxval(abs(w))
      if (sx <= huge(1.0_dp) / real(size(y), dp) / max(1.0_dp, sw) .and. &
          sw <= huge(1.0_dp) / real(size(y), dp) .and. &
          sx >= tiny(1.0_dp) * real(size(y), dp) / max(sw, tiny(1.0_dp))) then
         ! Keep ordinary input arithmetic identical to the R source; Mehran slopes depend on ties.
         sx = 1.0_dp
         sw = 1.0_dp
      else
         call normalize_sample(y, sx)
         call normalize_sample(w, sw)
      end if
      total_w = sum(w)
      income = y * w
      total_x = sum(income)
      if (total_w == 0.0_dp .or. total_x == 0.0_dp) return
      m = size(y)
      curve%p = [(0.0_dp, i = 0, m)]
      curve%l = curve%p
      curve%l_general = curve%p
      running_w = 0.0_dp
      running_x = 0.0_dp
      do i = 1, m
         running_w = running_w + w(i)
         running_x = running_x + income(i)
         curve%p(i + 1) = running_w / total_w
         curve%l(i + 1) = running_x / total_x
         curve%l_general(i + 1) = (running_x / total_w) * sx
      end do
      curve%status = ineq_ok
   end function lc

   pure function lc_mehran(x, n) result(curve)
      real(dp), intent(in) :: x(:) !! At least two class means; NaNs and their matching frequencies are removed.
      real(dp), intent(in) :: n(:) !! Class frequencies of length size(x); zero-width classes can make bounds undefined.
      type(lorenz_curve) :: curve
      type(lorenz_curve) :: lower
      real(dp), allocatable :: q(:), ordinate(:), beta1(:), beta2(:), beta(:), d(:), p(:), l(:)
      integer :: k, i
      allocate(curve%p(0), curve%l(0), curve%l_general(0))
      lower = lc(x, n)
      if (lower%status /= ineq_ok) return
      k = size(lower%p)
      if (k < 3) return
      allocate(q(k + 1), ordinate(k + 1), beta1(k - 1), beta2(k - 2), beta(k - 2), d(k - 2))
      q = 1.0_dp
      q(1) = 0.0_dp
      ordinate = 0.0_dp
      ordinate(k + 1) = 1.0_dp
      do i = k, 2, -1
         q(i) = 2.0_dp * lower%p(i) - q(i + 1)
      end do
      do i = 2, k
         ordinate(i) = 2.0_dp * lower%l(i - 1) - ordinate(i - 1)
      end do
      beta1 = (lower%l(2:k) - lower%l(1:k - 1)) / (lower%p(2:k) - lower%p(1:k - 1))
      beta2 = (ordinate(3:k) - ordinate(2:k - 1)) / (q(3:k) - q(2:k - 1))
      if (any(ieee_is_nan(beta1)) .or. any(ieee_is_nan(beta2))) return
      do i = 1, k - 2
         if (beta1(i) > beta2(i)) then
            beta(i) = beta1(i)
         else if (beta2(i) > beta1(i + 1)) then
            beta(i) = beta1(i + 1)
         else
            beta(i) = beta2(i)
         end if
      end do
      d = lower%l(2:k - 1) - beta * lower%p(2:k - 1)
      deallocate(q, ordinate)
      allocate(q(k - 2))
      do i = 1, k - 3
         q(i) = (d(i + 1) - d(i)) / (beta(i) - beta(i + 1))
      end do
      q(k - 2) = 1.0_dp
      ordinate = beta * q + d
      p = [0.0_dp, -d(1) / beta(1), q, 1.0_dp]
      l = [0.0_dp, 0.0_dp, ordinate, 1.0_dp]
      curve%p = pack(p, ieee_is_finite(p))
      curve%l = pack(l, ieee_is_finite(p))
      curve%status = ineq_ok
   end function lc_mehran

   pure function lasym(x, n, interval, na_rm) result(value)
      real(dp), intent(in) :: x(:) !! Finite incomes; NaNs are optionally omitted before sorting.
      real(dp), intent(in), optional :: n(:) !! Frequencies; with missing x, upstream uses the first retained-count weights.
      logical, intent(in), optional :: interval !! Return two bounds when a value equals the mean; defaults to false.
      logical, intent(in), optional :: na_rm !! Remove NaNs; defaults to true; false returns one or two NaNs.
      real(dp), allocatable :: value(:), y(:), w(:)
      integer, allocatable :: order(:)
      real(dp) :: sx, sw, mu, m, total_w, lm, ln, a, lma, delta, xm, xm1
      logical :: valid, bounds
      bounds = .false.
      if (present(interval)) bounds = interval
      value = [nan_value()]
      if (bounds) value = [nan_value(), nan_value()]
      call prepare_sample(x, y, valid, na_rm)
      if (.not. valid .or. size(y) == 0) return
      ! Match Lasym's historical na.omit(x); n[o] sequence, including its weight alignment.
      allocate(w(size(y)))
      w = 1.0_dp
      if (present(n)) then
         if (size(n) /= size(x)) return
         w = n(1:size(y))
      end if
      if (.not. all(ieee_is_finite(w))) return
      call r_order(y, order)
      y = y(order)
      w = w(order)
      ! Do not rescale frequencies: the upstream interpolation adds one unit, not one class weight.
      sw = maxval(abs(w))
      if (sw == 0.0_dp) return
      total_w = sum(w)
      if (total_w == 0.0_dp .or. .not. ieee_is_finite(total_w)) return
      if (maxval(abs(y)) <= huge(1.0_dp) / real(size(y), dp) / max(1.0_dp, sw)) then
         ! Preserve exact mean ties on ordinary inputs before dividing by the income scale.
         mu = sum(y * w) / total_w
         call normalize_sample(y, sx)
         mu = mu / sx
      else
         call normalize_sample(y, sx)
         mu = sum(y * (w / total_w))
      end if
      m = sum(w, mask=y < mu)
      lm = sum(w * y, mask=y < mu)
      ln = sum(w * y)
      if (ln == 0.0_dp) return
      if (any(y == mu)) then
         a = sum(w, mask=y == mu)
         lma = sum(w * y, mask=y <= mu)
         value = [m / total_w + lm / ln, (m + a) / total_w + lma / ln]
         if (.not. bounds) value = [sum(value) / 2.0_dp]
      else
         if (.not. any(y < mu) .or. .not. any(y > mu)) return
         xm = maxval(y, mask=y < mu)
         xm1 = minval(y, mask=y > mu)
         delta = (mu - xm) / (xm1 - xm)
         value = [(m + delta) / total_w + (lm + delta * xm1) / ln]
      end if
   end function lasym

   pure subroutine major(x, y, dominates, status)
      real(dp), intent(in) :: x(:) !! First vector; NaNs are removed as by R sort's default.
      real(dp), intent(in) :: y(:) !! Second vector; retained lengths and exact sums must agree.
      logical, intent(out) :: dominates !! True when every ascending cumulative sum of x is at most that of y.
      integer, intent(out), optional :: status !! Zero when comparable; two for unequal lengths or sums; one for infinities.
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: sum_a, sum_b, sx, sy, scale
      integer :: i
      logical :: valid_a, valid_b
      dominates = .false.
      if (present(status)) status = ineq_invalid
      call prepare_sample(x, a, valid_a, sort_values=.true.)
      call prepare_sample(y, b, valid_b, sort_values=.true.)
      if (.not. valid_a .or. .not. valid_b) return
      if (present(status)) status = ineq_incomparable
      if (size(a) /= size(b)) return
      if (size(a) == 0) then
         dominates = .true.
         if (present(status)) status = ineq_ok
         return
      end if
      scale = max(maxval(abs(a)), maxval(abs(b)))
      if (scale == 0.0_dp) scale = 1.0_dp
      ! Exact totals use original arithmetic when representable, preserving R's no-tolerance rule.
      sum_a = sum(a)
      sum_b = sum(b)
      if (ieee_is_finite(sum_a) .and. ieee_is_finite(sum_b)) then
         if (sum_a /= sum_b) return
      else
         if (sum(a / scale) /= sum(b / scale)) return
      end if
      sx = 0.0_dp
      sy = 0.0_dp
      dominates = .true.
      if (scale <= huge(1.0_dp) / real(size(a), dp)) scale = 1.0_dp
      do i = 1, size(a) - 1
         sx = sx + a(i) / scale
         sy = sy + b(i) / scale
         if (sx > sy) dominates = .false.
      end do
      if (present(status)) status = ineq_ok
   end subroutine major

   pure elemental real(dp) function lc_exp(p) result(value)
      real(dp), intent(in) :: p !! Population share in [0,1]; the removable endpoint singularity at one is filled with one.
      real(dp) :: term, contribution
      integer :: j
      value = nan_value()
      if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      if (p == 1.0_dp) then
         value = 1.0_dp
      else if (p < 1.0e-4_dp) then
         ! Positive series p**2/2 + p**3/6 + ... avoids subtracting nearly equal numbers.
         term = p * p
         value = 0.5_dp * term
         do j = 3, 32
            term = term * p
            contribution = term / real(j * (j - 1), dp)
            value = value + contribution
            if (contribution <= epsilon(value) * value) exit
         end do
      else
         value = p + (1.0_dp - p) * r_log1p(-p)
      end if
   end function lc_exp

   pure elemental real(dp) function lc_lognorm(p, parameter) result(value)
      real(dp), intent(in) :: p !! Population share in [0,1]; endpoints return zero and one.
      real(dp), intent(in), optional :: parameter !! Positive log-income spread, default one; nonpositive values use one.
      real(dp) :: sigma
      value = nan_value()
      if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      sigma = 1.0_dp
      if (present(parameter)) then
         if (ieee_is_nan(parameter)) return
         if (parameter > 0.0_dp) sigma = parameter
      end if
      if (p == 0.0_dp .or. p == 1.0_dp) then
         value = p
      else
         value = r_pnorm(r_qnorm(p) - sigma)
      end if
   end function lc_lognorm

   pure elemental real(dp) function lc_pareto(p, parameter) result(value)
      real(dp), intent(in) :: p !! Population share in [0,1].
      real(dp), intent(in), optional :: parameter !! Pareto shape greater than one, default two; inadmissible finite values use two.
      real(dp) :: shape, k
      value = nan_value()
      if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      shape = 2.0_dp
      if (present(parameter)) then
         if (.not. ieee_is_finite(parameter)) return
         if (parameter > 1.0_dp) shape = parameter
      end if
      k = (shape - 1.0_dp) / shape
      if (p == 1.0_dp) then
         value = 1.0_dp
      else
         value = -r_expm1(k * r_log1p(-p))
      end if
   end function lc_pareto

   pure real(dp) function lc_singh(p, parameter) result(value)
      real(dp), intent(in) :: p !! Population share in [0,1].
      real(dp), intent(in), optional :: parameter(:) !! Two parameters, default [2,2]; inadmissible values use the default.
      real(dp) :: b, d, z
      value = nan_value()
      if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      b = 1.0_dp
      d = 1.0_dp / 3.0_dp
      if (present(parameter)) then
         if (size(parameter) >= 2) then
            if (all(ieee_is_finite(parameter(1:2)))) then
               if (parameter(1) > 0.0_dp .and. parameter(1) < 2.0_dp + parameter(2)) then
                  b = parameter(1) - 1.0_dp
                  d = 1.0_dp / (parameter(2) + 1.0_dp)
               end if
            end if
         end if
      end if
      ! The original guard is weaker than positivity of the resulting beta shapes.
      if (.not. ieee_is_finite(d)) return
      if (p == 1.0_dp) then
         z = 1.0_dp
      else
         z = -r_expm1(b * r_log1p(-p))
      end if
      value = lorenz_beta_cdf(z, 1.0_dp + d, b - d)
   end function lc_singh

   pure real(dp) function lc_dagum(p, parameter) result(value)
      real(dp), intent(in) :: p !! Population share in [0,1].
      real(dp), intent(in), optional :: parameter(:) !! Two parameters, default [2,2]; first parameter at most one uses the default.
      real(dp) :: a, b
      value = nan_value()
      if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      a = 0.5_dp
      b = 2.0_dp
      if (present(parameter)) then
         if (size(parameter) >= 2) then
            if (all(ieee_is_finite(parameter(1:2)))) then
               if (parameter(1) > 1.0_dp) then
                  a = 1.0_dp / parameter(1)
                  b = parameter(2)
               end if
            end if
         end if
      end if
      if (b == 0.0_dp) then
         value = 1.0_dp
         return
      end if
      value = lorenz_beta_cdf(p ** b, a + 1.0_dp / b, 1.0_dp - a)
   end function lc_dagum

   pure real(dp) function theor_lc(p, type, parameter) result(value)
      real(dp), intent(in) :: p !! Population share at which to evaluate the selected theoretical Lorenz curve.
      character(*), intent(in), optional :: type !! Exact R distribution name; default Singh-Maddala, without partial matching.
      real(dp), intent(in), optional :: parameter(:) !! Distribution parameters; omitted values use the upstream defaults.
      character(:), allocatable :: method
      real(dp) :: a
      method = 'Singh-Maddala'
      if (present(type)) method = type
      a = 0.0_dp
      if (present(parameter)) then
         if (size(parameter) > 0) a = parameter(1)
      end if
      select case (method)
      case ('Singh-Maddala')
         value = lc_singh(p, parameter)
      case ('Dagum')
         value = lc_dagum(p, parameter)
      case ('lognorm')
         value = lc_lognorm(p, a)
      case ('Pareto')
         value = lc_pareto(p, a)
      case ('exponential')
         value = lc_exp(p)
      case default
         value = nan_value()
      end select
   end function theor_lc

   pure elemental real(dp) function lorenz_beta_cdf(z, a, b) result(value)
      real(dp), intent(in) :: z !! Beta-CDF argument, possibly outside [0,1] for upstream-accepted parameters.
      real(dp), intent(in) :: a !! First beta shape; zero denotes the R degenerate limit.
      real(dp), intent(in) :: b !! Second beta shape; zero denotes the R degenerate limit.
      value = nan_value()
      if (ieee_is_nan(z) .or. .not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) return
      if (a < 0.0_dp .or. b < 0.0_dp) return
      if (a == 0.0_dp) then
         value = 1.0_dp
         if (b == 0.0_dp .and. z < 1.0_dp) value = 0.5_dp
      else if (z >= 1.0_dp) then
         value = 1.0_dp
      else if (b == 0.0_dp .or. z <= 0.0_dp) then
         value = 0.0_dp
      else
         value = r_regularized_beta(z, a, b)
      end if
   end function lorenz_beta_cdf

   pure subroutine paired_sample(x, n, y, w, valid)
      real(dp), intent(in) :: x(:) !! Original class means or observations, with NaNs removable.
      real(dp), intent(in), optional :: n(:) !! Optional matching class frequencies.
      real(dp), allocatable, intent(out) :: y(:) !! Observations after paired removal of missing x.
      real(dp), allocatable, intent(out) :: w(:) !! Matching frequencies, or unit frequencies by default.
      logical, intent(out) :: valid !! False for a shape mismatch or nonfinite retained values or frequencies.
      call prepare_sample(x, y, valid)
      allocate(w(size(y)))
      w = 1.0_dp
      if (present(n)) then
         if (size(n) /= size(x)) then
            valid = .false.
            return
         end if
         w = pack(n, .not. ieee_is_nan(x))
      end if
      valid = valid .and. all(ieee_is_finite(w))
   end subroutine paired_sample
end module ineq_lorenz
