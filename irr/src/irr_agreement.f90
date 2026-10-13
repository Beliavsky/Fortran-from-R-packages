! SPDX-License-Identifier: GPL-2.0-or-later
module irr_agreement
   use irr_common
   use r_descriptive, only: r_correlation, r_mean, r_variance
   use r_distributions, only: r_pchisq, r_pf, r_pnorm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_positive_inf
   implicit none
   private
   public :: agree, maxwell, meancor, meanrho, kendall, finn, robinson, iota, kripp_alpha
contains
   pure function agree(ratings, tolerance) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater numeric ratings; NaN rows are removed.
      real(dp), intent(in), optional :: tolerance !! Maximum within-subject range counted as agreement, default zero.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :)
      real(dp) :: tol
      call complete_ratings(ratings, x, res)
      if (res%status /= irr_ok) return
      tol = 0.0_dp
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) then
         res%status = irr_invalid
         return
      end if
      res%value = 100.0_dp * real(count(maxval(x, dim = 2) - minval(x, dim = 2) <= tol), dp) / real(res%subjects, dp)
   end function

   pure function maxwell(ratings) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Two numeric category columns, at most two distinct categories; NaN rows removed.
      type(irr_result) :: res
      real(dp), allocatable :: tab(:, :)
      integer :: i
      call contingency(ratings, tab, res)
      if (res%status /= irr_ok) return
      if (size(tab, 1) > 2) then
         res%status = irr_invalid
         return
      end if
      res%value = 0.0_dp
      do i = 1, size(tab, 1)
         res%value = res%value + tab(i, i)
      end do
      res%value = 2.0_dp * res%value / real(res%subjects, dp) - 1.0_dp
   end function

   pure function meancor(ratings, fisher) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Numeric subject-by-rater ratings; NaN rows removed listwise.
      logical, intent(in), optional :: fisher !! Average Fisher-transformed correlations and provide inference; default true.
      type(irr_result) :: res
      res = mean_correlation(ratings, fisher, .false.)
   end function

   pure function meanrho(ratings, fisher) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Numeric subject-by-rater ratings; tied values receive average ranks.
      logical, intent(in), optional :: fisher !! Average Fisher-transformed rank correlations; default true.
      type(irr_result) :: res
      res = mean_correlation(ratings, fisher, .true.)
   end function

   pure function mean_correlation(ratings, fisher, ranked) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater observations, with NaN listwise deletion.
      logical, intent(in), optional :: fisher !! Whether to exclude perfect correlations and average on Fisher scale.
      logical, intent(in) :: ranked !! Use average ranks for Spearman correlation when true.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :), ranks(:)
      real(dp) :: cor, total, ties
      integer :: i, j, n
      logical :: use_fisher
      call complete_ratings(ratings, x, res, 2)
      if (res%status /= irr_ok) return
      use_fisher = .true.
      if (present(fisher)) use_fisher = fisher
      if (ranked) then
         allocate(ranks(size(x, 1)))
         do j = 1, size(x, 2)
            call average_ranks(x(:, j), ranks, ties)
            x(:, j) = ranks
         end do
      end if
      total = 0.0_dp
      n = 0
      do i = 1, size(x, 2) - 1
         do j = i + 1, size(x, 2)
            cor = r_correlation(x(:, i), x(:, j))
            if (ieee_is_nan(cor)) then
               res%status = irr_degenerate
               return
            end if
            if (use_fisher) then
               if (abs(cor) >= 1.0_dp) then
                  res%dropped = res%dropped + 1
                  cycle
               end if
               cor = atanh(cor)
            end if
            n = n + 1
            total = total + cor
         end do
      end do
      if (n == 0) then
         res%status = irr_degenerate
         return
      end if
      res%value = total / real(n, dp)
      if (use_fisher) then
         res%value = tanh(res%value)
         if (res%subjects < 4) then
            res%status = irr_degenerate
            return
         end if
         call finish_z(res, 1.0_dp / real(res%subjects - 3, dp))
      end if
   end function

   pure function kendall(ratings, correct) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater ratings; subjects are ranked separately within each rater.
      logical, intent(in), optional :: correct !! Apply the tie correction to Kendall W, default false as in R.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :), ranks(:), sums(:)
      real(dp) :: ties, ties_total, n, r, denominator
      integer :: j
      logical :: adjusted
      call complete_ratings(ratings, x, res, 2)
      if (res%status /= irr_ok) return
      adjusted = .false.
      if (present(correct)) adjusted = correct
      allocate(ranks(size(x, 1)), sums(size(x, 1)))
      sums = 0.0_dp
      ties_total = 0.0_dp
      do j = 1, size(x, 2)
         call average_ranks(x(:, j), ranks, ties)
         sums = sums + ranks
         ties_total = ties_total + ties
      end do
      n = real(res%subjects, dp)
      r = real(res%raters, dp)
      denominator = r ** 2 * (n ** 3 - n)
      if (adjusted) denominator = denominator - r * ties_total
      if (denominator <= 0.0_dp) then
         res%status = irr_degenerate
         return
      end if
      res%value = 12.0_dp * sum((sums - r_mean(sums)) ** 2) / denominator
      res%statistic = r * (n - 1.0_dp) * res%value
      res%df1 = n - 1.0_dp
      res%p_value = r_pchisq(res%statistic, res%df1, lower_tail = .false.)
   end function

   pure function finn(ratings, s_levels, model) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater ordinal or numeric scores, with listwise NaN deletion.
      integer, intent(in) :: s_levels !! Number of possible scale values; must exceed one.
      character(*), intent(in), optional :: model !! 'oneway' (default) or 'twoway' variance model.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :)
      real(dp) :: msr, msw, msc, mse, expected, error
      character(:), allocatable :: model_used
      call complete_ratings(ratings, x, res, 2)
      if (res%status /= irr_ok) return
      model_used = 'oneway'
      if (present(model)) model_used = model
      if (s_levels < 2 .or. (model_used /= 'oneway' .and. model_used /= 'twoway')) then
         res%status = irr_invalid
         return
      end if
      call mean_squares(x, msr, msw, msc, mse)
      expected = (real(s_levels, dp) ** 2 - 1.0_dp) / 12.0_dp
      error = msw
      if (model_used == 'twoway') error = mse
      res%value = 1.0_dp - error / expected
      res%statistic = divide(expected, error)
      res%df1 = ieee_value(0.0_dp, ieee_positive_inf)
      res%df2 = real(res%subjects * (res%raters - 1), dp)
      ! Infinite numerator degrees of freedom: F = df2 / chi-square(df2).
      res%p_value = r_pchisq(res%df2 * error / expected, res%df2)
   end function

   pure function robinson(ratings) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Complete-case quantitative ratings with subjects in rows and raters in columns.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :)
      real(dp) :: msr, msw, msc, mse
      call complete_ratings(ratings, x, res, 2)
      if (res%status /= irr_ok) return
      call mean_squares(x, msr, msw, msc, mse)
      res%value = divide(msr, msr + real(res%raters - 1, dp) * mse)
      if (ieee_is_nan(res%value)) res%status = irr_degenerate
   end function

   pure function iota(ratings, scaledata, standardize) result(res)
      real(dp), intent(in) :: ratings(:, :, :) !! Subject-by-rater-by-variable data; rows with any NaN are removed across all
      !! variables.
      character(*), intent(in), optional :: scaledata !! 'quantitative' (default) or 'nominal' category distances.
      logical, intent(in), optional :: standardize !! Scale quantitative variables to unit sample standard deviation; default false.
      type(irr_result) :: res
      real(dp), allocatable :: flat(:, :), clean(:, :), x(:, :), values(:)
      integer, allocatable :: codes(:, :)
      real(dp) :: observed, expected, msr, msw, msc, mse, total, scale
      integer :: n, r, v, j, k, status
      logical :: standard
      character(:), allocatable :: method
      call init_result(res)
      r = size(ratings, 2)
      v = size(ratings, 3)
      if (r < 2 .or. v < 1) return
      if (size(ratings, kind = 8) > int(max_cells, 8) .or. r > max_raters) then
         res%status = irr_resource
         return
      end if
      method = 'quantitative'
      if (present(scaledata)) method = scaledata
      if (method /= 'quantitative' .and. method /= 'nominal') return
      standard = .false.
      if (present(standardize)) standard = standardize
      flat = reshape(ratings, [size(ratings, 1), r * v])
      call complete_ratings(flat, clean, res, 2)
      res%raters = r
      if (res%status /= irr_ok) return
      n = res%subjects
      observed = 0.0_dp
      expected = 0.0_dp
      do j = 1, v
         x = clean(:, (j - 1) * r + 1:j * r)
         if (method == 'quantitative') then
            if (standard) then
               scale = sqrt(r_variance(reshape(x, [n * r])))
               if (scale <= 0.0_dp) then
                  res%status = irr_degenerate
                  return
               end if
               x = (x - sum(x) / real(n * r, dp)) / scale
            end if
            call contribution(x, total, msw, msc)
            observed = observed + msw
            expected = expected + real(r - 1, dp) * total + msc
         else
            call encode_categories(x, codes, values, status)
            if (status /= irr_ok) then
               res%status = status
               return
            end if
            do k = 1, size(values)
               x = merge(1.0_dp / sqrt(2.0_dp), 0.0_dp, codes == k)
               call contribution(x, total, msw, msc)
               observed = observed + msw
               expected = expected + real(r - 1, dp) * total + msc
            end do
         end if
      end do
      res%value = 1.0_dp - divide(real(r, dp) * observed, expected)
      if (ieee_is_nan(res%value)) res%status = irr_degenerate
   end function

   pure subroutine contribution(x, total, within, columns)
      real(dp), intent(in) :: x(:, :) !! A finite subject-by-rater quantitative variable or category indicator.
      real(dp), intent(out) :: total !! Total sum of squares around the grand mean.
      real(dp), intent(out) :: within !! Within-subject sum of squares.
      real(dp), intent(out) :: columns !! Between-rater sum of squares.
      real(dp) :: msr, msw, msc, mse
      call mean_squares(x, msr, msw, msc, mse)
      total = sum((x - sum(x) / real(size(x), dp)) ** 2)
      within = msw * real(size(x, 1) * (size(x, 2) - 1), dp)
      columns = msc * real(size(x, 2) - 1, dp)
   end subroutine

   pure function kripp_alpha(ratings, method, legacy_complete_weights) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Rater-by-subject matrix; pairwise available ratings are used within each subject.
      character(*), intent(in), optional :: method !! 'nominal' (default), 'ordinal', 'interval', or 'ratio' distance.
      logical, intent(in), optional :: legacy_complete_weights !! Default true reproduces R 0.85 complete-data weighting; false
      !! always normalizes by available raters minus one.
      type(irr_result) :: res
      real(dp), allocatable :: values(:, :), levels(:), margins(:)
      integer, allocatable :: codes(:, :)
      integer :: i, j, k, a, b, n, nr, nc, status, count_present
      real(dp) :: denominator, observed, expected, distance, weight
      logical :: legacy, missing
      character(:), allocatable :: distance_method
      call init_result(res)
      res%raters = size(ratings, 1)
      res%subjects = size(ratings, 2)
      nr = res%raters
      n = res%subjects
      if (nr < 2 .or. n < 1) return
      if (nr > max_raters .or. size(ratings, kind = 8) > int(max_cells, 8)) then
         res%status = irr_resource
         return
      end if
      distance_method = 'nominal'
      if (present(method)) distance_method = method
      select case (distance_method)
      case ('nominal', 'ordinal', 'interval', 'ratio')
      case default
         return
      end select
      legacy = .true.
      if (present(legacy_complete_weights)) legacy = legacy_complete_weights
      missing = any(ieee_is_nan(ratings))
      allocate(values(count(.not. ieee_is_nan(ratings)), 1))
      values(:, 1) = pack(ratings, .not. ieee_is_nan(ratings))
      if (size(values, 1) == 0) then
         res%status = irr_degenerate
         return
      end if
      call encode_categories(values, codes, levels, status)
      if (status /= irr_ok) then
         res%status = status
         return
      end if
      nc = size(levels)
      res%levels = levels
      allocate(res%coincidence(nc, nc))
      res%coincidence = 0.0_dp
      do k = 1, n
         count_present = count(.not. ieee_is_nan(ratings(:, k)))
         if (count_present < 2) cycle
         denominator = real(count_present - 1, dp)
         if (legacy .and. .not. missing) denominator = 1.0_dp
         do i = 1, nr - 1
            if (ieee_is_nan(ratings(i, k))) cycle
            a = findloc(levels, ratings(i, k), dim = 1)
            do j = i + 1, nr
               if (ieee_is_nan(ratings(j, k))) cycle
               b = findloc(levels, ratings(j, k), dim = 1)
               weight = 1.0_dp / denominator
               res%coincidence(a, b) = res%coincidence(a, b) + weight
               res%coincidence(b, a) = res%coincidence(b, a) + weight
            end do
         end do
      end do
      res%status = irr_ok
      if (nc == 1) then
         res%value = 1.0_dp
         return
      end if
      margins = sum(res%coincidence, dim = 1)
      observed = 0.0_dp
      expected = 0.0_dp
      do b = 2, nc
         do a = 1, b - 1
            select case (distance_method)
            case ('nominal')
               distance = 1.0_dp
            case ('ordinal')
               distance = (sum(margins(a:b)) - (margins(a) + margins(b)) / 2.0_dp) ** 2
            case ('interval')
               distance = (levels(a) - levels(b)) ** 2
            case ('ratio')
               if (levels(a) + levels(b) == 0.0_dp) then
                  res%status = irr_degenerate
                  return
               end if
               distance = ((levels(a) - levels(b)) / (levels(a) + levels(b))) ** 2
            end select
            observed = observed + res%coincidence(a, b) * distance
            expected = expected + margins(a) * margins(b) * distance
         end do
      end do
      res%value = 1.0_dp - divide((sum(margins) - 1.0_dp) * observed, expected)
      if (ieee_is_nan(res%value)) res%status = irr_degenerate
   end function
end module
