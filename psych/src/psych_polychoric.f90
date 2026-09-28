module psych_polychoric
   use r_kinds, only : dp, r_pi
   use psych_distributions, only : normal_cdf, normal_quantile
   implicit none
   private
   public :: bivariate_normal_cdf, tetrachoric_table, polychoric_categories, phi_to_tetrachoric
   public :: integer_contingency_table
contains

   pure subroutine integer_contingency_table(x, y, table)
      integer, intent(in) :: x(:) !! Integer-coded first variable; categories may start at zero, one, or a negative value.
      integer, intent(in) :: y(:) !! Integer-coded second variable with the same number of observations as `x`.
      integer, allocatable, intent(out) :: table(:, :) !! Two-way frequency table preserving unused lower category labels.
      integer :: base_x, base_y, i, ix, iy, max_x, max_y, min_x, min_y, n

      n = min(size(x), size(y))
      if (n == 0) then
         allocate(table(0, 0))
         return
      end if

      min_x = minval(x(1:n))
      max_x = maxval(x(1:n))
      min_y = minval(y(1:n))
      max_y = maxval(y(1:n))
      base_x = min(1, min_x)
      base_y = min(1, min_y)
      allocate(table(max_x - base_x + 1, max_y - base_y + 1))
      table = 0

      do i = 1, n
         ix = x(i) - base_x + 1
         iy = y(i) - base_y + 1
         table(ix, iy) = table(ix, iy) + 1
      end do
   end subroutine integer_contingency_table

   pure real(dp) function bivariate_normal_cdf(h, k, rho) result(value)
      real(dp), intent(in) :: h !! First standardized upper integration limit.
      real(dp), intent(in) :: k !! Second standardized upper integration limit.
      real(dp), intent(in) :: rho !! Correlation constrained to the open interval (-1,1).
      integer, parameter :: panels = 160
      real(dp) :: r, step, t, integrand, sum_value
      integer :: i
      r = min(max(rho, -0.999999_dp), 0.999999_dp)
      if (abs(r) <= epsilon(1.0_dp)) then
         value = normal_cdf(h) * normal_cdf(k)
         return
      end if
      step = r / real(panels, dp)
      sum_value = 0.0_dp
      do i = 0, panels
         t = real(i, dp) * step
         integrand = bvn_rho_integrand(t, h, k)
         if (i == 0 .or. i == panels) then
            sum_value = sum_value + integrand
         else if (mod(i, 2) == 0) then
            sum_value = sum_value + 2.0_dp * integrand
         else
            sum_value = sum_value + 4.0_dp * integrand
         end if
      end do
      value = normal_cdf(h) * normal_cdf(k) + step * sum_value / 3.0_dp
      value = min(max(value, 0.0_dp), 1.0_dp)
   end function bivariate_normal_cdf

   pure subroutine tetrachoric_table(table, rho, correction)
      real(dp), intent(in) :: table(2, 2) !! Nonnegative 2-by-2 frequency or probability table.
      real(dp), intent(out) :: rho !! Estimated latent bivariate-normal correlation.
      real(dp), intent(in), optional :: correction !! Continuity addition applied to zero cells; defaults to 0.5 for count tables.
      real(dp) :: work(2, 2), total, add, row_p, col_p, target, h, k, lo, hi, mid, predicted
      integer :: iter
      work = max(table, 0.0_dp)
      add = 0.5_dp
      if (present(correction)) add = correction
      if (add > 0.0_dp .and. maxval(work) > 1.0_dp) then
         where (work <= tiny(1.0_dp)) work = work + add
      end if
      total = sum(work)
      if (total <= tiny(1.0_dp)) then
         rho = 0.0_dp
         return
      end if
      work = work / total
      row_p = work(1, 1) + work(1, 2)
      col_p = work(1, 1) + work(2, 1)
      target = work(1, 1)
      row_p = min(max(row_p, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
      col_p = min(max(col_p, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
      h = normal_quantile(row_p)
      k = normal_quantile(col_p)
      lo = -0.999_dp
      hi = 0.999_dp
      do iter = 1, 80
         mid = 0.5_dp * (lo + hi)
         predicted = bivariate_normal_cdf(h, k, mid)
         if (predicted < target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      rho = 0.5_dp * (lo + hi)
   end subroutine tetrachoric_table

   pure subroutine polychoric_categories(x, y, rho, correction)
      integer, intent(in) :: x(:) !! Ordered integer categories for the first variable.
      integer, intent(in) :: y(:) !! Ordered integer categories for the second variable.
      real(dp), intent(out) :: rho !! Maximum-likelihood latent normal correlation estimate.
      real(dp), intent(in), optional :: correction !! Optional pseudocount added to empty cells before normalization.
      integer, allocatable :: counts(:, :)
      real(dp), allocatable :: table(:, :), row_cut(:), col_cut(:)
      real(dp) :: add, a, b, c, d, fc, fd, ratio
      integer :: i, nx, ny, xmin, ymin, iter
      if (size(x) /= size(y) .or. size(x) == 0) then
         rho = 0.0_dp
         return
      end if
      xmin = minval(x)
      ymin = minval(y)
      nx = maxval(x) - xmin + 1
      ny = maxval(y) - ymin + 1
      allocate(counts(nx, ny))
      counts = 0
      do i = 1, size(x)
         counts(x(i) - xmin + 1, y(i) - ymin + 1) = counts(x(i) - xmin + 1, y(i) - ymin + 1) + 1
      end do
      allocate(table(nx, ny))
      table = real(counts, dp)
      add = 0.5_dp
      if (present(correction)) add = correction
      if (add > 0.0_dp) then
         where (counts == 0) table = add
      end if
      table = table / sum(table)
      call category_cuts(sum(table, dim=2), row_cut)
      call category_cuts(sum(table, dim=1), col_cut)
      ratio = 0.6180339887498948482_dp
      a = -0.995_dp
      b = 0.995_dp
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = poly_objective(c, table, row_cut, col_cut)
      fd = poly_objective(d, table, row_cut, col_cut)
      do iter = 1, 80
         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = poly_objective(c, table, row_cut, col_cut)
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = poly_objective(d, table, row_cut, col_cut)
         end if
      end do
      rho = 0.5_dp * (a + b)
   end subroutine polychoric_categories

   pure subroutine phi_to_tetrachoric(phi, predictor_rate, criterion_rate, rho, sample_size, correction)
      real(dp), intent(in) :: phi !! Observed phi correlation for a binary 2-by-2 relation.
      real(dp), intent(in) :: predictor_rate !! Marginal success proportion or count for the predictor.
      real(dp), intent(in) :: criterion_rate !! Marginal success proportion or count for the criterion.
      real(dp), intent(out) :: rho !! Tetrachoric correlation implied by phi and the supplied margins.
      real(dp), intent(in), optional :: sample_size !! Optional count denominator when margins are supplied as counts.
      real(dp), intent(in), optional :: correction !! Optional continuity correction passed to tetrachoric estimation.
      real(dp) :: rp, cp, n, cross, table(2, 2)
      n = 1.0_dp
      if (present(sample_size)) n = sample_size
      rp = predictor_rate / n
      cp = criterion_rate / n
      rp = min(max(rp, 0.0_dp), 1.0_dp)
      cp = min(max(cp, 0.0_dp), 1.0_dp)
      cross = phi * sqrt(max(rp * (1.0_dp - rp) * cp * (1.0_dp - cp), 0.0_dp))
      table(1, 1) = rp * cp + cross
      table(1, 2) = rp * (1.0_dp - cp) - cross
      table(2, 1) = (1.0_dp - rp) * cp - cross
      table(2, 2) = (1.0_dp - rp) * (1.0_dp - cp) + cross
      table = max(table, 0.0_dp)
      if (present(sample_size)) table = table * n
      call tetrachoric_table(table, rho, correction)
   end subroutine phi_to_tetrachoric

   pure real(dp) function bvn_rho_integrand(rho, h, k) result(value)
      real(dp), intent(in) :: rho !! Intermediate correlation in Plackett's integral representation.
      real(dp), intent(in) :: h !! First standardized upper limit.
      real(dp), intent(in) :: k !! Second standardized upper limit.
      real(dp) :: denominator
      denominator = max(1.0_dp - rho * rho, tiny(1.0_dp))
      value = exp(-(h * h - 2.0_dp * rho * h * k + k * k) / (2.0_dp * denominator)) / &
         (2.0_dp * r_pi * sqrt(denominator))
   end function bvn_rho_integrand

   pure subroutine category_cuts(probabilities, cuts)
      real(dp), intent(in) :: probabilities(:) !! Ordered category probabilities summing to one.
      real(dp), allocatable, intent(out) :: cuts(:) !! Interior standard-normal thresholds between categories.
      real(dp) :: cumulative
      integer :: i
      allocate(cuts(max(size(probabilities) - 1, 0)))
      cumulative = 0.0_dp
      do i = 1, size(cuts)
         cumulative = cumulative + probabilities(i)
         cuts(i) = normal_quantile(min(max(cumulative, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp))
      end do
   end subroutine category_cuts

   pure real(dp) function poly_objective(rho, table, row_cut, col_cut) result(value)
      real(dp), intent(in) :: rho !! Candidate latent correlation.
      real(dp), intent(in) :: table(:, :) !! Normalized contingency probabilities.
      real(dp), intent(in) :: row_cut(:) !! Interior thresholds for row categories.
      real(dp), intent(in) :: col_cut(:) !! Interior thresholds for column categories.
      real(dp) :: lower_r, upper_r, lower_c, upper_c, probability
      integer :: i, j
      value = 0.0_dp
      do i = 1, size(table, 1)
         lower_r = category_lower(row_cut, i)
         upper_r = category_upper(row_cut, i)
         do j = 1, size(table, 2)
            lower_c = category_lower(col_cut, j)
            upper_c = category_upper(col_cut, j)
            probability = rectangle_probability(lower_r, upper_r, lower_c, upper_c, rho)
            probability = max(probability, tiny(1.0_dp))
            value = value - table(i, j) * log(probability)
         end do
      end do
   end function poly_objective

   pure real(dp) function rectangle_probability(lr, ur, lc, uc, rho) result(value)
      real(dp), intent(in) :: lr !! Lower first-dimension threshold, with +/-9 used as numerical infinity.
      real(dp), intent(in) :: ur !! Upper first-dimension threshold.
      real(dp), intent(in) :: lc !! Lower second-dimension threshold.
      real(dp), intent(in) :: uc !! Upper second-dimension threshold.
      real(dp), intent(in) :: rho !! Latent normal correlation.
      value = bivariate_normal_cdf(ur, uc, rho) - bivariate_normal_cdf(lr, uc, rho) - &
         bivariate_normal_cdf(ur, lc, rho) + bivariate_normal_cdf(lr, lc, rho)
      value = min(max(value, 0.0_dp), 1.0_dp)
   end function rectangle_probability

   pure real(dp) function category_lower(cuts, category) result(value)
      real(dp), intent(in) :: cuts(:) !! Interior thresholds for one ordered categorical variable.
      integer, intent(in) :: category !! One-based category index whose lower threshold is requested.
      if (category <= 1) then
         value = -9.0_dp
      else
         value = cuts(category - 1)
      end if
   end function category_lower

   pure real(dp) function category_upper(cuts, category) result(value)
      real(dp), intent(in) :: cuts(:) !! Interior thresholds for one ordered categorical variable.
      integer, intent(in) :: category !! One-based category index whose upper threshold is requested.
      if (category > size(cuts)) then
         value = 9.0_dp
      else
         value = cuts(category)
      end if
   end function category_upper
end module psych_polychoric
