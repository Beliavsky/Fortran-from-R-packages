! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Implements the computational preprocessing behind respect.unordered.factors="order".
module ranger_factor_order
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   implicit none
   private

   public :: order_factors_regression, order_factors_classification, order_factors_survival
   public :: apply_factor_map, identity_factor_map

contains

   subroutine identity_factor_map(ncat, mapping)
      integer, intent(in) :: ncat(:)
      integer, allocatable, intent(out) :: mapping(:,:)
      integer :: j, k, maxcat
      maxcat = max(1, maxval(ncat))
      allocate(mapping(maxcat, size(ncat)))
      mapping = 0
      do j = 1, size(ncat)
         do k = 1, ncat(j)
            mapping(k, j) = k
         end do
      end do
   end subroutine identity_factor_map

   subroutine apply_factor_map(x, ncat, mapping, transformed)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:), mapping(:,:)
      real(dp), allocatable, intent(out) :: transformed(:,:)
      integer :: i, j, code

      allocate(transformed(size(x, 1), size(x, 2)))
      transformed = x
      do j = 1, size(x, 2)
         if (ncat(j) <= 1) cycle
         do i = 1, size(x, 1)
            if (ieee_is_nan(x(i, j))) cycle
            code = nint(x(i, j))
            if (code >= 1 .and. code <= ncat(j)) transformed(i, j) = real(mapping(code, j), dp)
         end do
      end do
   end subroutine apply_factor_map

   subroutine order_factors_regression(x, y, ncat, mapping, transformed)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, intent(in) :: ncat(:)
      integer, allocatable, intent(out) :: mapping(:,:)
      real(dp), allocatable, intent(out) :: transformed(:,:)
      real(dp), allocatable :: score(:)
      integer, allocatable :: order(:)
      integer :: j, k, i, nobs

      call identity_factor_map(ncat, mapping)
      do j = 1, size(ncat)
         if (ncat(j) <= 1) cycle
         allocate(score(ncat(j)), order(ncat(j)))
         do k = 1, ncat(j)
            score(k) = 0.0_dp
            nobs = 0
            do i = 1, size(y)
               if (ieee_is_nan(x(i, j))) cycle
               if (nint(x(i, j)) == k) then
                  score(k) = score(k) + y(i)
                  nobs = nobs + 1
               end if
            end do
            if (nobs > 0) then
               score(k) = score(k) / real(nobs, dp)
            else
               score(k) = huge(1.0_dp)
            end if
         end do
         call order_by_score(score, order)
         do k = 1, ncat(j)
            mapping(order(k), j) = k
         end do
         deallocate(score, order)
      end do
      call apply_factor_map(x, ncat, mapping, transformed)
   end subroutine order_factors_regression

   subroutine order_factors_classification(x, y, nclass, ncat, mapping, transformed)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:), nclass, ncat(:)
      integer, allocatable, intent(out) :: mapping(:,:)
      real(dp), allocatable, intent(out) :: transformed(:,:)
      real(dp), allocatable :: score(:)
      integer, allocatable :: order(:)
      integer :: j, k, i, nobs

      call identity_factor_map(ncat, mapping)
      do j = 1, size(ncat)
         if (ncat(j) <= 1) cycle
         allocate(score(ncat(j)), order(ncat(j)))
         if (nclass <= 2) then
            do k = 1, ncat(j)
               score(k) = 0.0_dp
               nobs = 0
               do i = 1, size(y)
                  if (ieee_is_nan(x(i, j))) cycle
                  if (nint(x(i, j)) == k) then
                     score(k) = score(k) + real(y(i), dp)
                     nobs = nobs + 1
                  end if
               end do
               if (nobs > 0) then
                  score(k) = score(k) / real(nobs, dp)
               else
                  score(k) = huge(1.0_dp)
               end if
            end do
         else
            call multiclass_pca_scores(x(:, j), y, nclass, ncat(j), score)
         end if
         call order_by_score(score, order)
         do k = 1, ncat(j)
            mapping(order(k), j) = k
         end do
         deallocate(score, order)
      end do
      call apply_factor_map(x, ncat, mapping, transformed)
   end subroutine order_factors_classification

   subroutine multiclass_pca_scores(x, y, nclass, ncat, score)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: y(:), nclass, ncat
      real(dp), intent(out) :: score(ncat)
      real(dp), allocatable :: counts(:,:), prob(:,:), weights(:), meanp(:), cov(:,:), centered(:,:), gram(:,:), pc(:)
      real(dp) :: total, denom, normv
      integer :: i, k, c, iter

      allocate(counts(ncat, nclass), prob(ncat, nclass), weights(ncat), meanp(nclass))
      allocate(cov(nclass, nclass), centered(nclass, nclass), gram(nclass, nclass), pc(nclass))
      counts = 0.0_dp
      do i = 1, size(y)
         if (ieee_is_nan(x(i))) cycle
         k = nint(x(i))
         if (k >= 1 .and. k <= ncat .and. y(i) >= 1 .and. y(i) <= nclass) &
            counts(k, y(i)) = counts(k, y(i)) + 1.0_dp
      end do
      weights = sum(counts, dim=2)
      total = sum(weights)
      prob = 0.0_dp
      do k = 1, ncat
         if (weights(k) > 0.0_dp) prob(k, :) = counts(k, :) / weights(k)
      end do
      if (total <= 0.0_dp .or. count(weights > 0.0_dp) < 2) then
         score = [(real(k, dp), k = 1, ncat)]
         return
      end if
      weights = weights / total
      meanp = matmul(weights, prob)
      cov = 0.0_dp
      do k = 1, ncat
         if (weights(k) <= 0.0_dp) cycle
         do i = 1, nclass
            do c = 1, nclass
               cov(i, c) = cov(i, c) + weights(k) * (prob(k, i) - meanp(i)) * (prob(k, c) - meanp(c))
            end do
         end do
      end do
      denom = 1.0_dp - sum(weights ** 2)
      if (denom > 0.0_dp) cov = cov / denom

      ! pca.order() applies prcomp() to the covariance matrix S itself.
      ! prcomp centers the columns before finding the first right singular vector.
      centered = cov
      do c = 1, nclass
         centered(:, c) = centered(:, c) - sum(centered(:, c)) / real(nclass, dp)
      end do
      gram = matmul(transpose(centered), centered)
      pc = [(real(c, dp), c = 1, nclass)]
      normv = sqrt(sum(pc ** 2))
      if (normv > 0.0_dp) pc = pc / normv
      do iter = 1, 200
         meanp = matmul(gram, pc)
         normv = sqrt(sum(meanp ** 2))
         if (normv <= tiny(1.0_dp)) exit
         meanp = meanp / normv
         if (sqrt(sum((meanp - pc) ** 2)) <= 1.0e-12_dp .or. &
            sqrt(sum((meanp + pc) ** 2)) <= 1.0e-12_dp) then
            pc = meanp
            exit
         end if
         pc = meanp
      end do
      score = matmul(prob, pc)
      do k = 1, ncat
         if (sum(counts(k, :)) <= 0.0_dp) score(k) = huge(1.0_dp)
      end do
   end subroutine multiclass_pca_scores

   subroutine order_factors_survival(x, time, event, ncat, mapping, transformed)
      real(dp), intent(in) :: x(:,:), time(:)
      integer, intent(in) :: event(:), ncat(:)
      integer, allocatable, intent(out) :: mapping(:,:)
      real(dp), allocatable, intent(out) :: transformed(:,:)
      real(dp), allocatable :: minimum_survival(:), quantile_time(:), score(:)
      integer, allocatable :: order(:)
      real(dp) :: common_probability
      integer :: j, k

      call identity_factor_map(ncat, mapping)
      do j = 1, size(ncat)
         if (ncat(j) <= 1) cycle
         allocate(minimum_survival(ncat(j)), quantile_time(ncat(j)), score(ncat(j)), order(ncat(j)))
         do k = 1, ncat(j)
            call km_minimum_survival(x(:, j), k, time, event, minimum_survival(k))
         end do
         if (any(minimum_survival >= 0.0_dp)) then
            common_probability = min(0.5_dp, 1.0_dp - maxval(minimum_survival, mask=minimum_survival >= 0.0_dp))
            common_probability = max(0.0_dp, common_probability)
         else
            common_probability = 0.0_dp
         end if
         do k = 1, ncat(j)
            if (minimum_survival(k) < 0.0_dp) then
               quantile_time(k) = -huge(1.0_dp)
            else
               call km_quantile(x(:, j), k, time, event, common_probability, quantile_time(k))
            end if
         end do
         score = quantile_time
         call order_by_score(score, order)
         do k = 1, ncat(j)
            mapping(order(k), j) = k
         end do
         deallocate(minimum_survival, quantile_time, score, order)
      end do
      call apply_factor_map(x, ncat, mapping, transformed)
   end subroutine order_factors_survival

   subroutine km_minimum_survival(x, category, time, event, minimum_survival)
      real(dp), intent(in) :: x(:), time(:)
      integer, intent(in) :: category, event(:)
      real(dp), intent(out) :: minimum_survival
      real(dp), allocatable :: event_times(:)
      real(dp) :: survival, at_risk, deaths
      integer :: i, j, ntime, nobs

      call category_event_times(x, category, time, event, event_times, ntime)
      nobs = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            if (nint(x(i)) == category) nobs = nobs + 1
         end if
      end do
      ! summary.survfit() defaults to censored = FALSE.  A stratum with no
      ! events therefore contributes no row to largest.quantile() and is
      ! prepended later as a missing level by ranger.R.
      if (nobs == 0 .or. ntime == 0) then
         minimum_survival = -1.0_dp
         return
      end if
      survival = 1.0_dp
      do j = 1, ntime
         at_risk = 0.0_dp
         deaths = 0.0_dp
         do i = 1, size(time)
            if (ieee_is_nan(x(i)) .or. nint(x(i)) /= category) cycle
            if (time(i) >= event_times(j)) at_risk = at_risk + 1.0_dp
            if (event(i) == 1 .and. exact_equal(time(i), event_times(j))) deaths = deaths + 1.0_dp
         end do
         if (at_risk > 0.0_dp) survival = survival * (1.0_dp - deaths / at_risk)
      end do
      minimum_survival = survival
   end subroutine km_minimum_survival

   subroutine km_quantile(x, category, time, event, probability, quantile_time)
      real(dp), intent(in) :: x(:), time(:), probability
      integer, intent(in) :: category, event(:)
      real(dp), intent(out) :: quantile_time
      real(dp), allocatable :: event_times(:)
      real(dp) :: survival, previous_cdf, cdf, at_risk, deaths
      real(dp) :: previous_time, next_time, last_time, tolerance
      integer :: i, j, ntime

      call category_event_times(x, category, time, event, event_times, ntime)
      if (ntime == 0) then
         quantile_time = huge(1.0_dp)
         return
      end if
      if (probability <= 0.0_dp) then
         quantile_time = 0.0_dp
         return
      end if

      last_time = 0.0_dp
      do i = 1, size(time)
         if (ieee_is_nan(x(i)) .or. nint(x(i)) /= category) cycle
         last_time = max(last_time, time(i))
      end do
      tolerance = sqrt(epsilon(1.0_dp))
      survival = 1.0_dp
      previous_cdf = 0.0_dp
      previous_time = 0.0_dp
      quantile_time = huge(1.0_dp)

      do j = 1, ntime
         at_risk = 0.0_dp
         deaths = 0.0_dp
         do i = 1, size(time)
            if (ieee_is_nan(x(i)) .or. nint(x(i)) /= category) cycle
            if (time(i) >= event_times(j)) at_risk = at_risk + 1.0_dp
            if (event(i) == 1 .and. exact_equal(time(i), event_times(j))) deaths = deaths + 1.0_dp
         end do
         if (at_risk > 0.0_dp) survival = survival * (1.0_dp - deaths / at_risk)
         cdf = 1.0_dp - survival

         if (probability <= cdf + tolerance) then
            if (abs(probability - previous_cdf) < tolerance) then
               quantile_time = 0.5_dp * (previous_time + event_times(j))
            else if (abs(probability - cdf) < tolerance) then
               if (j < ntime) then
                  next_time = event_times(j + 1)
               else
                  next_time = last_time
               end if
               quantile_time = 0.5_dp * (event_times(j) + next_time)
            else
               quantile_time = event_times(j)
            end if
            return
         end if
         previous_cdf = cdf
         previous_time = event_times(j)
      end do
   end subroutine km_quantile

   subroutine category_event_times(x, category, time, event, values, nvalue)
      real(dp), intent(in) :: x(:), time(:)
      integer, intent(in) :: category, event(:)
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out) :: nvalue
      real(dp), allocatable :: work(:)
      real(dp) :: key
      integer :: i, j

      allocate(work(size(time)))
      nvalue = 0
      do i = 1, size(time)
         if (ieee_is_nan(x(i))) cycle
         if (nint(x(i)) == category .and. event(i) == 1) then
            nvalue = nvalue + 1
            work(nvalue) = time(i)
         end if
      end do
      do i = 2, nvalue
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      j = 0
      do i = 1, nvalue
         if (j == 0) then
            j = 1
            work(j) = work(i)
         else if (.not. exact_equal(work(i), work(j))) then
            j = j + 1
            work(j) = work(i)
         end if
      end do
      nvalue = j
      allocate(values(nvalue))
      if (nvalue > 0) values = work(1:nvalue)
   end subroutine category_event_times

   subroutine order_by_score(score, order)
      real(dp), intent(in) :: score(:)
      integer, intent(out) :: order(size(score))
      integer :: i, j, key
      order = [(i, i = 1, size(score))]
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (score(order(j)) <= score(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine order_by_score

   pure logical function exact_equal(a, b) result(equal)
      real(dp), intent(in) :: a, b
      equal = abs(a - b) <= 0.0_dp
   end function exact_equal

end module ranger_factor_order
