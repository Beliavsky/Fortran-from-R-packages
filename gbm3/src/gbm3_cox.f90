! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_cox
   use gbm3_kinds, only : dp
   use gbm3_constants, only : GBM_TIES_BRESLOW
   use gbm3_math, only : argsort_desc_real
   use gbm3_types, only : gbm_options, gbm_tree
   implicit none
   private
   public :: cox_working_response, cox_deviance, cox_fit_best_constant, cox_bag_improvement
   public :: cox_baseline_hazard

contains

   subroutine cox_working_response(surv, strata, offset, f, weight, in_bag, options, residual)
      real(dp), intent(in) :: surv(:, :), offset(:), f(:), weight(:)
      integer, intent(in) :: strata(:)
      logical, intent(in) :: in_bag(:)
      type(gbm_options), intent(in) :: options
      real(dp), intent(out) :: residual(size(f))
      real(dp), allocatable :: mart(:), eta(:)
      real(dp) :: loglik
      allocate(mart(size(f)), eta(size(f)))
      eta = f + offset
      call cox_loglik_martingale(surv, strata, eta, weight, in_bag, options%cox_ties, loglik, mart)
      residual = 0.0_dp
      where (in_bag) residual = weight * mart
   end subroutine cox_working_response

   real(dp) function cox_deviance(surv, strata, offset, f, weight, options) result(dev)
      real(dp), intent(in) :: surv(:, :), offset(:), f(:), weight(:)
      integer, intent(in) :: strata(:)
      type(gbm_options), intent(in) :: options
      logical, allocatable :: include(:)
      real(dp), allocatable :: mart(:), eta(:)
      real(dp) :: loglik
      allocate(include(size(f)), mart(size(f)), eta(size(f)))
      include = .true.
      eta = f + offset
      call cox_loglik_martingale(surv, strata, eta, weight, include, options%cox_ties, loglik, mart)
      dev = -loglik
   end function cox_deviance

   subroutine cox_fit_best_constant(surv, strata, offset, f, weight, in_bag, assignment, min_obs, options, tree)
      real(dp), intent(in) :: surv(:, :), offset(:), f(:), weight(:)
      integer, intent(in) :: strata(:), assignment(:), min_obs
      logical, intent(in) :: in_bag(:)
      type(gbm_options), intent(in) :: options
      type(gbm_tree), intent(inout) :: tree
      real(dp), allocatable :: mart(:), eta(:), expected(:), actual(:)
      real(dp) :: loglik, status, prior
      integer :: i, node, status_col

      status_col = size(surv, 2)
      allocate(mart(size(f)), eta(size(f)), expected(tree%n_nodes), actual(tree%n_nodes))
      eta = f + offset
      call cox_loglik_martingale(surv, strata, eta, weight, in_bag, options%cox_ties, loglik, mart)
      prior = 1.0_dp / options%cox_prior_node_coeff_var
      expected = prior
      actual = prior
      do i = 1, size(f)
         if (.not. in_bag(i)) cycle
         node = assignment(i)
         if (tree%nodes(node)%num_obs < min_obs) cycle
         status = surv(i, status_col)
         expected(node) = expected(node) + max(0.0_dp, status - mart(i))
         actual(node) = actual(node) + status
      end do
      do node = 1, tree%n_nodes
         if (.not. tree%nodes(node)%is_terminal) cycle
         tree%nodes(node)%prediction = log(actual(node) / expected(node))
      end do
   end subroutine cox_fit_best_constant

   real(dp) function cox_bag_improvement(surv, strata, offset, f, weight, delta, in_bag, options) result(improvement)
      real(dp), intent(in) :: surv(:, :), offset(:), f(:), weight(:), delta(:)
      integer, intent(in) :: strata(:)
      logical, intent(in) :: in_bag(:)
      type(gbm_options), intent(in) :: options
      logical, allocatable :: include(:)
      real(dp), allocatable :: eta(:), eta_adj(:), mart(:)
      real(dp) :: ll0, ll1
      allocate(include(size(f)), eta(size(f)), eta_adj(size(f)), mart(size(f)))
      include = .not. in_bag
      eta = f + offset
      eta_adj = eta
      where (include) eta_adj = eta + options%shrinkage * delta
      call cox_loglik_martingale(surv, strata, eta, weight, include, options%cox_ties, ll0, mart)
      call cox_loglik_martingale(surv, strata, eta_adj, weight, include, options%cox_ties, ll1, mart)
      improvement = ll1 - ll0
   end function cox_bag_improvement

   subroutine cox_loglik_martingale(surv, strata, eta, weight, include, ties, loglik, mart)
      real(dp), intent(in) :: surv(:, :), eta(:), weight(:)
      integer, intent(in) :: strata(:), ties
      logical, intent(in) :: include(:)
      real(dp), intent(out) :: loglik, mart(size(eta))
      integer :: n, ncol, status_col, start_col, stop_col
      integer, allocatable :: strata_values(:), death_rows(:), ord(:)
      real(dp), allocatable :: death_times(:), score(:)
      integer :: ns, sidx, s, i, j, k, ndeath, nrisk
      real(dp) :: center, t, denom, dden, deathwt, avgwt, frac, hazard, ehazard, d
      logical :: risk, isdeath, same_time

      n = size(eta)
      ncol = size(surv, 2)
      if (size(surv, 1) /= n .or. size(weight) /= n .or. size(strata) /= n .or. size(include) /= n) &
         error stop "cox_loglik_martingale: shape mismatch"
      if (ncol == 2) then
         start_col = 0
         stop_col = 1
         status_col = 2
      else if (ncol == 3) then
         start_col = 1
         stop_col = 2
         status_col = 3
      else
         error stop "Cox response must have 2 (time,status) or 3 (start,stop,status) columns"
      end if

      mart = 0.0_dp
      do i = 1, n
         if (include(i) .and. surv(i, status_col) > 0.5_dp) mart(i) = 1.0_dp
      end do
      loglik = 0.0_dp
      call unique_int(strata, strata_values)
      ns = size(strata_values)

      do sidx = 1, ns
         s = strata_values(sidx)
         nrisk = count(include .and. strata == s)
         if (nrisk == 0) cycle
         center = -huge(1.0_dp)
         do i = 1, n
            if (include(i) .and. strata(i) == s) center = max(center, eta(i))
         end do
         allocate(score(n))
         score = 0.0_dp
         do i = 1, n
            if (include(i) .and. strata(i) == s) score(i) = exp(eta(i) - center)
         end do

         ndeath = count(include .and. strata == s .and. surv(:, status_col) > 0.5_dp)
         if (ndeath == 0) then
            deallocate(score)
            cycle
         end if
         allocate(death_rows(ndeath), death_times(ndeath), ord(ndeath))
         j = 0
         do i = 1, n
            if (include(i) .and. strata(i) == s .and. surv(i, status_col) > 0.5_dp) then
               j = j + 1
               death_rows(j) = i
               death_times(j) = surv(i, stop_col)
            end if
         end do
         call argsort_desc_real(death_times, ord)

         j = 1
         do while (j <= ndeath)
            t = death_times(ord(j))
            k = j + 1
            do while (k <= ndeath)
               same_time = .not. (death_times(ord(k)) < t .or. death_times(ord(k)) > t)
               if (.not. same_time) exit
               k = k + 1
            end do

            denom = 0.0_dp
            dden = 0.0_dp
            deathwt = 0.0_dp
            d = 0.0_dp
            do i = 1, n
               if (.not. include(i) .or. strata(i) /= s) cycle
               if (start_col == 0) then
                  risk = surv(i, stop_col) >= t
               else
                  risk = surv(i, start_col) < t .and. surv(i, stop_col) >= t
               end if
               if (risk) denom = denom + weight(i) * score(i)
               isdeath = surv(i, status_col) > 0.5_dp .and. &
                         .not. (surv(i, stop_col) < t .or. surv(i, stop_col) > t)
               if (isdeath) then
                  deathwt = deathwt + weight(i)
                  dden = dden + weight(i) * score(i)
                  d = d + weight(i) * (eta(i) - center)
               end if
            end do
            if (denom <= 0.0_dp) then
               j = k
               cycle
            end if
            loglik = loglik + d
            if (ties == GBM_TIES_BRESLOW .or. k - j == 1) then
               loglik = loglik - deathwt * log(denom)
               hazard = deathwt / denom
               ehazard = hazard
            else
               hazard = 0.0_dp
               ehazard = 0.0_dp
               avgwt = deathwt / real(k - j, dp)
               do i = 0, k - j - 1
                  frac = real(i, dp) / real(k - j, dp)
                  loglik = loglik - avgwt * log(denom - frac * dden)
                  hazard = hazard + avgwt / (denom - frac * dden)
                  ehazard = ehazard + (1.0_dp - frac) * avgwt / (denom - frac * dden)
               end do
            end if

            do i = 1, n
               if (.not. include(i) .or. strata(i) /= s) cycle
               if (start_col == 0) then
                  risk = surv(i, stop_col) >= t
               else
                  risk = surv(i, start_col) < t .and. surv(i, stop_col) >= t
               end if
               if (.not. risk) cycle
               isdeath = surv(i, status_col) > 0.5_dp .and. &
                         .not. (surv(i, stop_col) < t .or. surv(i, stop_col) > t)
               if (isdeath) then
                  mart(i) = mart(i) - ehazard * score(i)
               else
                  mart(i) = mart(i) - hazard * score(i)
               end if
            end do
            j = k
         end do

         deallocate(score, death_rows, death_times, ord)
      end do
   end subroutine cox_loglik_martingale

   subroutine unique_int(x, values)
      integer, intent(in) :: x(:)
      integer, allocatable, intent(out) :: values(:)
      integer, allocatable :: tmp(:)
      integer :: i, n
      allocate(tmp(size(x)))
      n = 0
      do i = 1, size(x)
         if (n == 0) then
            n = 1
            tmp(n) = x(i)
         else if (.not. any(tmp(1:n) == x(i))) then
            n = n + 1
            tmp(n) = x(i)
         end if
      end do
      allocate(values(n))
      values = tmp(1:n)
   end subroutine unique_int

   subroutine cox_baseline_hazard(surv, strata, eta, weight, ties, times, hazard, cumulative)
      real(dp), intent(in) :: surv(:, :), eta(:), weight(:)
      integer, intent(in) :: strata(:), ties
      real(dp), allocatable, intent(out) :: times(:), hazard(:), cumulative(:)
      ! This utility currently supports a single stratum, matching the common
      ! baseline-hazard use case. Stratified hazards can be obtained by calling
      ! it separately on each stratum slice.
      integer :: n, stop_col, status_col, i, j, k, nd, ndeath
      integer, allocatable :: ord(:)
      real(dp), allocatable :: dt(:), score(:)
      real(dp) :: center, t, denom, dden, deathwt, avgwt, frac, h
      logical :: risk, isdeath, same_time
      associate(dummy => strata)
      end associate
      n = size(eta)
      if (size(surv, 2) == 2) then
         stop_col = 1
         status_col = 2
      else
         stop_col = 2
         status_col = 3
      end if
      nd = count(surv(:, status_col) > 0.5_dp)
      if (nd == 0) then
         allocate(times(0), hazard(0), cumulative(0))
         return
      end if
      allocate(dt(nd), ord(nd), score(n))
      j = 0
      do i = 1, n
         if (surv(i, status_col) > 0.5_dp) then
            j = j + 1
            dt(j) = surv(i, stop_col)
         end if
      end do
      call argsort_desc_real(dt, ord)
      ndeath = 1
      do j = 2, nd
         if (dt(ord(j)) < dt(ord(j - 1)) .or. dt(ord(j)) > dt(ord(j - 1))) ndeath = ndeath + 1
      end do
      allocate(times(ndeath), hazard(ndeath), cumulative(ndeath))
      center = maxval(eta)
      score = exp(eta - center)
      j = 1
      k = 0
      do while (j <= nd)
         t = dt(ord(j))
         i = j + 1
         do while (i <= nd)
            same_time = .not. (dt(ord(i)) < t .or. dt(ord(i)) > t)
            if (.not. same_time) exit
            i = i + 1
         end do
         denom = 0.0_dp
         dden = 0.0_dp
         deathwt = 0.0_dp
         do n = 1, size(eta)
            risk = surv(n, stop_col) >= t
            if (size(surv, 2) == 3) risk = risk .and. surv(n, 1) < t
            if (risk) denom = denom + weight(n) * score(n)
            isdeath = surv(n, status_col) > 0.5_dp .and. &
                      .not. (surv(n, stop_col) < t .or. surv(n, stop_col) > t)
            if (isdeath) then
               deathwt = deathwt + weight(n)
               dden = dden + weight(n) * score(n)
            end if
         end do
         h = 0.0_dp
         if (ties == GBM_TIES_BRESLOW .or. i - j == 1) then
            h = deathwt / denom
         else
            avgwt = deathwt / real(i - j, dp)
            do n = 0, i - j - 1
               frac = real(n, dp) / real(i - j, dp)
               h = h + avgwt / (denom - frac * dden)
            end do
         end if
         k = k + 1
         times(k) = t
         hazard(k) = h * exp(-center)
         j = i
      end do
      ! Return ascending time, as R-facing baseline-hazard utilities typically do.
      call reverse_real(times)
      call reverse_real(hazard)
      cumulative(1) = hazard(1)
      do i = 2, size(hazard)
         cumulative(i) = cumulative(i - 1) + hazard(i)
      end do
   end subroutine cox_baseline_hazard

   subroutine reverse_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, n
      real(dp) :: tmp
      n = size(x)
      do i = 1, n / 2
         tmp = x(i)
         x(i) = x(n + 1 - i)
         x(n + 1 - i) = tmp
      end do
   end subroutine reverse_real

end module gbm3_cox
