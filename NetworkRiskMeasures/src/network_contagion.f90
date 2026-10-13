! SPDX-License-Identifier: GPL-3.0-only
! Derived from NetworkRiskMeasures 0.1.7, Carlos Cinelli and Thiago Cristiano Silva.
module network_contagion
   use network_risk_base, only: dp, nan_value, nr_ok, nr_shape, nr_domain, nr_iterations, nr_numeric, impact_matrix
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   type, public :: contagion_result
      real(dp), allocatable :: initial(:, :) !! n-by-scenario input shocks; identity for default all-vertex shocks.
      real(dp), allocatable :: final(:, :) !! n-by-scenario raw final stress, which may exceed one.
      real(dp), allocatable :: summary(:, :) !! 5-by-scenario original/additional stress, original/additional losses, defaults.
      real(dp), allocatable :: history(:, :, :) !! Optional n-by-max_it-by-scenario history; unused entries are NaN.
      integer, allocatable :: iterations(:) !! Number of completed propagation steps for each scenario.
      logical, allocatable :: converged(:) !! True only when the final stress increment norm is below tolerance.
   end type contagion_result
   public :: contagion, summarize_contagion
contains
   pure function propagation(stress, propagated, threshold, single_hit, threshold_method) result(flow)
      real(dp), intent(in) :: stress(:) !! Current raw stress per vertex, nonnegative for valid model inputs.
      real(dp), intent(in) :: propagated(:) !! Cumulative fraction already propagated by each vertex.
      real(dp), intent(in) :: threshold !! Trigger level for the threshold method, strictly positive and at most one.
      logical, intent(in) :: single_hit !! Stop further DebtRank propagation after a vertex first transmits any positive stress.
      logical, intent(in) :: threshold_method !! True selects threshold/default cascades; false selects DebtRank.
      real(dp) :: flow(size(stress)), clipped(size(stress)), begin_at, stop_at, full_at
      clipped = min(stress, 1.0_dp)
      begin_at = 0.0_dp
      full_at = 1.0_dp
      stop_at = 1.0_dp
      if (single_hit .and. .not. threshold_method) stop_at = 0.0_dp
      if (threshold_method) then
         begin_at = threshold
         full_at = threshold
      end if
      flow = clipped - propagated
      where (clipped >= full_at) flow = 1.0_dp - propagated
      where (clipped < begin_at .or. propagated > stop_at) flow = 0.0_dp
   end function propagation

   pure subroutine summarize_contagion(initial, final, weights, buffer, summary, status, threshold, cap)
      real(dp), intent(in) :: initial(:, :) !! n-by-scenario original nonnegative stress values.
      real(dp), intent(in) :: final(:, :) !! n-by-scenario nonnegative final raw stress values.
      real(dp), intent(in) :: weights(:) !! Finite nonnegative node weights, normalized to sum one by this routine.
      real(dp), intent(in) :: buffer(:) !! Finite positive capital buffer per vertex, in monetary units.
      real(dp), intent(out) :: summary(:, :) !! 5-by-scenario stress, additional stress, losses, additional losses, new defaults.
      integer, intent(out) :: status !! Zero on success; explicit shape, domain or numerical status otherwise.
      real(dp), intent(in), optional :: threshold !! Strictly positive default threshold, default one; comparison uses 1e-12 slack.
      real(dp), intent(in), optional :: cap !! Nonnegative stress cap used only for weighted stress summaries; default one.
      real(dp) :: w(size(weights)), level, ceiling, scale
      integer :: n, scenarios, s
      summary = nan_value()
      status = nr_shape
      n = size(initial, 1)
      scenarios = size(initial, 2)
      if (n < 1 .or. scenarios < 1) return
      if (any(shape(final) /= shape(initial))) return
      if (any(shape(summary) /= [5, scenarios])) return
      if (size(weights) /= n .or. size(buffer) /= n) return
      status = nr_domain
      if (any(.not. ieee_is_finite(initial)) .or. any(.not. ieee_is_finite(final))) return
      if (any(initial < 0.0_dp) .or. any(final < 0.0_dp)) return
      if (any(.not. ieee_is_finite(weights)) .or. any(.not. ieee_is_finite(buffer))) return
      if (any(weights < 0.0_dp) .or. any(buffer <= 0.0_dp)) return
      scale = maxval(weights)
      if (scale <= 0.0_dp) return
      level = 1.0_dp
      if (present(threshold)) level = threshold
      ceiling = 1.0_dp
      if (present(cap)) ceiling = cap
      if (.not. ieee_is_finite(level) .or. .not. ieee_is_finite(ceiling)) return
      if (level <= 0.0_dp .or. ceiling < 0.0_dp) return
      w = weights / scale
      w = w / sum(w)
      do s = 1, scenarios
         summary(1, s) = dot_product(min(ceiling, initial(:, s)), w)
         summary(2, s) = dot_product(min(ceiling, final(:, s)), w) - summary(1, s)
         summary(3, s) = dot_product(initial(:, s), buffer)
         summary(4, s) = dot_product(final(:, s) - initial(:, s), buffer)
         summary(5, s) = real(count(final(:, s) > level - 1.0e-12_dp .and. initial(:, s) < level), dp)
      end do
      status = nr_ok
      if (any(.not. ieee_is_finite(summary))) status = nr_numeric
   end subroutine summarize_contagion

   pure subroutine contagion(exposures, result, status, buffer, shocks, weights, method, exposure_type, &
      single_hit, alpha, keep_history, abs_tol, max_it)
      real(dp), intent(in) :: exposures(:, :) !! Square finite nonnegative bilateral network in the selected orientation.
      type(contagion_result), intent(out) :: result !! Scenario stresses, summaries, convergence diagnostics and optional history.
      integer, intent(out) :: status !! nr_ok if every scenario converges; nr_iterations retains valid unfinished results.
      real(dp), intent(in), optional :: buffer(:) !! Positive vertex capital buffers; default ones, also used for loss summaries.
      real(dp), intent(in), optional :: shocks(:, :) !! n-by-scenario nonnegative initial stress; default identity.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative finite node weights with positive total; default equal weights.
      character(*), intent(in), optional :: method !! debtrank (default) or threshold; custom R callback environments are omitted.
      character(*), intent(in), optional :: exposure_type !! Orientation: assets (default), liabilities, impact or vulnerability.
      logical, intent(in), optional :: single_hit !! Single-transmission DebtRank; default false, ignored by threshold.
      real(dp), intent(in), optional :: alpha !! Threshold cascade trigger in (0,1], default one; ignored by DebtRank.
      logical, intent(in), optional :: keep_history !! Allocate and return step history when true; default false.
      real(dp), intent(in), optional :: abs_tol !! Positive stress-increment Euclidean norm tolerance, default epsilon**0.2.
      integer, intent(in), optional :: max_it !! Positive propagation step limit, default min(1000,10*n).
      real(dp) :: v(size(exposures, 1), size(exposures, 2)), b(size(exposures, 1)), w(size(exposures, 1))
      real(dp) :: st(size(exposures, 1)), previous(size(exposures, 1)), flow(size(exposures, 1))
      real(dp) :: propagated(size(exposures, 1)), tol, trigger, change
      character(:), allocatable :: choice
      logical :: hit, keep, cascade
      integer :: n, scenarios, limit, i, s, summary_status
      n = size(exposures, 1)
      status = nr_shape
      if (n < 1 .or. size(exposures, 2) /= n) return
      b = 1.0_dp
      w = 1.0_dp
      if (present(buffer)) then
         if (size(buffer) /= n) return
         b = buffer
      end if
      if (present(weights)) then
         if (size(weights) /= n) return
         w = weights
      end if
      scenarios = n
      if (present(shocks)) then
         if (size(shocks, 1) /= n .or. size(shocks, 2) < 1) return
         scenarios = size(shocks, 2)
      end if
      status = nr_domain
      if (any(.not. ieee_is_finite(w)) .or. any(w < 0.0_dp)) return
      if (maxval(w) <= 0.0_dp) return
      if (present(shocks)) then
         if (any(.not. ieee_is_finite(shocks)) .or. any(shocks < 0.0_dp)) return
      end if
      choice = 'debtrank'
      if (present(method)) choice = method
      if (choice /= 'debtrank' .and. choice /= 'threshold') return
      cascade = choice == 'threshold'
      hit = .false.
      if (present(single_hit)) hit = single_hit
      keep = .false.
      if (present(keep_history)) keep = keep_history
      trigger = 1.0_dp
      if (present(alpha)) trigger = alpha
      tol = epsilon(1.0_dp) ** 0.2_dp
      if (present(abs_tol)) tol = abs_tol
      if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp) return
      if (.not. ieee_is_finite(trigger)) return
      if (trigger <= 0.0_dp .or. trigger > 1.0_dp) return
      limit = min(1000, 10 * n)
      if (present(max_it)) limit = max_it
      if (limit < 1) return
      call impact_matrix(exposures, v, status, b, exposure_type=exposure_type)
      if (status /= nr_ok) return
      allocate(result%initial(n, scenarios), result%final(n, scenarios), result%summary(5, scenarios))
      allocate(result%iterations(scenarios), result%converged(scenarios))
      result%initial = 0.0_dp
      result%final = nan_value()
      result%summary = nan_value()
      result%iterations = 0
      result%converged = .false.
      if (present(shocks)) then
         result%initial = shocks
      else
         do i = 1, n
            result%initial(i, i) = 1.0_dp
         end do
      end if
      if (keep) then
         allocate(result%history(n, limit, scenarios))
         result%history = nan_value()
      end if
      do s = 1, scenarios
         st = result%initial(:, s)
         propagated = 0.0_dp
         flow = propagation(st, propagated, trigger, hit, cascade)
         previous = st
         do i = 1, limit
            propagated = propagated + flow
            st = previous + matmul(flow, v)
            result%iterations(s) = i
            if (keep) result%history(:, i, s) = st
            if (any(.not. ieee_is_finite(st))) then
               result%final(:, s) = st
               status = nr_numeric
               return
            end if
            change = norm2(st - previous)
            flow = propagation(st, propagated, trigger, hit, cascade)
            if (change < tol) then
               result%converged(s) = .true.
               exit
            end if
            previous = st
         end do
         result%final(:, s) = st
      end do
      call summarize_contagion(result%initial, result%final, w, b, result%summary, summary_status)
      status = summary_status
      if (status == nr_ok .and. .not. all(result%converged)) status = nr_iterations
   end subroutine contagion
end module network_contagion
