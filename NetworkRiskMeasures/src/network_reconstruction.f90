! SPDX-License-Identifier: GPL-3.0-only
! Derived from NetworkRiskMeasures 0.1.7, Carlos Cinelli and Thiago Cristiano Silva.
module network_reconstruction
   use iso_fortran_env, only: int64
   use network_risk_base, only: dp, nan_value, nr_ok, nr_shape, nr_domain, nr_iterations, &
      nr_infeasible, nr_random_exhausted, nr_numeric
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: max_ent, min_dens, matrix_estimation
contains
   pure subroutine check_margins(rowsums, colsums, n, total, status, abs_tol)
      real(dp), intent(in) :: rowsums(:) !! Target nonnegative finite row sums of a zero-diagonal square matrix.
      real(dp), intent(in) :: colsums(:) !! Target nonnegative finite column sums in the same vertex order.
      integer, intent(out) :: n !! Common number of vertices, including zero-margin vertices.
      real(dp), intent(out) :: total !! Sum of target row sums when data are valid.
      integer, intent(out) :: status !! Zero if dimensions, values, balance and zero-diagonal feasibility hold.
      real(dp), intent(in) :: abs_tol !! Positive absolute L1 stopping tolerance used for feasibility slack.
      real(dp) :: slack
      n = size(rowsums)
      total = 0.0_dp
      status = nr_shape
      if (n < 1 .or. size(colsums) /= n) return
      status = nr_domain
      if (.not. ieee_is_finite(abs_tol) .or. abs_tol <= 0.0_dp) return
      if (any(.not. ieee_is_finite(rowsums)) .or. any(.not. ieee_is_finite(colsums))) return
      if (any(rowsums < 0.0_dp) .or. any(colsums < 0.0_dp)) return
      total = sum(rowsums)
      if (.not. ieee_is_finite(total) .or. .not. ieee_is_finite(sum(colsums))) return
      slack = max(abs_tol * 0.25_dp, 32.0_dp * epsilon(1.0_dp) * total)
      status = nr_infeasible
      if (abs(total - sum(colsums)) > slack) return
      if (any(rowsums > total - colsums + slack)) return
      status = nr_ok
   end subroutine check_margins

   pure subroutine max_ent(rowsums, colsums, result, status, max_it, abs_tol, iterations, residual)
      real(dp), intent(in) :: rowsums(:) !! Nonnegative finite target row totals, in vertex order.
      real(dp), intent(in) :: colsums(:) !! Nonnegative finite column totals with the same length and total.
      real(dp), intent(out) :: result(:, :) !! Zero-diagonal maximum-entropy estimate, including last iterate on nonconvergence.
      integer, intent(out) :: status !! nr_ok, nr_iterations, nr_infeasible, nr_shape or nr_domain.
      integer, intent(in), optional :: max_it !! Maximum row/column scaling sweeps, default 100000 and at least one.
      real(dp), intent(in), optional :: abs_tol !! Positive L1 marginal-error tolerance, default 0.001.
      integer, intent(out), optional :: iterations !! Number of completed scaling sweeps; zero for invalid or zero-total inputs.
      real(dp), intent(out), optional :: residual !! L1 sum of row and column discrepancies, NaN for invalid inputs.
      integer :: n, limit, t, i, j
      real(dp) :: total, tol, current, error
      result = nan_value()
      if (present(iterations)) iterations = 0
      if (present(residual)) residual = nan_value()
      tol = 1.0e-3_dp
      if (present(abs_tol)) tol = abs_tol
      limit = 100000
      if (present(max_it)) limit = max_it
      status = nr_domain
      if (limit < 1) return
      call check_margins(rowsums, colsums, n, total, status, tol)
      if (status /= nr_ok) return
      status = nr_shape
      if (any(shape(result) /= [n, n])) return
      result = 0.0_dp
      status = nr_ok
      if (total == 0.0_dp) then
         if (present(residual)) residual = sum(colsums)
         return
      end if
      result = 1.0_dp
      do i = 1, n
         result(i, i) = 0.0_dp
      end do
      do t = 1, limit
         do i = 1, n
            current = sum(result(i, :))
            if (current > 0.0_dp) result(i, :) = (result(i, :) / current) * rowsums(i)
         end do
         do j = 1, n
            current = sum(result(:, j))
            if (current > 0.0_dp) result(:, j) = (result(:, j) / current) * colsums(j)
         end do
         error = sum(abs(rowsums - sum(result, dim=2))) + sum(abs(colsums - sum(result, dim=1)))
         if (error < tol) exit
      end do
      if (present(iterations)) iterations = min(t, limit)
      if (present(residual)) residual = error
      status = nr_ok
      if (error >= tol) status = nr_iterations
   end subroutine max_ent

   pure subroutine next_uniform(state, cursor, value, status, uniforms)
      integer(int64), intent(inout) :: state !! Park-Miller state in 1:2147483646, ignored for supplied uniforms.
      integer, intent(inout) :: cursor !! Number of already consumed uniforms; advanced on successful consumption.
      real(dp), intent(out) :: value !! Next draw in [0,1), or NaN if the supplied stream is exhausted.
      integer, intent(out) :: status !! nr_ok or nr_random_exhausted; no implicit reuse of supplied draws.
      real(dp), intent(in), optional :: uniforms(:) !! Explicit finite [0,1) draws, consumed sequentially when supplied.
      status = nr_ok
      value = nan_value()
      if (present(uniforms)) then
         if (cursor >= size(uniforms)) then
            status = nr_random_exhausted
            return
         end if
         cursor = cursor + 1
         value = uniforms(cursor)
      else
         state = modulo(16807_int64 * state, 2147483647_int64)
         cursor = cursor + 1
         value = real(state, dp) / 2147483647.0_dp
      end if
   end subroutine next_uniform

   pure subroutine min_dens(rowsums, colsums, result, status, cost, lambda, k, alpha, delta, theta, &
      remove_prob, max_it, abs_tol, seed, uniforms, iterations, residual, draws_used)
      real(dp), intent(in) :: rowsums(:) !! Nonnegative finite target row totals, in vertex order.
      real(dp), intent(in) :: colsums(:) !! Nonnegative finite column totals, with the same length and total.
      real(dp), intent(out) :: result(:, :) !! Sparse zero-diagonal reconstruction; last accepted state on nonconvergence.
      integer, intent(out) :: status !! nr_ok or explicit shape, domain, infeasible, iteration or RNG-exhaustion code.
      real(dp), intent(in), optional :: cost !! Nonnegative cost of an added positive link; default one.
      real(dp), intent(in), optional :: lambda !! Initial allocation fraction in [0,1]; default one.
      integer, intent(in), optional :: k !! Number of iterations using lambda before full allocation; default 100.
      real(dp), intent(in), optional :: alpha !! Nonnegative scalar row-deficit weight, squared as in R; default 1/total.
      real(dp), intent(in), optional :: delta !! Nonnegative scalar column-deficit weight, squared as in R; default 1/total.
      real(dp), intent(in), optional :: theta !! Nonnegative Metropolis scaling of the objective change; default one.
      real(dp), intent(in), optional :: remove_prob !! Probability in [0,1] of random link removal; default 0.01.
      integer, intent(in), optional :: max_it !! Maximum add/remove attempts; default 100000, at least one.
      real(dp), intent(in), optional :: abs_tol !! Positive L1 row-deficit stopping tolerance; default 0.001.
      integer, intent(in), optional :: seed !! Local deterministic RNG seed, default 1; no global RNG state is modified.
      real(dp), intent(in), optional :: uniforms(:) !! Optional explicit finite [0,1) stream replacing local RNG, never recycled.
      integer, intent(out), optional :: iterations !! Number of attempted iterations, zero for invalid or zero-total margins.
      real(dp), intent(out), optional :: residual !! L1 sum of row and column deficits in returned estimate.
      integer, intent(out), optional :: draws_used !! Number of random draws consumed, including rejected proposals.
      integer :: n, limit, warm, t, i, j, p, chosen, nm, nv, count_draws, iter_done
      integer, allocatable :: available(:), occupied(:)
      real(dp), allocatable :: ad(:), ld(:), logs(:), probabilities(:)
      real(dp) :: total, tol, link_cost, fraction, aw, dw, temperature, remove, u, allocation, change, wsum, high, cumulative
      integer(int64) :: state
      logical :: remove_link
      result = nan_value()
      if (present(iterations)) iterations = 0
      if (present(residual)) residual = nan_value()
      if (present(draws_used)) draws_used = 0
      tol = 1.0e-3_dp
      if (present(abs_tol)) tol = abs_tol
      limit = 100000
      if (present(max_it)) limit = max_it
      warm = 100
      if (present(k)) warm = k
      link_cost = 1.0_dp
      if (present(cost)) link_cost = cost
      fraction = 1.0_dp
      if (present(lambda)) fraction = lambda
      temperature = 1.0_dp
      if (present(theta)) temperature = theta
      remove = 0.01_dp
      if (present(remove_prob)) remove = remove_prob
      status = nr_domain
      if (limit < 1 .or. warm < 0) return
      if (.not. all(ieee_is_finite([link_cost, fraction, temperature, remove]))) return
      if (link_cost < 0.0_dp .or. temperature < 0.0_dp) return
      if (fraction < 0.0_dp .or. fraction > 1.0_dp .or. remove < 0.0_dp .or. remove > 1.0_dp) return
      if (present(uniforms)) then
         if (any(.not. ieee_is_finite(uniforms))) return
         if (any(uniforms < 0.0_dp) .or. any(uniforms >= 1.0_dp)) return
      end if
      call check_margins(rowsums, colsums, n, total, status, tol)
      if (status /= nr_ok) return
      status = nr_shape
      if (any(shape(result) /= [n, n])) return
      aw = 0.0_dp
      if (total > 0.0_dp) aw = 1.0_dp / total
      dw = aw
      if (present(alpha)) aw = alpha
      if (present(delta)) dw = delta
      status = nr_domain
      if (.not. all(ieee_is_finite([aw, dw]))) return
      if (aw < 0.0_dp .or. dw < 0.0_dp) return
      if (max(aw, dw) > sqrt(huge(1.0_dp) / 4.0_dp)) return
      result = 0.0_dp
      status = nr_ok
      if (total == 0.0_dp) then
         if (present(residual)) residual = sum(colsums)
         return
      end if
      state = 1_int64
      if (present(seed)) state = modulo(int(seed, int64) - 1_int64, 2147483646_int64) + 1_int64
      allocate(available(n * (n - 1)), occupied(n * (n - 1)))
      allocate(logs(n * (n - 1)), probabilities(n * (n - 1)))
      ad = rowsums
      ld = colsums
      nm = 0
      do j = 1, n
         do i = 1, n
            if (i == j) cycle
            nm = nm + 1
            available(nm) = i + (j - 1) * n
         end do
      end do
      nv = 0
      count_draws = 0
      iter_done = 0
      status = nr_iterations
      do t = 1, limit
         iter_done = t
         if (t > warm) fraction = 1.0_dp
         probabilities = 0.0_dp
         logs = -huge(1.0_dp)
         do p = 1, nm
            i = modulo(available(p) - 1, n) + 1
            j = (available(p) - 1) / n + 1
            if (ad(i) > 0.0_dp .and. ld(j) > 0.0_dp) logs(p) = abs(log(ad(i)) - log(ld(j)))
         end do
         if (nm > 0) then
            high = maxval(logs(:nm))
            if (high > -huge(1.0_dp)) then
               do p = 1, nm
                  if (logs(p) > -huge(1.0_dp)) probabilities(p) = exp(logs(p) - high)
               end do
            end if
         end if
         wsum = sum(probabilities)
         call next_uniform(state, count_draws, u, status, uniforms)
         if (status /= nr_ok) exit
         remove_link = (u < remove .and. t > 1 .and. nv > 0) .or. wsum == 0.0_dp
         if (remove_link) then
            if (nv == 0) then
               status = nr_infeasible
               exit
            end if
            call next_uniform(state, count_draws, u, status, uniforms)
            if (status /= nr_ok) exit
            chosen = min(nv, 1 + int(u * real(nv, dp)))
            p = occupied(chosen)
            i = modulo(p - 1, n) + 1
            j = (p - 1) / n + 1
            ad(i) = ad(i) + result(i, j)
            ld(j) = ld(j) + result(i, j)
            result(i, j) = 0.0_dp
            nm = nm + 1
            available(nm) = p
            occupied(chosen:nv - 1) = occupied(chosen + 1:nv)
            nv = nv - 1
         else
            call next_uniform(state, count_draws, u, status, uniforms)
            if (status /= nr_ok) exit
            cumulative = 0.0_dp
            chosen = 0
            do p = 1, nm
               if (probabilities(p) <= 0.0_dp) cycle
               chosen = p
               cumulative = cumulative + probabilities(p)
               if (u * wsum < cumulative) exit
            end do
            p = available(chosen)
            i = modulo(p - 1, n) + 1
            j = (p - 1) / n + 1
            allocation = fraction * min(ad(i), ld(j))
            change = (aw * aw + dw * dw) * allocation
            if (allocation > 0.0_dp) change = change - link_cost
            if (.not. ieee_is_finite(change)) then
               status = nr_numeric
               exit
            end if
            call next_uniform(state, count_draws, u, status, uniforms)
            if (status /= nr_ok) exit
            if (change > 0.0_dp .or. u == 0.0_dp) then
               remove_link = .true.
            else
               remove_link = log(u) < temperature * change
            end if
            if (remove_link) then
               result(i, j) = allocation
               ad(i) = max(0.0_dp, ad(i) - allocation)
               ld(j) = max(0.0_dp, ld(j) - allocation)
               nv = nv + 1
               occupied(nv) = p
               available(chosen:nm - 1) = available(chosen + 1:nm)
               nm = nm - 1
            end if
         end if
         if (sum(abs(ad)) < tol) then
            status = nr_ok
            exit
         end if
         status = nr_iterations
      end do
      if (present(iterations)) iterations = iter_done
      if (present(residual)) residual = sum(abs(rowsums - sum(result, dim=2))) + &
         sum(abs(colsums - sum(result, dim=1)))
      if (present(draws_used)) draws_used = count_draws
   end subroutine min_dens

   pure subroutine matrix_estimation(rowsums, colsums, result, status, method, max_it, abs_tol, &
      cost, lambda, k, alpha, delta, theta, remove_prob, seed, uniforms, iterations, residual, draws_used)
      real(dp), intent(in) :: rowsums(:) !! Nonnegative row marginals of a square zero-diagonal network.
      real(dp), intent(in) :: colsums(:) !! Nonnegative column marginals in identical vertex order.
      real(dp), intent(out) :: result(:, :) !! Reconstructed n-by-n matrix; last iterate is retained on nonconvergence.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* diagnostic code.
      character(*), intent(in), optional :: method !! Reconstruction method: me (default) or md.
      integer, intent(in), optional :: max_it !! Maximum scaling sweeps or stochastic attempts, default 100000.
      real(dp), intent(in), optional :: abs_tol !! Positive absolute marginal-error stopping tolerance, default 0.001.
      real(dp), intent(in), optional :: cost !! Minimum-density positive-link cost; ignored for me.
      real(dp), intent(in), optional :: lambda !! Minimum-density initial allocation fraction in [0,1]; ignored for me.
      integer, intent(in), optional :: k !! Minimum-density fractional-allocation iteration count; ignored for me.
      real(dp), intent(in), optional :: alpha !! Minimum-density row-deficit coefficient, squared in the objective.
      real(dp), intent(in), optional :: delta !! Minimum-density column-deficit coefficient, squared in the objective.
      real(dp), intent(in), optional :: theta !! Minimum-density nonnegative Metropolis inverse-temperature scale.
      real(dp), intent(in), optional :: remove_prob !! Minimum-density random-removal probability in [0,1].
      integer, intent(in), optional :: seed !! Local reproducible minimum-density seed, default one.
      real(dp), intent(in), optional :: uniforms(:) !! Optional deterministic [0,1) stream for minimum-density reconstruction.
      integer, intent(out), optional :: iterations !! Number of executed scaling sweeps or stochastic attempts.
      real(dp), intent(out), optional :: residual !! L1 row-plus-column marginal error of the returned matrix.
      integer, intent(out), optional :: draws_used !! Consumed random draws for md; zero for me.
      character(:), allocatable :: choice
      choice = 'me'
      if (present(method)) choice = method
      if (present(draws_used)) draws_used = 0
      select case (choice)
      case ('me')
         call max_ent(rowsums, colsums, result, status, max_it, abs_tol, iterations, residual)
      case ('md')
         call min_dens(rowsums, colsums, result, status, cost, lambda, k, alpha, delta, theta, &
            remove_prob, max_it, abs_tol, seed, uniforms, iterations, residual, draws_used)
      case default
         result = nan_value()
         status = nr_domain
         if (present(iterations)) iterations = 0
         if (present(residual)) residual = nan_value()
      end select
   end subroutine matrix_estimation
end module network_reconstruction
