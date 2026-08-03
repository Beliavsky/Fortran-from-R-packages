! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_mpin
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use pinstimation_kinds, only : dp, i8
   use pinstimation_types, only : trade_counts, mpin_parameters, mpin_result
   use pinstimation_math, only : logistic, logit, softplus, inv_softplus, log_sum_exp, poisson_log_pmf, &
      random_poisson, quantile_sorted
   use pinstimation_optimization, only : optimizer_result, minimize_nelder_mead, minimize_bfgs
   implicit none
   private
   public :: mpin_loglik, mpin_posteriors, mpin_value, fit_mpin_ml, fit_mpin_ecm
   public :: initial_mpin, simulate_mpin, detectlayers_e, detectlayers_eg, detectlayers_ecm

contains

   pure logical function valid_mpin_parameters(p) result(ok)
      type(mpin_parameters), intent(in) :: p
      integer :: j
      j = p%layers()
      ok = j > 0 .and. allocated(p%delta) .and. allocated(p%mu)
      if (.not. ok) return
      ok = size(p%delta) == j .and. size(p%mu) == j
      if (.not. ok) return
      ok = all(p%alpha >= 0.0_dp) .and. sum(p%alpha) < 1.0_dp .and. &
         all(p%delta >= 0.0_dp) .and. all(p%delta <= 1.0_dp) .and. &
         all(p%mu > 0.0_dp) .and. p%eps_b > 0.0_dp .and. p%eps_s > 0.0_dp
   end function valid_mpin_parameters

   pure real(dp) function mpin_value(p) result(value)
      type(mpin_parameters), intent(in) :: p
      real(dp) :: informed
      informed = sum(p%alpha*p%mu)
      value = informed/(informed + p%eps_b + p%eps_s)
   end function mpin_value

   real(dp) function mpin_loglik(data, p) result(value)
      type(trade_counts), intent(in) :: data
      type(mpin_parameters), intent(in) :: p
      real(dp), allocatable :: terms(:)
      real(dp) :: b, s, base
      integer :: i, j, layers
      if (.not. data%valid() .or. .not. valid_mpin_parameters(p)) then
         value = -huge(1.0_dp)
         return
      end if
      layers = p%layers()
      allocate(terms(2*layers + 1))
      value = 0.0_dp
      do i = 1, data%size()
         b = real(data%buys(i), dp)
         s = real(data%sells(i), dp)
         base = b*log(p%eps_b) + s*log(p%eps_s) - p%eps_b - p%eps_s - &
            log_gamma(b + 1.0_dp) - log_gamma(s + 1.0_dp)
         terms(1) = safe_log(1.0_dp - sum(p%alpha))
         do j = 1, layers
            terms(2*j) = safe_log(p%alpha(j)*(1.0_dp - p%delta(j))) - p%mu(j) + &
               b*log(1.0_dp + p%mu(j)/p%eps_b)
            terms(2*j + 1) = safe_log(p%alpha(j)*p%delta(j)) - p%mu(j) + &
               s*log(1.0_dp + p%mu(j)/p%eps_s)
         end do
         value = value + base + log_sum_exp(terms)
      end do
   end function mpin_loglik

   subroutine mpin_posteriors(data, p, posterior)
      type(trade_counts), intent(in) :: data
      type(mpin_parameters), intent(in) :: p
      real(dp), allocatable, intent(out) :: posterior(:,:)
      real(dp), allocatable :: terms(:)
      real(dp) :: denominator
      integer :: i, j, layers
      layers = p%layers()
      allocate(posterior(data%size(), 2*layers + 1), terms(2*layers + 1))
      do i = 1, data%size()
         terms(1) = safe_log(1.0_dp - sum(p%alpha)) + poisson_log_pmf(data%buys(i), p%eps_b) + &
            poisson_log_pmf(data%sells(i), p%eps_s)
         do j = 1, layers
            terms(2*j) = safe_log(p%alpha(j)*(1.0_dp - p%delta(j))) + &
               poisson_log_pmf(data%buys(i), p%eps_b + p%mu(j)) + poisson_log_pmf(data%sells(i), p%eps_s)
            terms(2*j + 1) = safe_log(p%alpha(j)*p%delta(j)) + poisson_log_pmf(data%buys(i), p%eps_b) + &
               poisson_log_pmf(data%sells(i), p%eps_s + p%mu(j))
         end do
         denominator = log_sum_exp(terms)
         posterior(i,:) = exp(terms - denominator)
      end do
   end subroutine mpin_posteriors

   function initial_mpin(data, layers) result(p)
      type(trade_counts), intent(in) :: data
      integer, intent(in) :: layers
      type(mpin_parameters) :: p
      real(dp), allocatable :: b(:), s(:), aoi(:)
      real(dp) :: mb, ms, delta0, q
      integer :: j, n
      n = data%size()
      allocate(b(n), s(n), aoi(n), p%alpha(layers), p%delta(layers), p%mu(layers))
      b = real(data%buys, dp)
      s = real(data%sells, dp)
      aoi = abs(b - s)
      mb = sum(b)/real(max(1,n), dp)
      ms = sum(s)/real(max(1,n), dp)
      delta0 = min(0.9_dp, max(0.1_dp, real(count(s > b), dp)/real(max(1,n), dp)))
      p%alpha = 0.55_dp/real(layers, dp)
      p%delta = delta0
      do j = 1, layers
         q = real(j, dp)/real(layers + 1, dp)
         p%mu(j) = max(0.5_dp, quantile_sorted(aoi, q))
         if (j > 1) p%mu(j) = max(p%mu(j), p%mu(j - 1) + 0.25_dp)
      end do
      p%eps_b = max(0.1_dp, mb - 0.5_dp*sum(p%alpha*p%mu))
      p%eps_s = max(0.1_dp, ms - 0.5_dp*sum(p%alpha*p%mu))
   end function initial_mpin

   subroutine fit_mpin_ml(data, layers, result, initial, method, max_iterations, tolerance)
      type(trade_counts), intent(in) :: data
      integer, intent(in) :: layers
      type(mpin_result), intent(out) :: result
      type(mpin_parameters), intent(in), optional :: initial
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(mpin_parameters) :: p0, p, starts(4)
      type(optimizer_result) :: opt, best
      real(dp), allocatable :: u0(:)
      character(len=16) :: solver
      integer :: i, maxit
      real(dp) :: tol
      if (layers < 1 .or. .not. data%valid()) then
         result%status = 2
         return
      end if
      solver = 'NELDER-MEAD'
      if (present(method)) solver = uppercase(trim(method))
      maxit = 3500; tol = 1.0e-8_dp
      if (present(max_iterations)) maxit = max_iterations
      if (present(tolerance)) tol = tolerance
      p0 = initial_mpin(data, layers)
      if (present(initial)) p0 = sorted_mpin(initial)
      starts(1) = p0
      starts(2) = p0; starts(2)%alpha = 0.35_dp/real(layers,dp)
      starts(3) = p0; starts(3)%alpha = 0.75_dp/real(layers,dp); starts(3)%mu = 0.7_dp*p0%mu
      starts(4) = p0; starts(4)%delta = 1.0_dp - p0%delta; starts(4)%mu = 1.4_dp*p0%mu
      best%objective = huge(1.0_dp)
      do i = 1, size(starts)
         call mpin_pack(starts(i), u0)
         if (solver(1:min(4,len_trim(solver))) == 'BFGS') then
            call minimize_bfgs(objective, u0, opt, tol, maxit)
         else
            call minimize_nelder_mead(objective, u0, opt, tol, maxit)
         end if
         if (opt%objective < best%objective) best = opt
      end do
      call mpin_unpack(best%parameters, layers, p)
      call fill_result(data, p, best, result)

   contains
      real(dp) function objective(u) result(v)
         real(dp), intent(in) :: u(:)
         type(mpin_parameters) :: pp
         call mpin_unpack(u, layers, pp)
         v = -mpin_loglik(data, pp)
         if (.not. ieee_is_finite(v)) v = huge(1.0_dp)/100.0_dp
      end function objective
   end subroutine fit_mpin_ml

   subroutine fit_mpin_ecm(data, layers, result, initial, max_iterations, tolerance)
      type(trade_counts), intent(in) :: data
      integer, intent(in) :: layers
      type(mpin_result), intent(out) :: result
      type(mpin_parameters), intent(in), optional :: initial
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(mpin_parameters) :: p
      type(optimizer_result) :: inner
      real(dp), allocatable :: posterior(:,:), rates0(:)
      real(dp) :: old_ll, new_ll, tol, good, bad
      integer :: iter, j, maxit
      if (layers < 1 .or. .not. data%valid()) then
         result%status = 2
         return
      end if
      p = initial_mpin(data, layers)
      if (present(initial)) p = sorted_mpin(initial)
      maxit = 200; tol = 1.0e-7_dp
      if (present(max_iterations)) maxit = max_iterations
      if (present(tolerance)) tol = tolerance
      old_ll = mpin_loglik(data, p)
      do iter = 1, maxit
         call mpin_posteriors(data, p, posterior)
         do j = 1, layers
            good = sum(posterior(:,2*j))
            bad = sum(posterior(:,2*j + 1))
            p%alpha(j) = max(1.0e-8_dp, (good + bad)/real(data%size(), dp))
            p%delta(j) = min(1.0_dp - 1.0e-8_dp, max(1.0e-8_dp, bad/max(good + bad, 1.0e-12_dp)))
         end do
         if (sum(p%alpha) >= 1.0_dp) p%alpha = 0.999_dp*p%alpha/sum(p%alpha)
         call pack_rates(p, rates0)
         call minimize_nelder_mead(q_objective, rates0, inner, 1.0e-7_dp, 700)
         call unpack_rates(inner%parameters, p)
         new_ll = mpin_loglik(data, p)
         if (abs(new_ll - old_ll) <= tol*(1.0_dp + abs(old_ll))) exit
         old_ll = new_ll
      end do
      result%parameters = p
      result%log_likelihood = mpin_loglik(data, p)
      result%iterations = iter
      result%evaluations = inner%evaluations
      result%converged = iter <= maxit
      result%status = merge(0, 1, result%converged)
      call finish_mpin_result(data, result)

   contains
      real(dp) function q_objective(u) result(v)
         real(dp), intent(in) :: u(:)
         type(mpin_parameters) :: pp
         integer :: i, jj
         real(dp) :: q
         pp = p
         call unpack_rates(u, pp)
         q = 0.0_dp
         do i = 1, data%size()
            q = q + posterior(i,1)*(poisson_log_pmf(data%buys(i), pp%eps_b) + &
               poisson_log_pmf(data%sells(i), pp%eps_s))
            do jj = 1, layers
               q = q + posterior(i,2*jj)*(poisson_log_pmf(data%buys(i), pp%eps_b + pp%mu(jj)) + &
                  poisson_log_pmf(data%sells(i), pp%eps_s))
               q = q + posterior(i,2*jj + 1)*(poisson_log_pmf(data%buys(i), pp%eps_b) + &
                  poisson_log_pmf(data%sells(i), pp%eps_s + pp%mu(jj)))
            end do
         end do
         v = -q
         if (.not. ieee_is_finite(v)) v = huge(1.0_dp)/100.0_dp
      end function q_objective
   end subroutine fit_mpin_ecm

   subroutine simulate_mpin(days, p, data, states, seed, status)
      integer, intent(in) :: days
      type(mpin_parameters), intent(in) :: p
      type(trade_counts), intent(out) :: data
      integer, allocatable, intent(out), optional :: states(:)
      integer, intent(in), optional :: seed
      integer, intent(out), optional :: status
      real(dp), allocatable :: cprob(:)
      real(dp) :: u, lambda_b, lambda_s
      integer :: i, j, state, stat1, stat2, layers
      if (present(status)) status = 0
      if (days < 1 .or. .not. valid_mpin_parameters(p)) then
         allocate(data%buys(0), data%sells(0))
         if (present(status)) status = 1
         return
      end if
      if (present(seed)) call seed_local(seed)
      layers = p%layers()
      allocate(cprob(2*layers + 1), data%buys(days), data%sells(days))
      if (present(states)) allocate(states(days))
      cprob(1) = 1.0_dp - sum(p%alpha)
      do j = 1, layers
         cprob(2*j) = p%alpha(j)*(1.0_dp - p%delta(j))
         cprob(2*j + 1) = p%alpha(j)*p%delta(j)
      end do
      do i = 2, size(cprob)
         cprob(i) = cprob(i) + cprob(i - 1)
      end do
      do i = 1, days
         call random_number(u)
         state = 1
         do while (state < size(cprob) .and. u > cprob(state))
            state = state + 1
         end do
         lambda_b = p%eps_b; lambda_s = p%eps_s
         if (state > 1) then
            j = state/2
            if (modulo(state,2) == 0) then
               lambda_b = lambda_b + p%mu(j)
            else
               lambda_s = lambda_s + p%mu(j)
            end if
         end if
         data%buys(i) = random_poisson(lambda_b, stat1)
         data%sells(i) = random_poisson(lambda_s, stat2)
         if (present(states)) states(i) = state - 1
         if (present(status)) status = max(status, max(stat1, stat2))
      end do
   end subroutine simulate_mpin

   integer function detectlayers_e(data, confidence, max_layers) result(layers)
      type(trade_counts), intent(in) :: data
      real(dp), intent(in), optional :: confidence
      integer, intent(in), optional :: max_layers
      real(dp), allocatable :: x(:), sorted(:), gaps(:)
      real(dp) :: threshold, conf
      integer :: maxl, n
      conf = 0.995_dp; maxl = 8
      if (present(confidence)) conf = confidence
      if (present(max_layers)) maxl = max_layers
      n = data%size()
      if (n < 3) then
         layers = 1
         return
      end if
      allocate(x(n), sorted(n), gaps(n-1))
      x = abs(real(data%buys - data%sells, dp))
      sorted = x
      call sort_local(sorted)
      gaps = sorted(2:) - sorted(:n-1)
      threshold = quantile_sorted(gaps, min(0.999_dp, max(0.5_dp, conf)))
      layers = 1 + count(gaps > max(threshold, 0.5_dp))
      layers = min(maxl, max(1, layers))
   end function detectlayers_e

   integer function detectlayers_eg(data, max_layers) result(layers)
      type(trade_counts), intent(in) :: data
      integer, intent(in), optional :: max_layers
      real(dp), allocatable :: x(:)
      real(dp) :: bic, best_bic
      integer :: k, maxl
      maxl = 6
      if (present(max_layers)) maxl = max_layers
      allocate(x(data%size()))
      x = abs(real(data%buys - data%sells, dp))
      best_bic = huge(1.0_dp); layers = 1
      do k = 1, min(maxl, max(1,data%size()/3))
         bic = kmeans_bic(x, k)
         if (bic < best_bic) then
            best_bic = bic
            layers = k
         end if
      end do
   end function detectlayers_eg

   integer function detectlayers_ecm(data, max_layers) result(layers)
      type(trade_counts), intent(in) :: data
      integer, intent(in), optional :: max_layers
      type(mpin_result) :: fit
      real(dp) :: bic, best_bic
      integer :: k, maxl, npar
      maxl = 4
      if (present(max_layers)) maxl = max_layers
      best_bic = huge(1.0_dp); layers = 1
      do k = 1, maxl
         call fit_mpin_ecm(data, k, fit, max_iterations=80, tolerance=1.0e-5_dp)
         npar = 3*k + 2
         bic = -2.0_dp*fit%log_likelihood + real(npar,dp)*log(real(max(1,data%size()),dp))
         if (bic < best_bic) then
            best_bic = bic
            layers = k
         end if
      end do
   end function detectlayers_ecm

   subroutine fill_result(data, p, opt, result)
      type(trade_counts), intent(in) :: data
      type(mpin_parameters), intent(in) :: p
      type(optimizer_result), intent(in) :: opt
      type(mpin_result), intent(out) :: result
      result%parameters = p
      result%log_likelihood = -opt%objective
      result%iterations = opt%iterations
      result%evaluations = opt%evaluations
      result%status = opt%status
      result%converged = opt%converged
      call finish_mpin_result(data, result)
   end subroutine fill_result

   subroutine finish_mpin_result(data, result)
      type(trade_counts), intent(in) :: data
      type(mpin_result), intent(inout) :: result
      integer :: layers
      real(dp) :: denominator
      layers = result%parameters%layers()
      allocate(result%mpin_layer(layers), result%good_layer(layers), result%bad_layer(layers))
      denominator = sum(result%parameters%alpha*result%parameters%mu) + result%parameters%eps_b + result%parameters%eps_s
      result%mpin_layer = result%parameters%alpha*result%parameters%mu/denominator
      result%good_layer = result%mpin_layer*(1.0_dp - result%parameters%delta)
      result%bad_layer = result%mpin_layer*result%parameters%delta
      result%mpin = sum(result%mpin_layer)
      call mpin_posteriors(data, result%parameters, result%posteriors)
   end subroutine finish_mpin_result

   subroutine mpin_pack(p_input, u)
      type(mpin_parameters), intent(in) :: p_input
      real(dp), allocatable, intent(out) :: u(:)
      type(mpin_parameters) :: p
      real(dp) :: noinfo, increment
      integer :: j, layers
      p = sorted_mpin(p_input)
      layers = p%layers()
      allocate(u(3*layers + 2))
      noinfo = max(1.0e-12_dp, 1.0_dp - sum(p%alpha))
      do j = 1, layers
         u(j) = log(max(p%alpha(j),1.0e-12_dp)/noinfo)
         u(layers + j) = logit(p%delta(j))
         if (j == 1) then
            increment = p%mu(j)
         else
            increment = p%mu(j) - p%mu(j - 1)
         end if
         u(2*layers + j) = inv_softplus(max(increment,1.0e-8_dp))
      end do
      u(3*layers + 1) = inv_softplus(p%eps_b)
      u(3*layers + 2) = inv_softplus(p%eps_s)
   end subroutine mpin_pack

   subroutine mpin_unpack(u, layers, p)
      real(dp), intent(in) :: u(:)
      integer, intent(in) :: layers
      type(mpin_parameters), intent(out) :: p
      real(dp), allocatable :: ex(:)
      real(dp) :: m, denom
      integer :: j
      allocate(p%alpha(layers), p%delta(layers), p%mu(layers), ex(layers))
      m = max(0.0_dp, maxval(u(1:layers)))
      ex = exp(u(1:layers) - m)
      denom = exp(-m) + sum(ex)
      p%alpha = ex/denom
      do j = 1, layers
         p%delta(j) = logistic(u(layers + j))
         if (j == 1) then
            p%mu(j) = softplus(u(2*layers + j)) + 1.0e-10_dp
         else
            p%mu(j) = p%mu(j - 1) + softplus(u(2*layers + j)) + 1.0e-10_dp
         end if
      end do
      p%eps_b = softplus(u(3*layers + 1)) + 1.0e-10_dp
      p%eps_s = softplus(u(3*layers + 2)) + 1.0e-10_dp
   end subroutine mpin_unpack

   subroutine pack_rates(p, u)
      type(mpin_parameters), intent(in) :: p
      real(dp), allocatable, intent(out) :: u(:)
      integer :: j, layers
      real(dp) :: increment
      layers = p%layers()
      allocate(u(layers + 2))
      do j = 1, layers
         if (j == 1) then
            increment = p%mu(j)
         else
            increment = p%mu(j) - p%mu(j - 1)
         end if
         u(j) = inv_softplus(max(increment, 1.0e-8_dp))
      end do
      u(layers+1) = inv_softplus(p%eps_b)
      u(layers+2) = inv_softplus(p%eps_s)
   end subroutine pack_rates

   subroutine unpack_rates(u, p)
      real(dp), intent(in) :: u(:)
      type(mpin_parameters), intent(inout) :: p
      integer :: j, layers
      layers = p%layers()
      do j = 1, layers
         if (j == 1) then
            p%mu(j) = softplus(u(j)) + 1.0e-10_dp
         else
            p%mu(j) = p%mu(j-1) + softplus(u(j)) + 1.0e-10_dp
         end if
      end do
      p%eps_b = softplus(u(layers+1)) + 1.0e-10_dp
      p%eps_s = softplus(u(layers+2)) + 1.0e-10_dp
   end subroutine unpack_rates

   function sorted_mpin(p_input) result(p)
      type(mpin_parameters), intent(in) :: p_input
      type(mpin_parameters) :: p
      integer, allocatable :: order(:)
      integer :: i, j, key, layers
      layers = p_input%layers()
      allocate(order(layers), p%alpha(layers), p%delta(layers), p%mu(layers))
      order = [(i, i=1,layers)]
      do i = 2, layers
         key = order(i); j = i - 1
         do while (j >= 1)
            if (p_input%mu(order(j)) <= p_input%mu(key)) exit
            order(j+1) = order(j); j = j - 1
         end do
         order(j+1) = key
      end do
      p%alpha = p_input%alpha(order)
      p%delta = p_input%delta(order)
      p%mu = p_input%mu(order)
      p%eps_b = p_input%eps_b; p%eps_s = p_input%eps_s
   end function sorted_mpin

   real(dp) function kmeans_bic(x, k) result(bic)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable :: center(:), newcenter(:), sums(:)
      integer, allocatable :: group(:), counts(:)
      real(dp) :: sse, best, d
      integer :: i, j, iter, n
      n = size(x)
      allocate(center(k), newcenter(k), sums(k), group(n), counts(k))
      do j = 1, k
         center(j) = quantile_sorted(x, (real(j,dp)-0.5_dp)/real(k,dp))
      end do
      do iter = 1, 100
         do i = 1, n
            group(i) = 1; best = abs(x(i)-center(1))
            do j = 2, k
               d = abs(x(i)-center(j))
               if (d < best) then
                  best = d; group(i) = j
               end if
            end do
         end do
         sums = 0.0_dp; counts = 0
         do i = 1, n
            sums(group(i)) = sums(group(i)) + x(i)
            counts(group(i)) = counts(group(i)) + 1
         end do
         newcenter = center
         do j = 1, k
            if (counts(j) > 0) newcenter(j) = sums(j)/real(counts(j),dp)
         end do
         if (maxval(abs(newcenter-center)) < 1.0e-8_dp) exit
         center = newcenter
      end do
      sse = 0.0_dp
      do i = 1, n
         sse = sse + (x(i)-newcenter(group(i)))**2
      end do
      bic = real(n,dp)*log(max(sse/real(max(1,n),dp),1.0e-12_dp)) + real(2*k,dp)*log(real(max(2,n),dp))
   end function kmeans_bic

   subroutine sort_local(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i); j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j); j = j - 1
         end do
         x(j+1) = key
      end do
   end subroutine sort_local

   pure real(dp) function safe_log(x) result(v)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         v = log(x)
      else
         v = -huge(1.0_dp)
      end if
   end function safe_log

   pure function uppercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) then
            out(i:i) = achar(code - 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function uppercase

   subroutine seed_local(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: n, i
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 130363*i, huge(1) - 1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_local

end module pinstimation_mpin
