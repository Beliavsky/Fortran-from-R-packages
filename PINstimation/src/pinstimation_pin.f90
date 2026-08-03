! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_pin
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use pinstimation_kinds, only : dp, i8, pi
   use pinstimation_types, only : trade_counts, pin_parameters, pin_result, bayes_pin_result
   use pinstimation_math, only : logistic, logit, inv_softplus, softplus, log_sum_exp, poisson_log_pmf, &
      random_poisson, quantile_sorted
   use pinstimation_optimization, only : optimizer_result, minimize_nelder_mead, minimize_bfgs
   implicit none
   private
   public :: pin_loglik, pin_loglik_e, pin_loglik_lk, pin_loglik_eho
   public :: fit_pin, pin_ea, pin_gwj, pin_yz, pin_posteriors, pin_value
   public :: simulate_pin, fit_pin_bayes, initial_pin_ea, initial_pin_gwj

contains

   pure logical function valid_pin_parameters(p) result(ok)
      type(pin_parameters), intent(in) :: p
      ok = p%alpha >= 0.0_dp .and. p%alpha <= 1.0_dp .and. &
         p%delta >= 0.0_dp .and. p%delta <= 1.0_dp .and. &
         p%mu > 0.0_dp .and. p%eps_b > 0.0_dp .and. p%eps_s > 0.0_dp
   end function valid_pin_parameters

   pure real(dp) function pin_value(p) result(value)
      type(pin_parameters), intent(in) :: p
      value = p%alpha*p%mu/(p%alpha*p%mu + p%eps_b + p%eps_s)
   end function pin_value

   real(dp) function pin_loglik(data, p) result(value)
      type(trade_counts), intent(in) :: data
      type(pin_parameters), intent(in) :: p
      integer :: i
      real(dp) :: terms(3)
      if (.not. data%valid() .or. .not. valid_pin_parameters(p)) then
         value = -huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, data%size()
         terms(1) = safe_log(p%alpha*(1.0_dp - p%delta)) + &
            poisson_log_pmf(data%buys(i), p%eps_b + p%mu) + poisson_log_pmf(data%sells(i), p%eps_s)
         terms(2) = safe_log(p%alpha*p%delta) + poisson_log_pmf(data%buys(i), p%eps_b) + &
            poisson_log_pmf(data%sells(i), p%eps_s + p%mu)
         terms(3) = safe_log(1.0_dp - p%alpha) + poisson_log_pmf(data%buys(i), p%eps_b) + &
            poisson_log_pmf(data%sells(i), p%eps_s)
         value = value + log_sum_exp(terms)
      end do
   end function pin_loglik

   real(dp) function pin_loglik_e(data, p) result(value)
      type(trade_counts), intent(in) :: data
      type(pin_parameters), intent(in) :: p
      integer :: i
      real(dp) :: b, s, g1, g2, gmax, mix
      if (.not. data%valid() .or. .not. valid_pin_parameters(p)) then
         value = -huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, data%size()
         b = real(data%buys(i), dp)
         s = real(data%sells(i), dp)
         g1 = -p%mu + b*log(1.0_dp + p%mu/p%eps_b)
         g2 = -p%mu + s*log(1.0_dp + p%mu/p%eps_s)
         gmax = max(0.0_dp, max(g1, g2))
         mix = p%alpha*(1.0_dp - p%delta)*exp(g1 - gmax) + &
            p%alpha*p%delta*exp(g2 - gmax) + (1.0_dp - p%alpha)*exp(-gmax)
         if (mix <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
         end if
         value = value + b*log(p%eps_b) + s*log(p%eps_s) - p%eps_b - p%eps_s + &
            gmax + log(mix) - log_gamma(b + 1.0_dp) - log_gamma(s + 1.0_dp)
      end do
   end function pin_loglik_e

   real(dp) function pin_loglik_lk(data, p) result(value)
      type(trade_counts), intent(in) :: data
      type(pin_parameters), intent(in) :: p
      integer :: i
      real(dp) :: b, s, e1, e2, e3, emax, mix
      if (.not. data%valid() .or. .not. valid_pin_parameters(p)) then
         value = -huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, data%size()
         b = real(data%buys(i), dp)
         s = real(data%sells(i), dp)
         e1 = -p%mu - b*log(1.0_dp + p%mu/p%eps_b)
         e2 = -p%mu - s*log(1.0_dp + p%mu/p%eps_s)
         e3 = -b*log(1.0_dp + p%mu/p%eps_b) - s*log(1.0_dp + p%mu/p%eps_s)
         emax = max(e1, max(e2, e3))
         mix = p%alpha*p%delta*exp(e1 - emax) + p%alpha*(1.0_dp - p%delta)*exp(e2 - emax) + &
            (1.0_dp - p%alpha)*exp(e3 - emax)
         if (mix <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
         end if
         value = value + b*log(p%eps_b + p%mu) + s*log(p%eps_s + p%mu) - &
            (p%eps_b + p%eps_s) + emax + log(mix) - log_gamma(b + 1.0_dp) - log_gamma(s + 1.0_dp)
      end do
   end function pin_loglik_lk

   real(dp) function pin_loglik_eho(data, p) result(value)
      type(trade_counts), intent(in) :: data
      type(pin_parameters), intent(in) :: p
      integer :: i
      real(dp) :: b, s, m, xb, xs, p1, mix
      if (.not. data%valid() .or. .not. valid_pin_parameters(p)) then
         value = -huge(1.0_dp)
         return
      end if
      xb = p%eps_b/(p%mu + p%eps_b)
      xs = p%eps_s/(p%mu + p%eps_s)
      value = 0.0_dp
      do i = 1, data%size()
         b = real(data%buys(i), dp)
         s = real(data%sells(i), dp)
         m = min(b, s) + max(b, s)/2.0_dp
         p1 = b*log(p%mu + p%eps_b) + s*log(p%mu + p%eps_s) + m*log(xb) + m*log(xs) - &
            (p%eps_b + p%eps_s)
         mix = p%alpha*p%delta*exp(-p%mu)*xb**(b - m)*xs**(-m) + &
            p%alpha*(1.0_dp - p%delta)*exp(-p%mu)*xb**(-m)*xs**(s - m) + &
            (1.0_dp - p%alpha)*xb**(b - m)*xs**(s - m)
         if (mix <= 0.0_dp .or. .not. ieee_is_finite(mix)) then
            value = pin_loglik_e(data, p)
            return
         end if
         value = value + p1 + log(mix) - log_gamma(b + 1.0_dp) - log_gamma(s + 1.0_dp)
      end do
   end function pin_loglik_eho

   subroutine pin_posteriors(data, p, posterior)
      type(trade_counts), intent(in) :: data
      type(pin_parameters), intent(in) :: p
      real(dp), allocatable, intent(out) :: posterior(:,:)
      real(dp) :: terms(3), denominator
      integer :: i
      allocate(posterior(data%size(), 3))
      do i = 1, data%size()
         terms(1) = safe_log(1.0_dp - p%alpha) + poisson_log_pmf(data%buys(i), p%eps_b) + &
            poisson_log_pmf(data%sells(i), p%eps_s)
         terms(2) = safe_log(p%alpha*(1.0_dp - p%delta)) + &
            poisson_log_pmf(data%buys(i), p%eps_b + p%mu) + poisson_log_pmf(data%sells(i), p%eps_s)
         terms(3) = safe_log(p%alpha*p%delta) + poisson_log_pmf(data%buys(i), p%eps_b) + &
            poisson_log_pmf(data%sells(i), p%eps_s + p%mu)
         denominator = log_sum_exp(terms)
         posterior(i,:) = exp(terms - denominator)
      end do
   end subroutine pin_posteriors

   function initial_pin_ea(data) result(p)
      type(trade_counts), intent(in) :: data
      type(pin_parameters) :: p
      real(dp), allocatable :: b(:), s(:), imbalance(:)
      real(dp) :: mb, ms
      integer :: n
      n = data%size()
      allocate(b(n), s(n), imbalance(n))
      b = real(data%buys, dp)
      s = real(data%sells, dp)
      mb = sum(b)/real(max(1,n), dp)
      ms = sum(s)/real(max(1,n), dp)
      imbalance = abs(b - s)
      p%alpha = 0.5_dp
      p%delta = real(count(s > b), dp)/real(max(1,n), dp)
      p%delta = min(0.9_dp, max(0.1_dp, p%delta))
      p%mu = max(0.5_dp, sum(imbalance)/real(max(1,n), dp))
      p%eps_b = max(0.1_dp, mb - 0.5_dp*p%alpha*p%mu)
      p%eps_s = max(0.1_dp, ms - 0.5_dp*p%alpha*p%mu)
   end function initial_pin_ea

   function initial_pin_gwj(data) result(p)
      type(trade_counts), intent(in) :: data
      type(pin_parameters) :: p
      real(dp), allocatable :: b(:), s(:)
      real(dp) :: q_b, q_s, mb, ms
      integer :: n
      n = data%size()
      allocate(b(n), s(n))
      b = real(data%buys, dp)
      s = real(data%sells, dp)
      q_b = quantile_sorted(b, 0.25_dp)
      q_s = quantile_sorted(s, 0.25_dp)
      mb = sum(b)/real(max(1,n), dp)
      ms = sum(s)/real(max(1,n), dp)
      p%eps_b = max(0.1_dp, q_b)
      p%eps_s = max(0.1_dp, q_s)
      p%mu = max(0.5_dp, (mb - p%eps_b) + (ms - p%eps_s))
      p%alpha = min(0.9_dp, max(0.1_dp, ((mb - p%eps_b) + (ms - p%eps_s))/max(p%mu, 0.1_dp)))
      p%delta = min(0.9_dp, max(0.1_dp, (ms - p%eps_s)/max((mb - p%eps_b) + (ms - p%eps_s), 0.1_dp)))
   end function initial_pin_gwj

   subroutine fit_pin(data, result, initial, factorization, method, max_iterations, tolerance)
      type(trade_counts), intent(in) :: data
      type(pin_result), intent(out) :: result
      type(pin_parameters), intent(in), optional :: initial
      character(len=*), intent(in), optional :: factorization, method
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(pin_parameters) :: starts(8), candidate
      type(optimizer_result) :: opt, best
      character(len=16) :: fact, solver
      real(dp) :: u0(5)
      integer :: i, maxit
      real(dp) :: tol

      fact = 'E'
      solver = 'NELDER-MEAD'
      if (present(factorization)) fact = uppercase(trim(factorization))
      if (present(method)) solver = uppercase(trim(method))
      maxit = 2500
      tol = 1.0e-8_dp
      if (present(max_iterations)) maxit = max_iterations
      if (present(tolerance)) tol = tolerance
      starts(1) = initial_pin_ea(data)
      starts(2) = initial_pin_gwj(data)
      starts(3) = starts(1); starts(3)%alpha = 0.2_dp; starts(3)%delta = 0.2_dp
      starts(4) = starts(1); starts(4)%alpha = 0.2_dp; starts(4)%delta = 0.8_dp
      starts(5) = starts(1); starts(5)%alpha = 0.8_dp; starts(5)%delta = 0.2_dp
      starts(6) = starts(1); starts(6)%alpha = 0.8_dp; starts(6)%delta = 0.8_dp
      starts(7) = starts(2); starts(7)%mu = max(0.1_dp, 0.5_dp*starts(2)%mu)
      starts(8) = starts(2); starts(8)%mu = 2.0_dp*starts(2)%mu
      if (present(initial)) starts(1) = initial
      best%objective = huge(1.0_dp)
      do i = 1, size(starts)
         call pin_pack(starts(i), u0)
         if (solver(1:min(len_trim(solver),4)) == 'BFGS') then
            call minimize_bfgs(objective, u0, opt, tol, maxit)
         else
            call minimize_nelder_mead(objective, u0, opt, tol, maxit)
         end if
         if (opt%objective < best%objective) best = opt
      end do
      call pin_unpack(best%parameters, candidate)
      result%parameters = candidate
      result%log_likelihood = -best%objective
      result%pin = pin_value(candidate)
      result%pin_good = result%pin*(1.0_dp - candidate%delta)
      result%pin_bad = result%pin*candidate%delta
      result%iterations = best%iterations
      result%evaluations = best%evaluations
      result%status = best%status
      result%converged = best%converged
      call pin_posteriors(data, candidate, result%posteriors)

   contains
      real(dp) function objective(u) result(v)
         real(dp), intent(in) :: u(:)
         type(pin_parameters) :: p
         call pin_unpack(u, p)
         select case (trim(fact))
         case ('LK')
            v = -pin_loglik_lk(data, p)
         case ('EHO')
            v = -pin_loglik_eho(data, p)
         case ('NONE','DIRECT')
            v = -pin_loglik(data, p)
         case default
            v = -pin_loglik_e(data, p)
         end select
         if (.not. ieee_is_finite(v)) v = huge(1.0_dp)/100.0_dp
      end function objective
   end subroutine fit_pin

   subroutine pin_ea(data, result, max_iterations)
      type(trade_counts), intent(in) :: data
      type(pin_result), intent(out) :: result
      integer, intent(in), optional :: max_iterations
      type(pin_parameters) :: p
      p = initial_pin_ea(data)
      call fit_pin(data, result, initial=p, factorization='E', max_iterations=max_iterations)
   end subroutine pin_ea

   subroutine pin_gwj(data, result, max_iterations)
      type(trade_counts), intent(in) :: data
      type(pin_result), intent(out) :: result
      integer, intent(in), optional :: max_iterations
      type(pin_parameters) :: p
      p = initial_pin_gwj(data)
      call fit_pin(data, result, initial=p, factorization='E', max_iterations=max_iterations)
   end subroutine pin_gwj

   subroutine pin_yz(data, result, grid_size, max_iterations)
      type(trade_counts), intent(in) :: data
      type(pin_result), intent(out) :: result
      integer, intent(in), optional :: grid_size, max_iterations
      type(pin_parameters) :: p
      integer :: g
      g = 5
      if (present(grid_size)) g = max(2, grid_size)
      p = initial_pin_ea(data)
      p%alpha = 1.0_dp/(2.0_dp*real(g, dp))
      p%delta = p%alpha
      call fit_pin(data, result, initial=p, factorization='E', max_iterations=max_iterations)
   end subroutine pin_yz

   subroutine simulate_pin(days, p, data, states, seed, status)
      integer, intent(in) :: days
      type(pin_parameters), intent(in) :: p
      type(trade_counts), intent(out) :: data
      integer, allocatable, intent(out), optional :: states(:)
      integer, intent(in), optional :: seed
      integer, intent(out), optional :: status
      integer :: i, local_state, stat1, stat2
      real(dp) :: u, lambda_b, lambda_s
      if (present(status)) status = 0
      if (days < 1 .or. .not. valid_pin_parameters(p)) then
         if (present(status)) status = 1
         allocate(data%buys(0), data%sells(0))
         return
      end if
      if (present(seed)) call seed_local(seed)
      allocate(data%buys(days), data%sells(days))
      if (present(states)) allocate(states(days))
      do i = 1, days
         call random_number(u)
         if (u < 1.0_dp - p%alpha) then
            local_state = 0
            lambda_b = p%eps_b
            lambda_s = p%eps_s
         else
            call random_number(u)
            if (u < 1.0_dp - p%delta) then
               local_state = 1
               lambda_b = p%eps_b + p%mu
               lambda_s = p%eps_s
            else
               local_state = 2
               lambda_b = p%eps_b
               lambda_s = p%eps_s + p%mu
            end if
         end if
         data%buys(i) = random_poisson(lambda_b, stat1)
         data%sells(i) = random_poisson(lambda_s, stat2)
         if (present(states)) states(i) = local_state
         if (present(status)) status = max(status, max(stat1, stat2))
      end do
   end subroutine simulate_pin

   subroutine fit_pin_bayes(data, result, sweeps, burnin, thin, initial, proposal_scale, seed)
      type(trade_counts), intent(in) :: data
      type(bayes_pin_result), intent(out) :: result
      integer, intent(in), optional :: sweeps, burnin, thin, seed
      type(pin_parameters), intent(in), optional :: initial
      real(dp), intent(in), optional :: proposal_scale
      type(pin_result) :: mle
      type(pin_parameters) :: p
      real(dp) :: current(5), proposal(5), lp_current, lp_proposal, scale, u
      integer :: nsweep, nburn, nthin, nkeep, i, kept, accepted
      nsweep = 2000; nburn = 500; nthin = 1; scale = 0.12_dp
      if (present(sweeps)) nsweep = sweeps
      if (present(burnin)) nburn = burnin
      if (present(thin)) nthin = max(1, thin)
      if (present(proposal_scale)) scale = proposal_scale
      if (present(seed)) call seed_local(seed)
      if (present(initial)) then
         p = initial
      else
         call fit_pin(data, mle, max_iterations=1200)
         p = mle%parameters
      end if
      call pin_pack(p, current)
      lp_current = log_posterior(current)
      nkeep = max(0, (nsweep - nburn)/nthin)
      allocate(result%draws(nkeep, 5))
      kept = 0; accepted = 0
      do i = 1, nsweep
         proposal = current + scale*normal_vector(5)
         lp_proposal = log_posterior(proposal)
         call random_number(u)
         if (log(max(u, tiny(1.0_dp))) < lp_proposal - lp_current) then
            current = proposal
            lp_current = lp_proposal
            accepted = accepted + 1
         end if
         if (i > nburn .and. modulo(i - nburn, nthin) == 0) then
            kept = kept + 1
            call pin_unpack(current, p)
            result%draws(kept,:) = [p%alpha, p%delta, p%mu, p%eps_b, p%eps_s]
         end if
      end do
      if (nkeep > 0) then
         result%posterior_mean%alpha = sum(result%draws(:,1))/real(nkeep, dp)
         result%posterior_mean%delta = sum(result%draws(:,2))/real(nkeep, dp)
         result%posterior_mean%mu = sum(result%draws(:,3))/real(nkeep, dp)
         result%posterior_mean%eps_b = sum(result%draws(:,4))/real(nkeep, dp)
         result%posterior_mean%eps_s = sum(result%draws(:,5))/real(nkeep, dp)
      end if
      result%acceptance_rate = real(accepted, dp)/real(max(1,nsweep), dp)

   contains
      real(dp) function log_posterior(uvec) result(v)
         real(dp), intent(in) :: uvec(:)
         type(pin_parameters) :: pp
         call pin_unpack(uvec, pp)
         v = pin_loglik_e(data, pp) - 0.5_dp*sum((uvec/5.0_dp)**2)
      end function log_posterior
   end subroutine fit_pin_bayes

   subroutine pin_pack(p, u)
      type(pin_parameters), intent(in) :: p
      real(dp), intent(out) :: u(5)
      u = [logit(p%alpha), logit(p%delta), inv_softplus(p%mu), inv_softplus(p%eps_b), inv_softplus(p%eps_s)]
   end subroutine pin_pack

   subroutine pin_unpack(u, p)
      real(dp), intent(in) :: u(:)
      type(pin_parameters), intent(out) :: p
      p%alpha = logistic(u(1))
      p%delta = logistic(u(2))
      p%mu = softplus(u(3)) + 1.0e-10_dp
      p%eps_b = softplus(u(4)) + 1.0e-10_dp
      p%eps_s = softplus(u(5)) + 1.0e-10_dp
   end subroutine pin_unpack

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
         put(i) = modulo(seed + 104729*i, huge(1) - 1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_local

   function normal_vector(n) result(z)
      integer, intent(in) :: n
      real(dp) :: z(n), u1, u2
      integer :: i
      i = 1
      do while (i <= n)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         z(i) = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
         if (i + 1 <= n) z(i + 1) = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
         i = i + 2
      end do
   end function normal_vector

end module pinstimation_pin
