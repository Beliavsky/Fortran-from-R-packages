! SPDX-License-Identifier: GPL-2.0-or-later
module ceoptim_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ceoptim_kinds, only : dp
   use ceoptim_rng, only : rng_state, rng_seed, rng_categorical
   use ceoptim_sampling, only : rtmvnorm, tmvn_result
   use ceoptim_types, only : ce_control, ce_continuous_control, ce_discrete_control, &
      ce_result, ce_state, ce_objective
   implicit none
   private
   public :: ce_optimize

contains

   subroutine ce_optimize(f, result, control, continuous, discrete)
      procedure(ce_objective) :: f
      type(ce_result), intent(out) :: result
      type(ce_control), intent(in), optional :: control
      type(ce_continuous_control), intent(in), optional :: continuous
      type(ce_discrete_control), intent(in), optional :: discrete

      type(ce_control) :: ctl
      type(rng_state) :: rng
      type(tmvn_result) :: tmv
      type(ce_state), allocatable :: state_work(:)
      real(dp), allocatable :: xc(:, :), y(:), mu(:), sd(:), sigma(:, :)
      real(dp), allocatable :: tau(:, :), tau0(:, :), empirical(:, :), ceprocess(:)
      integer, allocatable :: xd(:, :), order(:), categories(:), best_xd(:)
      real(dp), allocatable :: best_xc(:)
      integer :: p, q, maxcat, nelite, iter, idx, state_count
      integer :: actual_nfe
      real(dp) :: alpha, beta, gamma, eps_sd, eta, sgn
      real(dp) :: optimum, gammat, diffopt, max_sd, max_prob_dev
      logical :: has_c, has_d, stagnated

      ctl = ce_control()
      if (present(control)) ctl = control
      result%status = 0
      result%message = 'ok'

      has_c = present(continuous)
      has_d = present(discrete)
      p = 0
      q = 0
      maxcat = 0

      if (has_c) then
         if (.not. allocated(continuous%mean) .or. .not. allocated(continuous%sd)) then
            call fail(result, 1, 'continuous mean and sd must both be supplied')
            return
         end if
         p = size(continuous%mean)
         if (size(continuous%sd) /= p .or. p < 1) then
            call fail(result, 2, 'continuous mean and sd have invalid dimensions')
            return
         end if
         if (any(.not. ieee_is_finite(continuous%mean)) .or. &
             any(.not. ieee_is_finite(continuous%sd)) .or. any(continuous%sd < 0.0_dp)) then
            call fail(result, 3, 'continuous mean/sd must be finite and sd nonnegative')
            return
         end if
         if (allocated(continuous%con_mat) .neqv. allocated(continuous%con_vec)) then
            call fail(result, 4, 'con_mat and con_vec must be supplied together')
            return
         end if
         if (allocated(continuous%con_mat)) then
            if (size(continuous%con_mat, 2) /= p .or. &
                size(continuous%con_mat, 1) /= size(continuous%con_vec)) then
               call fail(result, 5, 'continuous constraint dimensions are not conformable')
               return
            end if
         end if
         alpha = continuous%smooth_mean
         beta = continuous%smooth_sd
         eps_sd = continuous%sd_thr
      else
         alpha = 1.0_dp
         beta = 1.0_dp
         eps_sd = 0.001_dp
      end if

      if (has_d) then
         if (.not. allocated(discrete%categories)) then
            call fail(result, 6, 'discrete categories must be supplied')
            return
         end if
         q = size(discrete%categories)
         if (q < 1 .or. any(discrete%categories < 1)) then
            call fail(result, 7, 'all category counts must be positive')
            return
         end if
         maxcat = maxval(discrete%categories)
         allocate(categories(q))
         categories = discrete%categories
         gamma = discrete%smooth_prob
         eta = discrete%prob_thr
      else
         gamma = 1.0_dp
         eta = 0.001_dp
         allocate(categories(0))
      end if

      if (p + q == 0) then
         call fail(result, 8, 'at least one continuous or discrete variable is required')
         return
      end if
      if (ctl%n < 1 .or. ctl%iter_thr < 0 .or. ctl%no_improve_thr < 1) then
         call fail(result, 9, 'N, iter_thr, or no_improve_thr is invalid')
         return
      end if
      if (ctl%rho <= 0.0_dp .or. ctl%rho >= 1.0_dp) then
         call fail(result, 10, 'rho must lie strictly between zero and one')
         return
      end if
      if (min(alpha, beta, gamma) < 0.0_dp .or. max(alpha, beta, gamma) > 1.0_dp) then
         call fail(result, 11, 'smoothing parameters must lie in [0,1]')
         return
      end if

      nelite = nint(real(ctl%n, dp) * ctl%rho)
      if (nelite < 1 .or. nelite > ctl%n) then
         call fail(result, 12, 'N*rho rounds to an invalid elite sample size')
         return
      end if
      if (p > 0 .and. nelite < 2) then
         call fail(result, 13, 'continuous optimization requires at least two elite samples for sd estimation')
         return
      end if

      call rng_seed(rng, ctl%seed)
      allocate(xc(ctl%n, p), xd(ctl%n, q), y(ctl%n), order(ctl%n))
      allocate(best_xc(p), best_xd(q))
      if (p > 0) then
         allocate(mu(p), sd(p), sigma(p, p))
         mu = continuous%mean
         sd = continuous%sd
         call diagonal_covariance(sd, sigma)
         call sample_continuous(ctl%n, mu, sigma, continuous, rng, tmv)
         if (tmv%status /= 0) then
            call fail(result, 20 + tmv%status, 'initial truncated-normal sampling failed: '//tmv%message)
            return
         end if
         xc = tmv%x
      end if

      if (q > 0) then
         allocate(tau(maxcat, q), tau0(maxcat, q), empirical(maxcat, q))
         call initialize_probabilities(discrete, categories, tau0, result)
         if (result%status /= 0) return
         tau = tau0
         call sample_discrete(ctl%n, categories, tau, rng, xd)
      end if

      sgn = 1.0_dp
      if (ctl%maximize) sgn = -1.0_dp
      call evaluate_population(f, xc, xd, sgn, y)
      actual_nfe = ctl%n
      call argsort(y, order)

      if (p > 0) then
         call elite_continuous_stats(xc, order, nelite, mu, sd)
         mu = alpha * mu + (1.0_dp - alpha) * continuous%mean
         sd = beta * sd + (1.0_dp - beta) * continuous%sd
         call diagonal_covariance(sd, sigma)
      end if
      if (q > 0) then
         call elite_discrete_probs(xd, order, nelite, categories, empirical)
         tau = gamma * empirical + (1.0_dp - gamma) * tau0
      end if

      idx = order(1)
      if (p > 0) best_xc = xc(idx, :)
      if (q > 0) best_xd = xd(idx, :)
      optimum = y(idx)
      gammat = y(order(nelite))

      allocate(state_work(max(0, ctl%iter_thr)), ceprocess(max(1, ctl%iter_thr)))
      iter = 0
      state_count = 0
      diffopt = huge(1.0_dp)
      stagnated = .false.
      max_sd = current_max_sd(sd, p)
      max_prob_dev = current_max_prob_dev(tau, categories, q)

      do while (iter < ctl%iter_thr .and. .not. stagnated .and. &
                distributions_active(p, q, max_sd, eps_sd, max_prob_dev, eta))
         state_count = state_count + 1
         call store_state(state_work(state_count), iter, sgn * optimum, sgn * gammat, &
            mu, tau, p, q, max_sd, max_prob_dev)

         if (ctl%verbose) call print_progress(iter, sgn * optimum, max_sd, max_prob_dev, p, q)

         if (p > 0) then
            call sample_continuous(ctl%n, mu, sigma, continuous, rng, tmv)
            if (tmv%status /= 0) then
               call fail(result, 40 + tmv%status, 'truncated-normal sampling failed: '//tmv%message)
               return
            end if
            xc = tmv%x
         end if
         if (q > 0) call sample_discrete(ctl%n, categories, tau, rng, xd)

         call evaluate_population(f, xc, xd, sgn, y)
         actual_nfe = actual_nfe + ctl%n
         call argsort(y, order)

         if (y(order(1)) < optimum) then
            if (p > 0) best_xc = xc(order(1), :)
            if (q > 0) best_xd = xd(order(1), :)
            optimum = y(order(1))
            if (y(order(nelite)) < gammat) gammat = y(order(nelite))
         end if

         ceprocess(iter + 1) = optimum
         if (iter > ctl%no_improve_thr) then
            diffopt = sum(abs(ceprocess(iter - ctl%no_improve_thr:iter - 1) - optimum))
            stagnated = diffopt <= tiny(1.0_dp)
         end if

         if (p > 0) then
            call elite_continuous_stats(xc, order, nelite, mu0=mu, sd0=sd, &
               alpha=alpha, beta=beta)
            call diagonal_covariance(sd, sigma)
         end if
         if (q > 0) then
            tau0 = tau
            call elite_discrete_probs(xd, order, nelite, categories, empirical)
            tau = gamma * empirical + (1.0_dp - gamma) * tau0
         end if

         iter = iter + 1
         max_sd = current_max_sd(sd, p)
         max_prob_dev = current_max_prob_dev(tau, categories, q)
      end do

      allocate(result%continuous(p), result%discrete(q), result%categories(q))
      result%continuous = best_xc
      result%discrete = best_xd
      result%categories = categories
      result%optimum = sgn * optimum
      result%niter = iter
      result%nfe = iter * ctl%n
      result%actual_nfe = actual_nfe
      if (iter == ctl%iter_thr) then
         result%convergence = 'Not converged'
      else if (stagnated) then
         result%convergence = 'Optimum did not change for the requested number of iterations'
      else
         result%convergence = 'Variance converged'
      end if
      allocate(result%states(state_count))
      if (state_count > 0) result%states = state_work(1:state_count)
      allocate(result%final_mean(p), result%final_sd(p))
      if (p > 0) then
         result%final_mean = mu
         result%final_sd = sd
      end if
      if (q > 0) then
         allocate(result%final_probs(maxcat, q))
         result%final_probs = tau
      else
         allocate(result%final_probs(0, 0))
      end if
      result%status = 0
      result%message = 'ok'
   end subroutine ce_optimize

   subroutine initialize_probabilities(discrete, categories, tau, result)
      type(ce_discrete_control), intent(in) :: discrete
      integer, intent(in) :: categories(:)
      real(dp), intent(out) :: tau(:, :)
      type(ce_result), intent(inout) :: result
      integer :: j, k
      real(dp) :: total

      tau = 0.0_dp
      if (allocated(discrete%probs)) then
         if (size(discrete%probs, 2) /= size(categories) .or. &
             size(discrete%probs, 1) < maxval(categories)) then
            call fail(result, 14, 'discrete probs has invalid dimensions')
            return
         end if
         do j = 1, size(categories)
            k = categories(j)
            if (any(discrete%probs(1:k, j) < 0.0_dp)) then
               call fail(result, 15, 'discrete probabilities must be nonnegative')
               return
            end if
            total = sum(discrete%probs(1:k, j))
            if (total <= 0.0_dp) then
               call fail(result, 16, 'each discrete probability column must have positive mass')
               return
            end if
            tau(1:k, j) = discrete%probs(1:k, j) / total
         end do
      else
         do j = 1, size(categories)
            k = categories(j)
            tau(1:k, j) = 1.0_dp / real(k, dp)
         end do
      end if
   end subroutine initialize_probabilities

   subroutine sample_continuous(n, mu, sigma, continuous, rng, tmv)
      integer, intent(in) :: n
      real(dp), intent(in) :: mu(:), sigma(:, :)
      type(ce_continuous_control), intent(in) :: continuous
      type(rng_state), intent(inout) :: rng
      type(tmvn_result), intent(out) :: tmv
      if (allocated(continuous%con_mat)) then
         call rtmvnorm(n, mu, sigma, rng, tmv, a=continuous%con_mat, b=continuous%con_vec)
      else
         call rtmvnorm(n, mu, sigma, rng, tmv)
      end if
   end subroutine sample_continuous

   subroutine sample_discrete(n, categories, tau, rng, xd)
      integer, intent(in) :: n, categories(:)
      real(dp), intent(in) :: tau(:, :)
      type(rng_state), intent(inout) :: rng
      integer, intent(out) :: xd(:, :)
      integer :: i, j
      do j = 1, size(categories)
         do i = 1, n
            xd(i, j) = rng_categorical(rng, tau(:, j), categories(j))
         end do
      end do
   end subroutine sample_discrete

   subroutine evaluate_population(f, xc, xd, sgn, y)
      procedure(ce_objective) :: f
      real(dp), intent(in) :: xc(:, :)
      integer, intent(in) :: xd(:, :)
      real(dp), intent(in) :: sgn
      real(dp), intent(out) :: y(:)
      integer :: i
      real(dp) :: v
      do i = 1, size(y)
         v = f(xc(i, :), xd(i, :))
         if (ieee_is_finite(v)) then
            y(i) = sgn * v
         else
            y(i) = huge(1.0_dp)
         end if
      end do
   end subroutine evaluate_population

   subroutine elite_continuous_stats(xc, order, nelite, mu, sd, mu0, sd0, alpha, beta)
      real(dp), intent(in) :: xc(:, :)
      integer, intent(in) :: order(:), nelite
      real(dp), intent(out), optional :: mu(:), sd(:)
      real(dp), intent(inout), optional :: mu0(:), sd0(:)
      real(dp), intent(in), optional :: alpha, beta
      real(dp), allocatable :: m(:), s(:)
      real(dp) :: d
      integer :: p, j, ii

      p = size(xc, 2)
      allocate(m(p), s(p))
      m = 0.0_dp
      do ii = 1, nelite
         m = m + xc(order(ii), :)
      end do
      m = m / real(nelite, dp)
      s = 0.0_dp
      do ii = 1, nelite
         do j = 1, p
            d = xc(order(ii), j) - m(j)
            s(j) = s(j) + d*d
         end do
      end do
      s = sqrt(s / real(nelite - 1, dp))

      if (present(mu)) mu = m
      if (present(sd)) sd = s
      if (present(mu0)) mu0 = alpha * m + (1.0_dp - alpha) * mu0
      if (present(sd0)) sd0 = beta * s + (1.0_dp - beta) * sd0
   end subroutine elite_continuous_stats

   subroutine elite_discrete_probs(xd, order, nelite, categories, probs)
      integer, intent(in) :: xd(:, :), order(:), nelite, categories(:)
      real(dp), intent(out) :: probs(:, :)
      integer :: j, ii, k
      probs = 0.0_dp
      do j = 1, size(categories)
         do ii = 1, nelite
            k = xd(order(ii), j) + 1
            if (k >= 1 .and. k <= categories(j)) probs(k, j) = probs(k, j) + 1.0_dp
         end do
         probs(1:categories(j), j) = probs(1:categories(j), j) / real(nelite, dp)
      end do
   end subroutine elite_discrete_probs

   subroutine diagonal_covariance(sd, sigma)
      real(dp), intent(in) :: sd(:)
      real(dp), intent(out) :: sigma(:, :)
      integer :: j
      sigma = 0.0_dp
      do j = 1, size(sd)
         sigma(j, j) = sd(j)*sd(j)
      end do
   end subroutine diagonal_covariance

   subroutine argsort(x, order)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: order(:)
      integer :: i, j, key
      do i = 1, size(x)
         order(i) = i
      end do
      do i = 2, size(x)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (x(order(j)) <= x(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine argsort

   pure real(dp) function current_max_sd(sd, p) result(v)
      real(dp), allocatable, intent(in) :: sd(:)
      integer, intent(in) :: p
      if (p > 0) then
         v = maxval(sd)
      else
         v = 0.0_dp
      end if
   end function current_max_sd

   pure real(dp) function current_max_prob_dev(tau, categories, q) result(v)
      real(dp), allocatable, intent(in) :: tau(:, :)
      integer, intent(in) :: categories(:), q
      integer :: j
      v = 0.0_dp
      if (q > 0) then
         do j = 1, q
            v = max(v, 1.0_dp - maxval(tau(1:categories(j), j)))
         end do
      end if
   end function current_max_prob_dev

   pure logical function distributions_active(p, q, max_sd, eps_sd, max_prob, eta) result(active)
      integer, intent(in) :: p, q
      real(dp), intent(in) :: max_sd, eps_sd, max_prob, eta
      active = (p > 0 .and. max_sd > eps_sd) .or. (q > 0 .and. max_prob > eta)
   end function distributions_active

   subroutine store_state(state, iter, optimum, gammat, mu, tau, p, q, max_sd, max_prob)
      type(ce_state), intent(inout) :: state
      integer, intent(in) :: iter, p, q
      real(dp), intent(in) :: optimum, gammat, max_sd, max_prob
      real(dp), allocatable, intent(in) :: mu(:), tau(:, :)
      state%iter = iter
      state%optimum = optimum
      state%gammat = gammat
      state%max_sd = max_sd
      state%max_prob_dev = max_prob
      allocate(state%mean(p))
      if (p > 0) state%mean = mu
      if (q > 0) then
         allocate(state%probs(size(tau, 1), q))
         state%probs = tau
      else
         allocate(state%probs(0, 0))
      end if
   end subroutine store_state

   subroutine print_progress(iter, optimum, max_sd, max_prob, p, q)
      integer, intent(in) :: iter, p, q
      real(dp), intent(in) :: optimum, max_sd, max_prob
      if (p > 0 .and. q > 0) then
         write(*, '(a,i0,a,es14.6,a,es12.4,a,es12.4)') &
            'iter: ', iter, ' opt: ', optimum, ' maxSd: ', max_sd, ' maxProbs: ', max_prob
      else if (p > 0) then
         write(*, '(a,i0,a,es14.6,a,es12.4)') 'iter: ', iter, ' opt: ', optimum, ' maxSd: ', max_sd
      else
         write(*, '(a,i0,a,es14.6,a,es12.4)') 'iter: ', iter, ' opt: ', optimum, ' maxProbs: ', max_prob
      end if
   end subroutine print_progress

   subroutine fail(result, code, message)
      type(ce_result), intent(inout) :: result
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      result%status = code
      result%message = message
      result%convergence = 'Error'
   end subroutine fail

end module ceoptim_core
