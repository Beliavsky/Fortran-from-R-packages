! SPDX-License-Identifier: GPL-3.0-or-later
module loo_evaluation
   use loo_types
   use loo_importance, only: importance_sampling, importance_weights
   use matrixstats, only: log_sum_exp, col_vars
   use posterior_diagnostics, only: ess_mean
   use r_transforms, only: r_expm1, r_log1mexp, r_log1pexp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_negative_inf
   implicit none
   private
   public :: loo, loo_matrix, loo_array, loo_function, loo_i
   public :: waic, waic_matrix, waic_array, waic_function, elpd
   public :: loo_approximate_posterior, loo_approximate_posterior_matrix
   public :: loo_approximate_posterior_array, loo_approximate_posterior_function
   public :: relative_eff, relative_eff_function, compute_point_estimate, thin_draws
   interface loo
      module procedure loo_matrix, loo_array
   end interface loo
   interface waic
      module procedure waic_matrix, waic_array
   end interface waic
   interface elpd
      module procedure elpd_matrix, elpd_array
   end interface elpd
   interface loo_approximate_posterior
      module procedure loo_approximate_posterior_matrix, loo_approximate_posterior_array
   end interface loo_approximate_posterior
   interface relative_eff
      module procedure relative_eff_matrix, relative_eff_array, relative_eff_vector
   end interface relative_eff
contains
   pure function loo_matrix(log_lik, r_eff, is_method) result(out)
      real(dp), intent(in) :: log_lik(:, :) !! Finite log likelihoods, posterior draws by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiencies, scalar array or one per observation; default one.
      character(*), intent(in), optional :: is_method !! 'psis' (default), 'sis', or 'tis'.
      type(loo_result) :: out
      character(len = 4) :: method
      method = 'psis'
      if (present(is_method)) method = is_method
      if (any(.not. ieee_is_finite(log_lik))) return
      out%importance = importance_sampling(-log_lik, trim(method), r_eff)
      call finish_loo(log_lik, out)
   end function loo_matrix

   pure subroutine finish_loo(log_lik, out)
      real(dp), intent(in) :: log_lik(:, :) !! Finite draw-by-observation log likelihoods matching stored importance weights.
      type(loo_result), intent(inout) :: out !! Result containing successful importance weights; populated with pointwise estimates.
      real(dp), allocatable :: lw(:, :)
      real(dp) :: elpd_i, lpd_i
      integer :: j, n, s
      if (out%importance%status /= 0) return
      s = size(log_lik, 1)
      n = size(log_lik, 2)
      lw = importance_weights(out%importance)
      allocate(out%pointwise(n, 5))
      do j = 1, n
         elpd_i = log_sum_exp(log_lik(:, j) + lw(:, j))
         lpd_i = log_sum_exp(log_lik(:, j)) - log(real(s, dp))
         out%pointwise(j, 1) = elpd_i
         out%pointwise(j, 2) = mcse_elpd(log_lik(:, j), lw(:, j), elpd_i, out%importance%r_eff(j))
         out%pointwise(j, 3) = lpd_i - elpd_i
         out%pointwise(j, 4) = -2.0_dp * elpd_i
         out%pointwise(j, 5) = out%importance%pareto_k(j)
      end do
      out%estimates = summarize_pointwise(out%pointwise(:, [1, 3, 4]))
      out%n_draws = s
      out%method = 'loo'
      out%status = 0
   end subroutine finish_loo

   pure function mcse_elpd(ll, lw, elpd_value, r_eff) result(mcse)
      real(dp), intent(in) :: ll(:) !! Finite log likelihoods for a single observation.
      real(dp), intent(in) :: lw(:) !! Normalized log importance weights matching likelihood draws.
      real(dp), intent(in) :: elpd_value !! Log of the importance-weighted predictive density.
      real(dp), intent(in) :: r_eff !! Positive relative MCMC efficiency for the observation.
      real(dp) :: mcse, d, log_abs, terms(size(ll)), log_variance
      integer :: i
      terms = ieee_value(0.0_dp, ieee_negative_inf)
      do i = 1, size(ll)
         d = ll(i) - elpd_value
         if (d == 0.0_dp) cycle
         if (d > 0.0_dp) then
            log_abs = d + r_log1mexp(-d)
         else
            log_abs = log(-r_expm1(d))
         end if
         terms(i) = 2.0_dp * (lw(i) + log_abs)
      end do
      log_variance = log_sum_exp(terms) - log(r_eff)
      mcse = sqrt(r_log1pexp(log_variance))
   end function mcse_elpd

   pure function waic_matrix(log_lik) result(out)
      real(dp), intent(in) :: log_lik(:, :) !! Finite posterior-draw by observation log likelihoods; at least two draws.
      type(loo_result) :: out
      integer :: n, s, j
      n = size(log_lik, 2)
      s = size(log_lik, 1)
      if (n < 1 .or. s < 2) return
      if (any(.not. ieee_is_finite(log_lik))) return
      allocate(out%pointwise(n, 3))
      out%pointwise(:, 2) = col_vars(log_lik)
      do j = 1, n
         out%pointwise(j, 1) = log_sum_exp(log_lik(:, j)) - log(real(s, dp)) - out%pointwise(j, 2)
      end do
      out%pointwise(:, 3) = -2.0_dp * out%pointwise(:, 1)
      out%estimates = summarize_pointwise(out%pointwise)
      out%n_draws = s
      out%method = 'waic'
      out%status = 0
   end function waic_matrix

   pure function elpd_matrix(log_lik) result(out)
      real(dp), intent(in) :: log_lik(:, :) !! Finite posterior-draw by observation log likelihoods.
      type(loo_result) :: out
      integer :: n, s, j
      n = size(log_lik, 2)
      s = size(log_lik, 1)
      if (n < 1 .or. s < 1) return
      if (any(.not. ieee_is_finite(log_lik))) return
      allocate(out%pointwise(n, 2))
      do j = 1, n
         out%pointwise(j, 1) = log_sum_exp(log_lik(:, j)) - log(real(s, dp))
      end do
      out%pointwise(:, 2) = -2.0_dp * out%pointwise(:, 1)
      out%estimates = summarize_pointwise(out%pointwise)
      out%n_draws = s
      out%method = 'elpd'
      out%status = 0
   end function elpd_matrix

   pure function loo_approximate_posterior_matrix(log_lik, log_p, log_g) result(out)
      real(dp), intent(in) :: log_lik(:, :) !! Finite log likelihoods, approximate-posterior draws by observations.
      real(dp), intent(in) :: log_p(:) !! Finite target log posterior density for each draw, up to an additive constant.
      real(dp), intent(in) :: log_g(:) !! Finite approximating log density for each draw, up to an additive constant.
      type(loo_result) :: out
      real(dp), allocatable :: ratios(:, :)
      integer :: j, s
      s = size(log_lik, 1)
      if (s < 1 .or. size(log_lik, 2) < 1) return
      if (size(log_p) /= s .or. size(log_g) /= s) return
      if (any(.not. ieee_is_finite(log_lik))) return
      if (any(.not. ieee_is_finite(log_p)) .or. any(.not. ieee_is_finite(log_g))) return
      ratios = -log_lik
      do j = 1, size(log_lik, 2)
         ratios(:, j) = ratios(:, j) + log_p - log_g
         if (any(.not. ieee_is_finite(ratios(:, j)))) return
         ratios(:, j) = ratios(:, j) - maxval(ratios(:, j))
      end do
      out%importance = importance_sampling(ratios, 'psis')
      call finish_loo(log_lik, out)
   end function loo_approximate_posterior_matrix

   pure function loo_array(log_lik, r_eff, is_method) result(out)
      real(dp), intent(in) :: log_lik(:, :, :) !! Log likelihoods with shape (iterations, chains, observations).
      real(dp), intent(in), optional :: r_eff(:) !! Scalar array or per-observation relative efficiency; default one.
      character(*), intent(in), optional :: is_method !! 'psis', 'tis', or 'sis', default 'psis'.
      type(loo_result) :: out
      out = loo_matrix(reshape(log_lik, &
         [size(log_lik, 1) * size(log_lik, 2), size(log_lik, 3)]), r_eff, is_method)
   end function loo_array

   pure function waic_array(log_lik) result(out)
      real(dp), intent(in) :: log_lik(:, :, :) !! Log likelihoods with shape (iterations, chains, observations).
      type(loo_result) :: out
      out = waic_matrix(reshape(log_lik, &
         [size(log_lik, 1) * size(log_lik, 2), size(log_lik, 3)]))
   end function waic_array

   pure function elpd_array(log_lik) result(out)
      real(dp), intent(in) :: log_lik(:, :, :) !! Log likelihoods with shape (iterations, chains, observations).
      type(loo_result) :: out
      out = elpd_matrix(reshape(log_lik, &
         [size(log_lik, 1) * size(log_lik, 2), size(log_lik, 3)]))
   end function elpd_array

   pure function loo_approximate_posterior_array(log_lik, log_p, log_g) result(out)
      real(dp), intent(in) :: log_lik(:, :, :) !! Log likelihoods with shape (iterations, chains, observations).
      real(dp), intent(in) :: log_p(:) !! Target log density in chain-major draw order.
      real(dp), intent(in) :: log_g(:) !! Approximating log density in chain-major draw order.
      type(loo_result) :: out
      out = loo_approximate_posterior_matrix(reshape(log_lik, &
         [size(log_lik, 1) * size(log_lik, 2), size(log_lik, 3)]), log_p, log_g)
   end function loo_approximate_posterior_array

   function loo_function(llfun, n_draws, n_observations, r_eff, is_method) result(out)
      procedure(log_likelihood_callback), pointer, intent(in) :: llfun !! Callback filling a draw vector for a one-based observation index.
      integer, intent(in) :: n_draws !! Positive number of posterior draws returned by each callback.
      integer, intent(in) :: n_observations !! Positive number of observations to evaluate.
      real(dp), intent(in), optional :: r_eff(:) !! Scalar array or per-observation relative efficiency; default one.
      character(*), intent(in), optional :: is_method !! Importance method, 'psis' by default.
      type(loo_result) :: out
      real(dp), allocatable :: log_lik(:, :)
      integer :: j
      if (.not. associated(llfun)) return
      if (n_draws < 1 .or. n_observations < 1) return
      allocate(log_lik(n_draws, n_observations))
      do j = 1, n_observations
         call llfun(j, log_lik(:, j))
      end do
      out = loo_matrix(log_lik, r_eff, is_method)
   end function loo_function

   function waic_function(llfun, n_draws, n_observations) result(out)
      procedure(log_likelihood_callback), pointer, intent(in) :: llfun !! Callback filling a draw vector for a one-based observation index.
      integer, intent(in) :: n_draws !! Positive number of posterior draws returned by each callback.
      integer, intent(in) :: n_observations !! Positive number of observations to evaluate.
      type(loo_result) :: out
      real(dp), allocatable :: log_lik(:, :)
      integer :: j
      if (.not. associated(llfun)) return
      if (n_draws < 1 .or. n_observations < 1) return
      allocate(log_lik(n_draws, n_observations))
      do j = 1, n_observations
         call llfun(j, log_lik(:, j))
      end do
      out = waic_matrix(log_lik)
   end function waic_function

   function loo_approximate_posterior_function(llfun, n_draws, n_observations, log_p, log_g) result(out)
      procedure(log_likelihood_callback), pointer, intent(in) :: llfun !! Callback filling a draw vector for a one-based observation index.
      integer, intent(in) :: n_draws !! Positive number of posterior draws returned by each callback.
      integer, intent(in) :: n_observations !! Positive number of observations to evaluate.
      real(dp), intent(in) :: log_p(:) !! Target log density for each of n_draws draws.
      real(dp), intent(in) :: log_g(:) !! Approximating log density for the same draws.
      type(loo_result) :: out
      real(dp), allocatable :: log_lik(:, :)
      integer :: j
      if (.not. associated(llfun)) return
      if (n_draws < 1 .or. n_observations < 1) return
      allocate(log_lik(n_draws, n_observations))
      do j = 1, n_observations
         call llfun(j, log_lik(:, j))
      end do
      out = loo_approximate_posterior_matrix(log_lik, log_p, log_g)
   end function loo_approximate_posterior_function

   function loo_i(observation, llfun, n_draws, r_eff, is_method) result(out)
      integer, intent(in) :: observation !! Positive observation index sent unchanged to the callback.
      procedure(log_likelihood_callback), pointer, intent(in) :: llfun !! Callback returning log likelihoods for the requested observation.
      integer, intent(in) :: n_draws !! Positive number of posterior draws in the callback output.
      real(dp), intent(in), optional :: r_eff(:) !! One-element relative-efficiency array; default one.
      character(*), intent(in), optional :: is_method !! 'psis', 'tis', or 'sis', default 'psis'.
      type(loo_result) :: out
      real(dp), allocatable :: ll(:, :)
      if (.not. associated(llfun)) return
      if (n_draws < 1 .or. observation < 1) return
      allocate(ll(n_draws, 1))
      call llfun(observation, ll(:, 1))
      out = loo_matrix(ll, r_eff, is_method)
   end function loo_i

   pure function relative_eff_array(likelihood) result(eff)
      real(dp), intent(in) :: likelihood(:, :, :) !! Unlogged likelihood values with shape (iterations, chains, observations).
      real(dp) :: eff(size(likelihood, 3))
      integer :: j, s
      s = size(likelihood, 1) * size(likelihood, 2)
      eff = loo_nan()
      if (s < 1) return
      do j = 1, size(likelihood, 3)
         eff(j) = ess_mean(likelihood(:, :, j)) / real(s, dp)
      end do
   end function relative_eff_array

   pure function relative_eff_matrix(likelihood, chain_id) result(eff)
      real(dp), intent(in) :: likelihood(:, :) !! Unlogged likelihood matrix with posterior draws in rows.
      integer, intent(in) :: chain_id(:) !! Contiguous positive chain IDs 1 through C, equal counts; interleaving is allowed.
      real(dp) :: eff(size(likelihood, 2))
      real(dp), allocatable :: chains(:, :, :)
      integer :: c, n, iter, i, at
      eff = loo_nan()
      n = size(likelihood, 1)
      if (size(chain_id) /= n .or. n < 1) return
      if (any(chain_id < 1)) return
      c = maxval(chain_id)
      if (c > n) return
      if (mod(n, c) /= 0) return
      iter = n / c
      do i = 1, c
         if (count(chain_id == i) /= iter) return
      end do
      allocate(chains(iter, c, size(likelihood, 2)))
      do i = 1, c
         at = 0
         do n = 1, size(chain_id)
            if (chain_id(n) /= i) cycle
            at = at + 1
            chains(at, i, :) = likelihood(n, :)
         end do
      end do
      eff = relative_eff_array(chains)
   end function relative_eff_matrix

   pure function relative_eff_vector(likelihood, chain_id) result(eff)
      real(dp), intent(in) :: likelihood(:) !! Unlogged likelihood draws for one observation.
      integer, intent(in) :: chain_id(:) !! Chain IDs corresponding to likelihood draws, with equal chain lengths.
      real(dp) :: eff, values(1)
      values = relative_eff_matrix(reshape(likelihood, [size(likelihood), 1]), chain_id)
      eff = values(1)
   end function relative_eff_vector

   function relative_eff_function(likelihood_fun, n_draws, n_observations, chain_id) result(eff)
      procedure(log_likelihood_callback), pointer, intent(in) :: likelihood_fun !! Callback filling unlogged likelihoods despite the generic callback type name.
      integer, intent(in) :: n_draws !! Positive number of likelihood draws per observation.
      integer, intent(in) :: n_observations !! Nonnegative number of observations to evaluate.
      integer, intent(in) :: chain_id(:) !! Chain ID for every draw, using equally sized chains numbered from one.
      real(dp) :: eff(max(0, n_observations))
      real(dp), allocatable :: values(:)
      integer :: j
      eff = loo_nan()
      if (.not. associated(likelihood_fun)) return
      if (n_draws < 1 .or. n_observations < 1) return
      allocate(values(n_draws))
      do j = 1, n_observations
         call likelihood_fun(j, values)
         eff(j) = relative_eff_vector(values, chain_id)
      end do
   end function relative_eff_function

   pure function compute_point_estimate(draws) result(estimate)
      real(dp), intent(in) :: draws(:, :) !! Posterior draws by model parameters; nonempty first dimension required.
      real(dp) :: estimate(size(draws, 2))
      estimate = loo_nan()
      if (size(draws, 1) > 0) estimate = sum(draws, dim = 1) / real(size(draws, 1), dp)
   end function compute_point_estimate

   pure function thin_draws(draws, n_keep) result(thinned)
      real(dp), intent(in) :: draws(:, :) !! Posterior draws by model parameters, retained in original order.
      integer, intent(in) :: n_keep !! Number of equally spaced retained draws, from one through the original draw count.
      real(dp), allocatable :: thinned(:, :)
      integer :: i, index
      allocate(thinned(0, size(draws, 2)))
      if (n_keep < 1 .or. n_keep > size(draws, 1)) return
      deallocate(thinned)
      allocate(thinned(n_keep, size(draws, 2)))
      do i = 1, n_keep
         index = i * (size(draws, 1) / n_keep)
         thinned(i, :) = draws(index, :)
      end do
   end function thin_draws
end module loo_evaluation
