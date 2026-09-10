! SPDX-License-Identifier: GPL-2.0-only
module marss_constrained_em
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_fit_result, marss_kf_result, marss_hatyt_result, marss_constraints
   use marss_kalman, only : marss_kfss
   use marss_analysis, only : marss_hatyt
   use marss_constraints_mod, only : marss_constraints_valid, marss_constraints_parameter_count
   use marss_constraints_mod, only : marss_constraint_start_vector, marss_apply_constraints
   use marss_parameters, only : marss_b_at, marss_u_at, marss_q_at, marss_z_at, marss_a_at, marss_r_at
   use marss_parameters, only : marss_v0_effective
   use marss_covariance, only : psd_inverse_logdet
   use marss_utils, only : identity_matrix
   implicit none
   private
   public :: marss_kem_linear

contains

   pure subroutine marss_kem_linear(start_model, constraints, max_iter, tol, fit, mstep_iter)
      type(marss_model), intent(in) :: start_model !! Template model supplying data and unconstrained blocks.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta constraints estimated by generalized EM.
      integer, intent(in) :: max_iter !! Maximum number of generalized-EM iterations; must be positive.
      real(dp), intent(in) :: tol !! Relative observed-log-likelihood convergence tolerance; must be positive.
      type(marss_fit_result), intent(out) :: fit !! Constrained EM result, beta vector, and convergence diagnostics.
      integer, intent(in), optional :: mstep_iter !! Maximum BFGS iterations per generalized-EM M-step; default is 25.
      type(marss_model) :: current
      type(marss_model) :: candidate
      type(marss_kf_result) :: kf
      type(marss_kf_result) :: candidate_kf
      type(marss_hatyt_result) :: hat
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: theta_new(:)
      real(dp), allocatable :: trial_theta(:)
      real(dp) :: alpha
      real(dp) :: old_loglik
      real(dp) :: new_loglik
      integer :: backtrack
      integer :: info
      integer :: iter
      integer :: inner_max
      integer :: p

      fit%converged = .false.
      fit%iterations = 0
      fit%info = 0
      if (start_model%diffuse) then
         fit%info = 5
         return
      end if
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         fit%info = 1
         return
      end if
      inner_max = 25
      if (present(mstep_iter)) inner_max = mstep_iter
      if (inner_max < 1) then
         fit%info = 2
         return
      end if
      if (.not. marss_constraints_valid(start_model, constraints)) then
         fit%info = 3
         return
      end if
      call marss_constraint_start_vector(constraints, theta)
      p = size(theta)
      call marss_apply_constraints(start_model, constraints, theta, current, info)
      if (info /= 0) then
         fit%info = 10 + info
         return
      end if
      call marss_kfss(current, kf)
      if (.not. kf%ok) then
         fit%info = 20 + kf%info
         return
      end if
      old_loglik = kf%loglik
      if (p == 0) then
         fit%model = current
         fit%kf = kf
         fit%loglik = old_loglik
         fit%free_parameters = theta
         fit%converged = .true.
         return
      end if
      allocate(theta_new(p), trial_theta(p))

      do iter = 1, max_iter
         call marss_hatyt(current, kf, hat)
         call maximize_expected_loglik(start_model, constraints, theta, current, kf, hat, inner_max, tol, theta_new, info)
         if (info /= 0) then
            fit%info = 30 + info
            return
         end if

         call marss_apply_constraints(start_model, constraints, theta_new, candidate, info)
         if (info == 0) call marss_kfss(candidate, candidate_kf)
         if (info /= 0 .or. .not. candidate_kf%ok) then
            new_loglik = -huge(1.0_dp)
         else
            new_loglik = candidate_kf%loglik
         end if

         if (new_loglik < old_loglik - 100.0_dp * epsilon(1.0_dp) * (1.0_dp + abs(old_loglik))) then
            alpha = 0.5_dp
            do backtrack = 1, 40
               trial_theta = theta + alpha * (theta_new - theta)
               call marss_apply_constraints(start_model, constraints, trial_theta, candidate, info)
               if (info == 0) then
                  call marss_kfss(candidate, candidate_kf)
                  if (candidate_kf%ok) then
                     new_loglik = candidate_kf%loglik
                     if (new_loglik >= old_loglik - 100.0_dp * epsilon(1.0_dp) * &
                        (1.0_dp + abs(old_loglik))) exit
                  end if
               end if
               alpha = 0.5_dp * alpha
            end do
            if (backtrack > 40) then
               theta_new = theta
               candidate = current
               candidate_kf = kf
               new_loglik = old_loglik
            else
               theta_new = trial_theta
            end if
         end if

         if (abs(new_loglik - old_loglik) <= tol * (1.0_dp + abs(old_loglik)) .and. &
            maxval(abs(theta_new - theta)) <= sqrt(tol) * (1.0_dp + maxval(abs(theta)))) then
            fit%converged = .true.
         end if
         theta = theta_new
         current = candidate
         kf = candidate_kf
         old_loglik = new_loglik
         if (fit%converged) exit
      end do

      fit%iterations = min(iter, max_iter)
      fit%model = current
      fit%kf = kf
      fit%loglik = old_loglik
      fit%free_parameters = theta
      fit%info = 0
   end subroutine marss_kem_linear

   pure subroutine maximize_expected_loglik(template, constraints, theta_start, e_model, kf, hat, &
      max_iter, tol, theta, info)
      type(marss_model), intent(in) :: template !! Template model used to reconstruct affine-constrained M-step candidates.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining the free beta coordinates.
      real(dp), intent(in) :: theta_start(:) !! Beta coordinates at the start of the conditional maximization.
      type(marss_model), intent(in) :: e_model !! Model at which the E-step conditional moments were computed.
      type(marss_kf_result), intent(in) :: kf !! Smoothed state moments held fixed during this M-step.
      type(marss_hatyt_result), intent(in) :: hat !! Conditional observation moments held fixed during this M-step.
      integer, intent(in) :: max_iter !! Maximum BFGS iterations in this conditional maximization.
      real(dp), intent(in) :: tol !! Relative objective and gradient convergence tolerance.
      real(dp), intent(out) :: theta(:) !! Beta coordinates after maximizing the expected complete-data log likelihood.
      integer, intent(out) :: info !! Zero on success; nonzero when all finite-difference or line-search trial points are invalid.
      real(dp), allocatable :: direction(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: gradient_new(:)
      real(dp), allocatable :: h_inv(:, :)
      real(dp), allocatable :: left(:, :)
      real(dp), allocatable :: right(:, :)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: theta_new(:)
      real(dp), allocatable :: y(:)
      real(dp) :: alpha
      real(dp) :: f
      real(dp) :: f_new
      real(dp) :: gtd
      real(dp) :: rho
      real(dp) :: ys
      integer :: iter
      integer :: p

      p = size(theta_start)
      if (size(theta) /= p) then
         info = 1
         return
      end if
      if (p /= marss_constraints_parameter_count(constraints)) then
         info = 2
         return
      end if
      theta = theta_start
      if (p == 0) then
         info = 0
         return
      end if
      allocate(direction(p), gradient(p), gradient_new(p), h_inv(p, p), left(p, p), right(p, p))
      allocate(s(p), theta_new(p), y(p))
      call expected_objective(template, constraints, theta, e_model, kf, hat, f, info)
      if (info /= 0) return
      call expected_gradient(template, constraints, theta, e_model, kf, hat, f, gradient, info)
      if (info /= 0) return
      h_inv = identity_matrix(p)

      do iter = 1, max_iter
         if (maxval(abs(gradient)) <= tol) exit
         direction = -matmul(h_inv, gradient)
         gtd = dot_product(gradient, direction)
         if (gtd >= -epsilon(1.0_dp)) then
            h_inv = identity_matrix(p)
            direction = -gradient
         end if
         call expected_armijo(template, constraints, theta, e_model, kf, hat, f, gradient, direction, &
            theta_new, f_new, alpha, info)
         if (info /= 0) return
         call expected_gradient(template, constraints, theta_new, e_model, kf, hat, f_new, gradient_new, info)
         if (info /= 0) return
         s = theta_new - theta
         y = gradient_new - gradient
         ys = dot_product(y, s)
         if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, sqrt(dot_product(y, y) * dot_product(s, s)))) then
            rho = 1.0_dp / ys
            left = identity_matrix(p) - rho * outer_product(s, y)
            right = identity_matrix(p) - rho * outer_product(y, s)
            h_inv = matmul(matmul(left, h_inv), right) + rho * outer_product(s, s)
         else
            h_inv = identity_matrix(p)
         end if
         theta = theta_new
         gradient = gradient_new
         if (abs(f_new - f) <= tol * (1.0_dp + abs(f))) exit
         f = f_new
      end do
      info = 0
   end subroutine maximize_expected_loglik

   pure subroutine expected_objective(template, constraints, theta, e_model, kf, hat, value, info)
      type(marss_model), intent(in) :: template !! Template model for reconstructing a constrained M-step candidate.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining candidate beta coordinates.
      real(dp), intent(in) :: theta(:) !! Candidate beta coordinates.
      type(marss_model), intent(in) :: e_model !! Model defining the E-step transition convention and data horizon.
      type(marss_kf_result), intent(in) :: kf !! Smoothed moments from the current E-step.
      type(marss_hatyt_result), intent(in) :: hat !! Conditional observation moments from the current E-step.
      real(dp), intent(out) :: value !! Negative expected complete-data log likelihood up to constants independent of the candidate.
      integer, intent(out) :: info !! Zero on success; nonzero if reconstruction or a covariance-support check fails.
      type(marss_model) :: candidate
      real(dp), allocatable :: cross(:, :)
      real(dp), allocatable :: cov_prev(:, :)
      real(dp), allocatable :: mean_prev(:)
      real(dp) :: bt(size(e_model%b, 1), size(e_model%b, 2))
      real(dp) :: qt(size(e_model%q, 1), size(e_model%q, 2))
      real(dp) :: rt(size(e_model%r, 1), size(e_model%r, 2))
      real(dp) :: ut(size(e_model%u))
      real(dp) :: zt(size(e_model%z, 1), size(e_model%z, 2))
      real(dp) :: at(size(e_model%a))
      real(dp), allocatable :: ess(:, :)
      real(dp), allocatable :: ett(:, :)
      real(dp), allocatable :: exs(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: theta_mat(:, :)
      real(dp) :: initial_mean(size(e_model%x0))
      real(dp) :: v0eff(size(e_model%v0, 1), size(e_model%v0, 2))
      real(dp) :: delta(size(e_model%x0))
      real(dp) :: term
      integer :: first_transition
      integer :: m
      integer :: n
      integer :: t
      integer :: tt

      call marss_apply_constraints(template, constraints, theta, candidate, info)
      if (info /= 0) then
         value = huge(1.0_dp)
         return
      end if
      m = size(e_model%b, 1)
      n = size(e_model%z, 1)
      tt = size(e_model%y, 2)
      first_transition = merge(2, 1, e_model%tinitx == 1)
      value = 0.0_dp

      if (e_model%tinitx == 0) then
         initial_mean = kf%x0_smooth
         allocate(residual(m, m))
         delta = initial_mean - candidate%x0
         residual = kf%v0_smooth + outer_product(delta, delta)
      else
         initial_mean = kf%x_smooth(:, 1)
         allocate(residual(m, m))
         delta = initial_mean - candidate%x0
         residual = kf%p_smooth(:, :, 1) + outer_product(delta, delta)
      end if
      v0eff = marss_v0_effective(candidate)
      call covariance_expected_term(v0eff, residual, term, info)
      if (info /= 0) then
         value = huge(1.0_dp)
         return
      end if
      value = value + term
      deallocate(residual)

      do t = first_transition, tt
         call transition_moments(e_model, kf, t, mean_prev, cov_prev, cross)
         allocate(ess(m + 1, m + 1), ett(m, m), exs(m, m + 1), residual(m, m), theta_mat(m, m + 1))
         ess = 0.0_dp
         ess(1:m, 1:m) = cov_prev + outer_product(mean_prev, mean_prev)
         ess(1:m, m + 1) = mean_prev
         ess(m + 1, 1:m) = mean_prev
         ess(m + 1, m + 1) = 1.0_dp
         ett = kf%p_smooth(:, :, t) + outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
         exs(:, 1:m) = cross + outer_product(kf%x_smooth(:, t), mean_prev)
         exs(:, m + 1) = kf%x_smooth(:, t)
         bt = marss_b_at(candidate, t)
         ut = marss_u_at(candidate, t)
         theta_mat(:, 1:m) = bt
         theta_mat(:, m + 1) = ut
         residual = ett - matmul(theta_mat, transpose(exs)) - matmul(exs, transpose(theta_mat)) + &
            matmul(matmul(theta_mat, ess), transpose(theta_mat))
         residual = 0.5_dp * (residual + transpose(residual))
         qt = marss_q_at(candidate, t)
         call covariance_expected_term(qt, residual, term, info)
         if (info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         value = value + term
         deallocate(mean_prev, cov_prev, cross, ess, ett, exs, residual, theta_mat)
      end do

      do t = 1, tt
         allocate(ess(m + 1, m + 1), exs(n, m + 1), residual(n, n), theta_mat(n, m + 1))
         ess = 0.0_dp
         ess(1:m, 1:m) = kf%p_smooth(:, :, t) + outer_product(kf%x_smooth(:, t), kf%x_smooth(:, t))
         ess(1:m, m + 1) = kf%x_smooth(:, t)
         ess(m + 1, 1:m) = kf%x_smooth(:, t)
         ess(m + 1, m + 1) = 1.0_dp
         exs(:, 1:m) = hat%yxt(:, :, t)
         exs(:, m + 1) = hat%yt(:, t)
         zt = marss_z_at(candidate, t)
         at = marss_a_at(candidate, t)
         theta_mat(:, 1:m) = zt
         theta_mat(:, m + 1) = at
         residual = hat%ot(:, :, t) - matmul(theta_mat, transpose(exs)) - matmul(exs, transpose(theta_mat)) + &
            matmul(matmul(theta_mat, ess), transpose(theta_mat))
         residual = 0.5_dp * (residual + transpose(residual))
         rt = marss_r_at(candidate, t)
         call covariance_expected_term(rt, residual, term, info)
         if (info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         value = value + term
         deallocate(ess, exs, residual, theta_mat)
      end do
      info = 0
   end subroutine expected_objective

   pure subroutine covariance_expected_term(covariance, expected_outer, value, info)
      real(dp), intent(in) :: covariance(:, :) !! Candidate covariance matrix for one complete-data Gaussian residual block.
      real(dp), intent(in) :: expected_outer(:, :) !! E[e e'] for that residual under the current E-step distribution.
      real(dp), intent(out) :: value !! Log-pseudodeterminant plus expected quadratic contribution for this block.
      integer, intent(out) :: info !! Zero when covariance is PSD and expected residual mass lies in its numerical support.
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: support_error
      real(dp) :: scale
      real(dp) :: logdet
      integer :: rank

      if (any(shape(covariance) /= shape(expected_outer))) then
         value = huge(1.0_dp)
         info = 1
         return
      end if
      call psd_inverse_logdet(covariance, inverse, logdet, rank, info)
      if (info /= 0) then
         value = huge(1.0_dp)
         return
      end if
      support_error = trace_matrix(expected_outer - matmul(matmul(covariance, inverse), expected_outer))
      scale = max(1.0_dp, abs(trace_matrix(expected_outer)))
      if (support_error > 1.0e-8_dp * scale) then
         value = huge(1.0_dp)
         info = 2
         return
      end if
      value = logdet + trace_matrix(matmul(inverse, expected_outer))
      info = 0
   end subroutine covariance_expected_term

   pure subroutine expected_gradient(template, constraints, theta, e_model, kf, hat, f0, gradient, info)
      type(marss_model), intent(in) :: template !! Template model for constrained expected-likelihood evaluations.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining the beta coordinates.
      real(dp), intent(in) :: theta(:) !! Current M-step beta coordinates.
      type(marss_model), intent(in) :: e_model !! Model used for the current E-step moments.
      type(marss_kf_result), intent(in) :: kf !! Smoothed state moments held fixed during differencing.
      type(marss_hatyt_result), intent(in) :: hat !! Conditional observation moments held fixed during differencing.
      real(dp), intent(in) :: f0 !! Expected objective at theta, used by one-sided fallback differences.
      real(dp), intent(out) :: gradient(:) !! Finite-difference gradient of the negative expected complete-data log likelihood.
      integer, intent(out) :: info !! Zero on success; nonzero when both perturbation directions are invalid for a coordinate.
      real(dp), allocatable :: work(:)
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      integer :: im
      integer :: ip
      integer :: j

      if (size(gradient) /= size(theta)) then
         info = 1
         return
      end if
      allocate(work(size(theta)))
      do j = 1, size(theta)
         h = 1.0e-5_dp * (1.0_dp + abs(theta(j)))
         work = theta
         work(j) = theta(j) + h
         call expected_objective(template, constraints, work, e_model, kf, hat, fp, ip)
         work(j) = theta(j) - h
         call expected_objective(template, constraints, work, e_model, kf, hat, fm, im)
         if (ip == 0 .and. im == 0) then
            gradient(j) = (fp - fm) / (2.0_dp * h)
         else if (ip == 0) then
            gradient(j) = (fp - f0) / h
         else if (im == 0) then
            gradient(j) = (f0 - fm) / h
         else
            info = 10 + j
            return
         end if
      end do
      info = 0
   end subroutine expected_gradient

   pure subroutine expected_armijo(template, constraints, theta, e_model, kf, hat, f, gradient, direction, &
      theta_new, f_new, alpha, info)
      type(marss_model), intent(in) :: template !! Template model for constrained expected-likelihood line-search trials.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining beta coordinates.
      real(dp), intent(in) :: theta(:) !! Current M-step beta coordinates.
      type(marss_model), intent(in) :: e_model !! Model supplying the current E-step moments.
      type(marss_kf_result), intent(in) :: kf !! Smoothed state moments held fixed during the line search.
      type(marss_hatyt_result), intent(in) :: hat !! Conditional observation moments held fixed during the line search.
      real(dp), intent(in) :: f !! Current negative expected complete-data log likelihood.
      real(dp), intent(in) :: gradient(:) !! Current expected-objective gradient.
      real(dp), intent(in) :: direction(:) !! Descent direction proposed by the conditional BFGS update.
      real(dp), intent(out) :: theta_new(:) !! Accepted beta coordinates.
      real(dp), intent(out) :: f_new !! Expected objective at the accepted coordinates.
      real(dp), intent(out) :: alpha !! Accepted line-search step length.
      integer, intent(out) :: info !! Zero on success; nonzero when no covariance-valid Armijo step is found.
      real(dp), parameter :: c1 = 1.0e-4_dp
      real(dp) :: slope
      integer :: eval_info
      integer :: trial

      slope = dot_product(gradient, direction)
      alpha = 1.0_dp
      do trial = 1, 40
         theta_new = theta + alpha * direction
         call expected_objective(template, constraints, theta_new, e_model, kf, hat, f_new, eval_info)
         if (eval_info == 0) then
            if (f_new <= f + c1 * alpha * slope) then
               info = 0
               return
            end if
         end if
         alpha = 0.5_dp * alpha
      end do
      theta_new = theta
      f_new = f
      info = 1
   end subroutine expected_armijo

   pure subroutine transition_moments(model, kf, t, mean_prev, cov_prev, cross)
      type(marss_model), intent(in) :: model !! Model providing the tinitx convention for the predecessor state.
      type(marss_kf_result), intent(in) :: kf !! Smoothed moments from which transition sufficient statistics are extracted.
      integer, intent(in) :: t !! Current state time whose predecessor moment is requested.
      real(dp), allocatable, intent(out) :: mean_prev(:) !! Smoothed predecessor-state mean, including x0 when t=1 and tinitx=0.
      real(dp), allocatable, intent(out) :: cov_prev(:, :) !! Smoothed predecessor-state covariance.
      real(dp), allocatable, intent(out) :: cross(:, :) !! Smoothed covariance Cov[x(t),x(t-1)|Y].
      integer :: m

      m = size(model%b, 1)
      allocate(mean_prev(m), cov_prev(m, m), cross(m, m))
      if (t == 1) then
         mean_prev = kf%x0_smooth
         cov_prev = kf%v0_smooth
      else
         mean_prev = kf%x_smooth(:, t - 1)
         cov_prev = kf%p_smooth(:, :, t - 1)
      end if
      cross = kf%p_lag(:, :, t)
   end subroutine transition_moments

   pure real(dp) function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose trace is returned.
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         value = value + a(i, i)
      end do
   end function trace_matrix

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the rank-one product.
      real(dp), intent(in) :: y(:) !! Right vector in the rank-one product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module marss_constrained_em
