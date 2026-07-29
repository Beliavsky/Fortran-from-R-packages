! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_optimize
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fracdiff_kinds, only : dp
   use fracdiff_status, only : fd_ok, fd_invalid_input, fd_optimization_failed, fd_iteration_limit
   use fracdiff_linalg, only : solve_linear_system, symmetric_rank_k, vector_norm2
   use fracdiff_filter, only : haslett_raftery_filter, arma_residuals, arma_residual_jacobian
   implicit none
   private

   type, public :: profile_evaluation
      integer :: status = fd_ok
      integer :: function_evaluations = 0
      integer :: gradient_evaluations = 0
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: white_noise_variance = huge(1.0_dp)
      real(dp) :: residual_norm = huge(1.0_dp)
      real(dp) :: estimated_mean = 0.0_dp
      real(dp) :: sum_log_v = 0.0_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
   end type profile_evaluation

   public :: levenberg_marquardt_arma
   public :: evaluate_profile_likelihood, evaluate_fixed_likelihood
   public :: brent_profile_d

contains

   subroutine levenberg_marquardt_arma(y, ar, ma, tolerance, maximum_iterations, &
                                       residual_norm, function_evaluations, &
                                       gradient_evaluations, status)
      real(dp), intent(in) :: y(:)
      real(dp), intent(inout) :: ar(:), ma(:)
      real(dp), intent(in) :: tolerance
      integer, intent(in) :: maximum_iterations
      real(dp), intent(out) :: residual_norm
      integer, intent(out) :: function_evaluations, gradient_evaluations, status

      real(dp), allocatable :: residuals(:), trial_residuals(:), jacobian(:,:)
      real(dp), allocatable :: normal(:,:), damped(:,:), gradient(:), step(:)
      real(dp), allocatable :: parameters(:), trial_parameters(:), trial_ar(:), trial_ma(:)
      real(dp) :: sse, trial_sse, lambda, max_diagonal, gradient_scale
      real(dp) :: relative_reduction, step_scale
      integer :: n, p, q, pq, nm, maxpq, iteration, solve_status, filter_status, i

      n = size(y)
      p = size(ar)
      q = size(ma)
      pq = p + q
      maxpq = max(p, q)
      nm = n - maxpq
      function_evaluations = 0
      gradient_evaluations = 0
      status = fd_ok
      residual_norm = huge(1.0_dp)

      if (pq == 0) then
         residual_norm = vector_norm2(y)
         return
      end if
      if (nm <= pq .or. tolerance <= 0.0_dp .or. maximum_iterations < 1) then
         status = fd_invalid_input
         return
      end if

      allocate(residuals(nm), trial_residuals(nm), jacobian(nm,pq))
      allocate(normal(pq,pq), damped(pq,pq), gradient(pq), step(pq))
      allocate(parameters(pq), trial_parameters(pq), trial_ar(p), trial_ma(q))
      if (q > 0) parameters(1:q) = ma
      if (p > 0) parameters(q+1:pq) = ar

      call arma_residual_jacobian(y, ar, ma, residuals, jacobian, filter_status)
      if (filter_status /= fd_ok) then
         status = filter_status
         return
      end if
      function_evaluations = 1
      gradient_evaluations = 1
      sse = dot_product(residuals, residuals)
      call symmetric_rank_k(jacobian, normal)
      gradient = matmul(transpose(jacobian), residuals)
      max_diagonal = max(1.0_dp, maxval([(abs(normal(i,i)), i=1,pq)]))
      lambda = 1.0e-3_dp*max_diagonal

      do iteration = 1, maximum_iterations
         gradient_scale = maxval(abs(gradient))
         if (gradient_scale <= tolerance*(1.0_dp + sqrt(max(sse,0.0_dp)))) exit

         damped = normal
         do i = 1, pq
            damped(i,i) = damped(i,i) + lambda*max(1.0_dp, normal(i,i))
         end do
         call solve_linear_system(damped, -gradient, step, solve_status)
         if (solve_status /= 0) then
            lambda = lambda*10.0_dp
            if (lambda > 1.0e20_dp) then
               status = fd_optimization_failed
               exit
            end if
            cycle
         end if

         trial_parameters = parameters + step
         if (maxval(abs(trial_parameters)) > 100.0_dp) then
            lambda = lambda*10.0_dp
            cycle
         end if
         if (q > 0) trial_ma = trial_parameters(1:q)
         if (p > 0) trial_ar = trial_parameters(q+1:pq)
         call arma_residuals(y, trial_ar, trial_ma, trial_residuals, filter_status)
         function_evaluations = function_evaluations + 1
         if (filter_status /= fd_ok) then
            lambda = lambda*10.0_dp
            cycle
         end if
         trial_sse = dot_product(trial_residuals, trial_residuals)

         if (ieee_is_finite(trial_sse) .and. trial_sse < sse) then
            relative_reduction = (sse - trial_sse)/max(1.0_dp, sse)
            step_scale = vector_norm2(step)/max(1.0_dp, vector_norm2(parameters))
            parameters = trial_parameters
            residuals = trial_residuals
            sse = trial_sse
            if (q > 0) ma = parameters(1:q)
            if (p > 0) ar = parameters(q+1:pq)
            call arma_residual_jacobian(y, ar, ma, residuals, jacobian, filter_status)
            function_evaluations = function_evaluations + 1
            gradient_evaluations = gradient_evaluations + 1
            if (filter_status /= fd_ok) then
               status = filter_status
               exit
            end if
            call symmetric_rank_k(jacobian, normal)
            gradient = matmul(transpose(jacobian), residuals)
            lambda = max(lambda/3.0_dp, 1.0e-15_dp)
            if (relative_reduction <= tolerance .and. step_scale <= sqrt(tolerance)) exit
         else
            lambda = lambda*10.0_dp
            if (lambda > 1.0e20_dp) then
               status = fd_optimization_failed
               exit
            end if
         end if
      end do

      if (iteration > maximum_iterations .and. status == fd_ok) status = fd_iteration_limit
      residual_norm = sqrt(max(sse, 0.0_dp))
   end subroutine levenberg_marquardt_arma

   subroutine evaluate_profile_likelihood(x, d, m_terms, ar_start, ma_start, tolerance, evaluation)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: d
      integer, intent(in) :: m_terms
      real(dp), intent(in) :: ar_start(:), ma_start(:)
      real(dp), intent(in) :: tolerance
      type(profile_evaluation), intent(out) :: evaluation

      real(dp), allocatable :: filtered(:)
      real(dp) :: n_real, residual_norm
      integer :: n, p, q, maxpq, nm, filter_status, optimizer_status
      integer :: nfun, ngrad

      n = size(x)
      p = size(ar_start)
      q = size(ma_start)
      maxpq = max(p,q)
      nm = n - maxpq
      evaluation%status = fd_ok
      allocate(evaluation%ar(p), evaluation%ma(q), filtered(n))
      evaluation%ar = ar_start
      evaluation%ma = ma_start

      call haslett_raftery_filter(x, d, m_terms, p + q > 0, filtered, &
                                  evaluation%sum_log_v, evaluation%estimated_mean, filter_status)
      if (filter_status /= fd_ok) then
         evaluation%status = filter_status
         return
      end if

      n_real = real(n, dp)
      if (p + q == 0) then
         evaluation%residual_norm = sqrt(dot_product(filtered, filtered))
         evaluation%white_noise_variance = dot_product(filtered, filtered)/n_real
      else
         call levenberg_marquardt_arma(filtered, evaluation%ar, evaluation%ma, tolerance, 100, &
                                      residual_norm, nfun, ngrad, optimizer_status)
         evaluation%function_evaluations = nfun
         evaluation%gradient_evaluations = ngrad
         if (optimizer_status /= fd_ok .and. optimizer_status /= fd_iteration_limit) then
            evaluation%status = optimizer_status
            return
         end if
         evaluation%residual_norm = residual_norm
         evaluation%white_noise_variance = residual_norm**2/real(nm - 1, dp)
         if (optimizer_status == fd_iteration_limit) evaluation%status = optimizer_status
      end if

      if (evaluation%white_noise_variance <= 0.0_dp .or. &
          .not. ieee_is_finite(evaluation%white_noise_variance)) then
         evaluation%status = fd_optimization_failed
         return
      end if
      evaluation%objective = 0.5_dp*(n_real*(log(evaluation%white_noise_variance) + 2.8378_dp) + &
                                      evaluation%sum_log_v)
      evaluation%log_likelihood = -evaluation%objective
   end subroutine evaluate_profile_likelihood

   subroutine evaluate_fixed_likelihood(x, d, m_terms, ar, ma, log_likelihood, &
                                        white_noise_variance, estimated_mean, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: d
      integer, intent(in) :: m_terms
      real(dp), intent(in) :: ar(:), ma(:)
      real(dp), intent(out) :: log_likelihood, white_noise_variance, estimated_mean
      integer, intent(out) :: status

      real(dp), allocatable :: filtered(:), residuals(:)
      real(dp) :: sum_log_v
      integer :: n, p, q, maxpq, nm

      n = size(x)
      p = size(ar)
      q = size(ma)
      maxpq = max(p,q)
      nm = n - maxpq
      allocate(filtered(n))
      call haslett_raftery_filter(x, d, m_terms, p + q > 0, filtered, sum_log_v, estimated_mean, status)
      if (status /= fd_ok) then
         log_likelihood = -huge(1.0_dp)
         white_noise_variance = huge(1.0_dp)
         return
      end if

      if (p + q == 0) then
         white_noise_variance = dot_product(filtered, filtered)/real(n,dp)
      else
         allocate(residuals(nm))
         call arma_residuals(filtered, ar, ma, residuals, status)
         if (status /= fd_ok) then
            log_likelihood = -huge(1.0_dp)
            white_noise_variance = huge(1.0_dp)
            return
         end if
         white_noise_variance = dot_product(residuals,residuals)/real(nm - 1,dp)
      end if
      if (white_noise_variance <= 0.0_dp .or. .not. ieee_is_finite(white_noise_variance)) then
         status = fd_optimization_failed
         log_likelihood = -huge(1.0_dp)
         return
      end if
      log_likelihood = -0.5_dp*(real(n,dp)*(log(white_noise_variance) + 2.8378_dp) + sum_log_v)
   end subroutine evaluate_fixed_likelihood

   subroutine brent_profile_d(x, m_terms, d_initial, d_range, d_tolerance, ar, ma, &
                              d_optimum, evaluation, outer_iterations, trace)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: m_terms
      real(dp), intent(in) :: d_initial, d_range(2), d_tolerance
      real(dp), intent(inout) :: ar(:), ma(:)
      real(dp), intent(out) :: d_optimum
      type(profile_evaluation), intent(out) :: evaluation
      integer, intent(out) :: outer_iterations
      logical, intent(in), optional :: trace

      real(dp), parameter :: golden = 0.38196601125011_dp
      real(dp) :: a, b, dstep, estep, midpoint, fu, fv, fw, fx
      real(dp) :: pfit, qfit, rfit, tol1, tol2, u, v, w, xpoint
      real(dp) :: sqrt_eps
      real(dp), allocatable :: current_ar(:), current_ma(:)
      type(profile_evaluation) :: trial
      logical :: do_trace
      integer :: iteration

      do_trace = .false.
      if (present(trace)) do_trace = trace
      a = d_range(1)
      b = d_range(2)
      sqrt_eps = sqrt(epsilon(1.0_dp))
      if (d_initial > a + d_tolerance .and. d_initial < b - d_tolerance) then
         xpoint = d_initial
      else
         xpoint = a + golden*(b - a)
      end if
      w = xpoint
      v = xpoint
      dstep = 0.0_dp
      estep = 0.0_dp
      allocate(current_ar(size(ar)), current_ma(size(ma)))
      current_ar = ar
      current_ma = ma
      call evaluate_profile_likelihood(x, xpoint, m_terms, current_ar, current_ma, &
                                       max(d_tolerance/10.0_dp, epsilon(1.0_dp)**0.3_dp), evaluation)
      if (allocated(evaluation%ar)) current_ar = evaluation%ar
      if (allocated(evaluation%ma)) current_ma = evaluation%ma
      fx = evaluation%objective
      fw = fx
      fv = fx

      if (do_trace) write(*,'(a)') " iteration              d        objective"
      if (do_trace) write(*,'(i10,2es16.7)') 1, xpoint, fx

      do iteration = 2, 100
         midpoint = 0.5_dp*(a + b)
         tol1 = sqrt_eps*(abs(xpoint) + 1.0_dp) + d_tolerance/3.0_dp
         tol2 = 2.0_dp*tol1
         if (abs(xpoint - midpoint) + 0.5_dp*(b - a) <= tol2) exit

         pfit = 0.0_dp
         qfit = 0.0_dp
         rfit = 0.0_dp
         if (abs(estep) > tol1) then
            rfit = (xpoint - w)*(fx - fv)
            qfit = (xpoint - v)*(fx - fw)
            pfit = (xpoint - v)*qfit - (xpoint - w)*rfit
            qfit = 2.0_dp*(qfit - rfit)
            if (qfit > 0.0_dp) pfit = -pfit
            qfit = abs(qfit)
            rfit = estep
            estep = dstep
         end if

         if (abs(pfit) < abs(0.5_dp*qfit*rfit) .and. pfit > qfit*(a - xpoint) .and. &
             pfit < qfit*(b - xpoint)) then
            dstep = pfit/qfit
            u = xpoint + dstep
            if (u - a < tol2 .or. b - u < tol2) dstep = merge(tol1, -tol1, xpoint < midpoint)
         else
            estep = merge(a - xpoint, b - xpoint, xpoint >= midpoint)
            dstep = golden*estep
         end if

         if (abs(dstep) >= tol1) then
            u = xpoint + dstep
         else
            u = xpoint + merge(tol1, -tol1, dstep > 0.0_dp)
         end if

         call evaluate_profile_likelihood(x, u, m_terms, current_ar, current_ma, &
                                          max(d_tolerance/10.0_dp, epsilon(1.0_dp)**0.3_dp), trial)
         if (allocated(trial%ar)) current_ar = trial%ar
         if (allocated(trial%ma)) current_ma = trial%ma
         fu = trial%objective
         if (do_trace) write(*,'(i10,2es16.7)') iteration, u, fu

         if (fu <= fx) then
            if (u >= xpoint) then
               a = xpoint
            else
               b = xpoint
            end if
            v = w
            fv = fw
            w = xpoint
            fw = fx
            xpoint = u
            fx = fu
            evaluation = trial
         else
            if (u >= xpoint) then
               b = u
            else
               a = u
            end if
            if (fu <= fw .or. nearly_equal(w, xpoint)) then
               v = w
               fv = fw
               w = u
               fw = fu
            else if (fu <= fv .or. nearly_equal(v, xpoint) .or. nearly_equal(v, w)) then
               v = u
               fv = fu
            end if
         end if
      end do

      outer_iterations = min(iteration, 100)
      d_optimum = xpoint
      call evaluate_profile_likelihood(x, d_optimum, m_terms, evaluation%ar, evaluation%ma, &
                                       max(d_tolerance/10.0_dp, epsilon(1.0_dp)**0.3_dp), trial)
      evaluation = trial
      ar = evaluation%ar
      ma = evaluation%ma
      if (iteration > 100 .and. evaluation%status == fd_ok) evaluation%status = fd_iteration_limit
   end subroutine brent_profile_d

   pure function nearly_equal(a, b) result(equal)
      real(dp), intent(in) :: a, b
      logical :: equal

      equal = abs(a - b) <= epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))
   end function nearly_equal

end module fracdiff_optimize
