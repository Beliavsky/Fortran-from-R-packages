! Profile-likelihood and confidence-interval kernels based on ordinal/R/clm.profile.R.
! Copyright (C) 2010-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_profile
   use ordinal_kinds, only : dp
   use ordinal_numerics, only : cholesky_factor, solve_cholesky
   use ordinal_clm, only : clm_problem, clm_fit_result, clm_nll, clm_analytic_gradient, clm_analytic_hessian
   implicit none
   private

   public :: clm_profile_likelihood
   public :: clm_profile_confidence_interval
   public :: clm_wald_confidence_interval
contains
   subroutine pack_clm_fit(problem, fit, par, status)
      type(clm_problem), intent(in) :: problem !! Model specification used to interpret the fitted coefficient blocks.
      type(clm_fit_result), intent(in) :: fit !! Fitted CLM supplying threshold, location, and scale coefficients.
      real(dp), allocatable, intent(out) :: par(:) !! Packed fitted coefficients in the objective's parameter order.
      integer, intent(out) :: status !! Zero on success; nonzero when flexible-link profiling or dimensions are unsupported.
      integer :: nalpha, nbeta, nzeta, npar

      status = 0
      if (problem%estimate_lambda) then
         status = 1
         allocate(par(0))
         return
      end if
      nalpha = size(fit%alpha)
      nbeta = size(fit%beta)
      nzeta = size(fit%zeta)
      npar = nalpha + nbeta + nzeta
      allocate(par(npar))
      if (nalpha > 0) par(:nalpha) = fit%alpha
      if (nbeta > 0) par(nalpha + 1:nalpha + nbeta) = fit%beta
      if (nzeta > 0) par(nalpha + nbeta + 1:) = fit%zeta
      if (size(fit%vcov, 1) /= npar .or. size(fit%vcov, 2) /= npar) status = 2
   end subroutine pack_clm_fit

   subroutine profile_one(problem, mle_par, parameter_index, fixed_value, loglik, status, max_iter, grad_tol)
      type(clm_problem), intent(in) :: problem !! Fixed-link CLM whose likelihood is profiled.
      real(dp), intent(in) :: mle_par(:) !! Full packed MLE used to initialize all nuisance coefficients.
      integer, intent(in) :: parameter_index !! One-based packed coefficient index to hold fixed.
      real(dp), intent(in) :: fixed_value !! Value at which the selected coefficient is constrained.
      real(dp), intent(out) :: loglik !! Maximized log-likelihood subject to the fixed coefficient value.
      integer, intent(out) :: status !! Zero on successful constrained Newton optimization; nonzero otherwise.
      integer, intent(in) :: max_iter !! Maximum constrained Newton iterations.
      real(dp), intent(in) :: grad_tol !! Infinity-norm convergence tolerance for nuisance coefficients.
      real(dp), allocatable :: free_par(:), trial_free(:), full_par(:), gradient(:), free_gradient(:)
      real(dp), allocatable :: hessian(:, :), free_hessian(:, :), work_hessian(:, :), chol(:, :), step(:)
      real(dp) :: fval, trial_f, alpha_step, inflation, scale
      integer :: npar, nfree, iter, ls, st, chst, solst, i

      npar = size(mle_par)
      loglik = -huge(1.0_dp)
      if (parameter_index < 1 .or. parameter_index > npar) then
         status = 10
         return
      end if
      allocate(full_par(npar), gradient(npar), hessian(npar, npar))
      if (npar == 1) then
         full_par(1) = fixed_value
         fval = clm_nll(full_par, problem)
         if (fval >= huge(1.0_dp)/100.0_dp) then
            status = 11
            return
         end if
         loglik = -fval
         status = 0
         return
      end if

      nfree = npar - 1
      allocate(free_par(nfree), trial_free(nfree), free_gradient(nfree), free_hessian(nfree, nfree), &
               work_hessian(nfree, nfree), chol(nfree, nfree), step(nfree))
      call remove_vector_entry(mle_par, parameter_index, free_par)
      status = 1
      do iter = 1, max_iter
         call insert_vector_entry(free_par, parameter_index, fixed_value, full_par)
         fval = clm_nll(full_par, problem)
         if (fval >= huge(1.0_dp)/100.0_dp) then
            status = 11
            return
         end if
         call clm_analytic_gradient(full_par, problem, gradient, st)
         if (st /= 0) then
            status = 12
            return
         end if
         call clm_analytic_hessian(full_par, problem, hessian, st)
         if (st /= 0) then
            status = 13
            return
         end if
         call remove_vector_entry(gradient, parameter_index, free_gradient)
         call remove_matrix_row_column(hessian, parameter_index, free_hessian)
         if (maxval(abs(free_gradient)) <= grad_tol) then
            status = 0
            exit
         end if

         scale = max(1.0_dp, maxval(abs(free_hessian)))
         inflation = 0.0_dp
         do i = 1, 14
            work_hessian = free_hessian
            if (inflation > 0.0_dp) then
               work_hessian = work_hessian + diagonal_matrix(nfree, inflation)
            end if
            call cholesky_factor(work_hessian, chol, chst)
            if (chst == 0) exit
            if (inflation <= 0.0_dp) then
               inflation = 1.0e-8_dp*scale
            else
               inflation = 10.0_dp*inflation
            end if
         end do
         if (chst /= 0) then
            status = 14
            return
         end if
         call solve_cholesky(chol, free_gradient, step, solst)
         if (solst /= 0) then
            status = 15
            return
         end if
         step = -step
         alpha_step = 1.0_dp
         do ls = 1, 50
            trial_free = free_par + alpha_step*step
            call insert_vector_entry(trial_free, parameter_index, fixed_value, full_par)
            trial_f = clm_nll(full_par, problem)
            if (trial_f < huge(1.0_dp)/100.0_dp .and. trial_f < fval) exit
            alpha_step = 0.5_dp*alpha_step
         end do
         if (ls > 50) then
            if (maxval(abs(alpha_step*step)) <= 1.0e-10_dp*max(1.0_dp, maxval(abs(free_par)))) then
               status = 0
               exit
            end if
            status = 16
            return
         end if
         free_par = trial_free
      end do
      call insert_vector_entry(free_par, parameter_index, fixed_value, full_par)
      fval = clm_nll(full_par, problem)
      if (fval >= huge(1.0_dp)/100.0_dp) then
         status = 17
         return
      end if
      if (status /= 0) then
         call clm_analytic_gradient(full_par, problem, gradient, st)
         if (st == 0) then
            call remove_vector_entry(gradient, parameter_index, free_gradient)
            if (maxval(abs(free_gradient)) <= 10.0_dp*grad_tol) status = 0
         end if
      end if
      loglik = -fval
   end subroutine profile_one

   pure subroutine remove_vector_entry(full, index, reduced)
      real(dp), intent(in) :: full(:) !! Full vector from which one indexed entry is removed.
      integer, intent(in) :: index !! One-based entry to omit from the full vector.
      real(dp), intent(out) :: reduced(:) !! Vector containing all full entries except the omitted one.

      if (index > 1) reduced(:index - 1) = full(:index - 1)
      if (index < size(full)) reduced(index:) = full(index + 1:)
   end subroutine remove_vector_entry

   pure subroutine insert_vector_entry(reduced, index, value, full)
      real(dp), intent(in) :: reduced(:) !! Vector missing one fixed coefficient.
      integer, intent(in) :: index !! One-based position at which the fixed coefficient is inserted.
      real(dp), intent(in) :: value !! Fixed coefficient value inserted into the full vector.
      real(dp), intent(out) :: full(:) !! Reconstructed full coefficient vector.

      if (index > 1) full(:index - 1) = reduced(:index - 1)
      full(index) = value
      if (index < size(full)) full(index + 1:) = reduced(index:)
   end subroutine insert_vector_entry

   pure subroutine remove_matrix_row_column(full, index, reduced)
      real(dp), intent(in) :: full(:, :) !! Full square matrix from which one matching row and column are removed.
      integer, intent(in) :: index !! One-based row and column to omit.
      real(dp), intent(out) :: reduced(:, :) !! Principal submatrix after removing the selected row and column.
      integer :: i, j, ii, jj

      ii = 0
      do i = 1, size(full, 1)
         if (i == index) cycle
         ii = ii + 1
         jj = 0
         do j = 1, size(full, 2)
            if (j == index) cycle
            jj = jj + 1
            reduced(ii, jj) = full(i, j)
         end do
      end do
   end subroutine remove_matrix_row_column

   pure function diagonal_matrix(n, value) result(matrix)
      integer, intent(in) :: n !! Order of the square diagonal matrix.
      real(dp), intent(in) :: value !! Constant placed on every diagonal entry.
      real(dp) :: matrix(n, n)
      integer :: i

      matrix = 0.0_dp
      do i = 1, n
         matrix(i, i) = value
      end do
   end function diagonal_matrix

   subroutine clm_profile_likelihood(problem, fit, parameter_index, values, loglik, status, max_iter, grad_tol)
      type(clm_problem), intent(in) :: problem !! Fixed-link CLM to profile; estimated flexible-link shapes are unsupported.
      type(clm_fit_result), intent(in) :: fit !! Converged CLM fit providing the MLE and nuisance-parameter starts.
      integer, intent(in) :: parameter_index !! One-based index in packed alpha, beta, then zeta coefficient order.
      real(dp), intent(in) :: values(:) !! Fixed coefficient values at which the profile likelihood is requested.
      real(dp), intent(out) :: loglik(:) !! Profiled log-likelihood values corresponding one-for-one with values.
      integer, intent(out) :: status(:) !! Per-value status: zero on successful nuisance optimization.
      integer, intent(in), optional :: max_iter !! Maximum nuisance BFGS iterations per profile point; defaults to 200.
      real(dp), intent(in), optional :: grad_tol !! Nuisance optimizer gradient tolerance; defaults to 1e-7.
      real(dp), allocatable :: mle_par(:)
      integer :: i, limit, pack_status
      real(dp) :: tol

      if (size(loglik) /= size(values) .or. size(status) /= size(values)) then
         if (size(status) > 0) status = 20
         if (size(loglik) > 0) loglik = -huge(1.0_dp)
         return
      end if
      call pack_clm_fit(problem, fit, mle_par, pack_status)
      if (pack_status /= 0 .or. parameter_index < 1 .or. parameter_index > size(mle_par)) then
         status = 21 + pack_status
         loglik = -huge(1.0_dp)
         return
      end if
      limit = 200
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-7_dp
      if (present(grad_tol)) tol = grad_tol
      do i = 1, size(values)
         call profile_one(problem, mle_par, parameter_index, values(i), loglik(i), status(i), limit, tol)
      end do
   end subroutine clm_profile_likelihood

   subroutine clm_wald_confidence_interval(fit, parameter_index, level, lower, upper, status)
      type(clm_fit_result), intent(in) :: fit !! Fitted CLM with an invertible Hessian covariance estimate.
      integer, intent(in) :: parameter_index !! One-based packed alpha, beta, then zeta coefficient index.
      real(dp), intent(in) :: level !! Two-sided confidence level strictly between zero and one.
      real(dp), intent(out) :: lower !! Lower endpoint of the normal-theory Wald interval.
      real(dp), intent(out) :: upper !! Upper endpoint of the normal-theory Wald interval.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid level, index, or variance.
      real(dp), allocatable :: par(:)
      real(dp) :: z, se
      integer :: nalpha, nbeta, nzeta, npar

      lower = 0.0_dp
      upper = 0.0_dp
      if (level <= 0.0_dp .or. level >= 1.0_dp) then
         status = 1
         return
      end if
      nalpha = size(fit%alpha)
      nbeta = size(fit%beta)
      nzeta = size(fit%zeta)
      npar = nalpha + nbeta + nzeta
      if (parameter_index < 1 .or. parameter_index > npar .or. &
          size(fit%vcov, 1) /= npar .or. size(fit%vcov, 2) /= npar) then
         status = 2
         return
      end if
      if (fit%vcov(parameter_index, parameter_index) <= 0.0_dp) then
         status = 3
         return
      end if
      allocate(par(npar))
      if (nalpha > 0) par(:nalpha) = fit%alpha
      if (nbeta > 0) par(nalpha + 1:nalpha + nbeta) = fit%beta
      if (nzeta > 0) par(nalpha + nbeta + 1:) = fit%zeta
      z = normal_quantile(0.5_dp*(1.0_dp + level))
      se = sqrt(fit%vcov(parameter_index, parameter_index))
      lower = par(parameter_index) - z*se
      upper = par(parameter_index) + z*se
      status = 0
   end subroutine clm_wald_confidence_interval

   subroutine clm_profile_confidence_interval(problem, fit, parameter_index, level, lower, upper, status, max_iter, grad_tol)
      type(clm_problem), intent(in) :: problem !! Fixed-link CLM whose coefficient profile defines the interval.
      type(clm_fit_result), intent(in) :: fit !! Converged CLM fit supplying the MLE, covariance, and reference log-likelihood.
      integer, intent(in) :: parameter_index !! One-based packed alpha, beta, then zeta coefficient index to profile.
      real(dp), intent(in) :: level !! Two-sided profile-likelihood confidence level strictly between zero and one.
      real(dp), intent(out) :: lower !! Lower fixed-coefficient value where the likelihood-root reaches the requested cutoff.
      real(dp), intent(out) :: upper !! Upper fixed-coefficient value where the likelihood-root reaches the requested cutoff.
      integer, intent(out) :: status !! Zero on success; nonzero for unsupported input, failed bracketing, or optimization.
      integer, intent(in), optional :: max_iter !! Maximum nuisance BFGS iterations per profile evaluation; defaults to 200.
      real(dp), intent(in), optional :: grad_tol !! Nuisance optimizer gradient tolerance; defaults to 1e-7.
      real(dp), allocatable :: mle_par(:)
      real(dp) :: target, z, se, initial_step, tol, endpoint
      integer :: pack_status, limit, side_status

      lower = 0.0_dp
      upper = 0.0_dp
      if (level <= 0.0_dp .or. level >= 1.0_dp) then
         status = 1
         return
      end if
      call pack_clm_fit(problem, fit, mle_par, pack_status)
      if (pack_status /= 0 .or. parameter_index < 1 .or. parameter_index > size(mle_par)) then
         status = 2 + pack_status
         return
      end if
      if (.not. fit%converged) then
         status = 5
         return
      end if
      z = normal_quantile(0.5_dp*(1.0_dp + level))
      target = fit%loglik - 0.5_dp*z*z
      if (fit%vcov(parameter_index, parameter_index) > 0.0_dp) then
         se = sqrt(fit%vcov(parameter_index, parameter_index))
      else
         se = 0.0_dp
      end if
      initial_step = max(0.25_dp*max(se, 0.0_dp), 0.025_dp*max(1.0_dp, abs(mle_par(parameter_index))))
      limit = 200
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-7_dp
      if (present(grad_tol)) tol = grad_tol

      call profile_endpoint(problem, mle_par, fit%loglik, parameter_index, -1, target, initial_step, endpoint, &
                            side_status, limit, tol)
      if (side_status /= 0) then
         status = 10 + side_status
         return
      end if
      lower = endpoint
      call profile_endpoint(problem, mle_par, fit%loglik, parameter_index, 1, target, initial_step, endpoint, &
                            side_status, limit, tol)
      if (side_status /= 0) then
         status = 20 + side_status
         return
      end if
      upper = endpoint
      status = 0
   end subroutine clm_profile_confidence_interval

   subroutine profile_endpoint(problem, mle_par, mle_loglik, parameter_index, direction, target, initial_step, endpoint, &
                               status, max_iter, grad_tol)
      type(clm_problem), intent(in) :: problem !! Fixed-link CLM whose nuisance coefficients are reoptimized at each trial point.
      real(dp), intent(in) :: mle_par(:) !! Packed unconstrained MLE defining the starting coefficient value.
      real(dp), intent(in) :: mle_loglik !! Unconstrained maximized log-likelihood used as the profile peak.
      integer, intent(in) :: parameter_index !! One-based packed coefficient index being profiled.
      integer, intent(in) :: direction !! Search direction, -1 for the lower endpoint and +1 for the upper endpoint.
      real(dp), intent(in) :: target !! Profile log-likelihood cutoff defining the confidence endpoint.
      real(dp), intent(in) :: initial_step !! Positive initial displacement from the MLE used for bracketing.
      real(dp), intent(out) :: endpoint !! Fixed coefficient value whose profile log-likelihood is approximately target.
      integer, intent(out) :: status !! Zero on success; nonzero when the endpoint cannot be bracketed or evaluated.
      integer, intent(in) :: max_iter !! Maximum nuisance BFGS iterations at each trial coefficient value.
      real(dp), intent(in) :: grad_tol !! Nuisance optimizer gradient tolerance at each trial coefficient value.
      real(dp) :: center, inside, outside, trial, loglik, step
      integer :: i, profile_status

      endpoint = mle_par(parameter_index)
      if (direction /= -1 .and. direction /= 1) then
         status = 1
         return
      end if
      center = mle_par(parameter_index)
      inside = center
      step = max(initial_step, sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(center)))
      outside = center + real(direction, dp)*step
      do i = 1, 32
         call profile_one(problem, mle_par, parameter_index, outside, loglik, profile_status, max_iter, grad_tol)
         if (profile_status /= 0 .and. profile_status /= 1 .and. profile_status /= 2) then
            loglik = -huge(1.0_dp)
         end if
         loglik = min(loglik, mle_loglik)
         if (loglik <= target) exit
         inside = outside
         step = 1.75_dp*step
         outside = center + real(direction, dp)*step
      end do
      if (i > 32) then
         status = 2
         return
      end if
      do i = 1, 36
         trial = 0.5_dp*(inside + outside)
         call profile_one(problem, mle_par, parameter_index, trial, loglik, profile_status, max_iter, grad_tol)
         if (profile_status /= 0 .and. profile_status /= 1 .and. profile_status /= 2) then
            loglik = -huge(1.0_dp)
         end if
         loglik = min(loglik, mle_loglik)
         if (loglik > target) then
            inside = trial
         else
            outside = trial
         end if
         if (abs(outside - inside) <= 1.0e-7_dp*max(1.0_dp, abs(trial))) exit
      end do
      endpoint = 0.5_dp*(inside + outside)
      status = 0
   end subroutine profile_endpoint

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p !! Standard-normal lower-tail probability strictly between zero and one.
      real(dp), parameter :: a1 = -3.969683028665376e1_dp
      real(dp), parameter :: a2 = 2.209460984245205e2_dp
      real(dp), parameter :: a3 = -2.759285104469687e2_dp
      real(dp), parameter :: a4 = 1.383577518672690e2_dp
      real(dp), parameter :: a5 = -3.066479806614716e1_dp
      real(dp), parameter :: a6 = 2.506628277459239_dp
      real(dp), parameter :: b1 = -5.447609879822406e1_dp
      real(dp), parameter :: b2 = 1.615858368580409e2_dp
      real(dp), parameter :: b3 = -1.556989798598866e2_dp
      real(dp), parameter :: b4 = 6.680131188771972e1_dp
      real(dp), parameter :: b5 = -1.328068155288572e1_dp
      real(dp), parameter :: c1 = -7.784894002430293e-3_dp
      real(dp), parameter :: c2 = -3.223964580411365e-1_dp
      real(dp), parameter :: c3 = -2.400758277161838_dp
      real(dp), parameter :: c4 = -2.549732539343734_dp
      real(dp), parameter :: c5 = 4.374664141464968_dp
      real(dp), parameter :: c6 = 2.938163982698783_dp
      real(dp), parameter :: d1 = 7.784695709041462e-3_dp
      real(dp), parameter :: d2 = 3.224671290700398e-1_dp
      real(dp), parameter :: d3 = 2.445134137142996_dp
      real(dp), parameter :: d4 = 3.754408661907416_dp
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q, r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6)/((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      else if (p > phigh) then
         q = sqrt(-2.0_dp*log(1.0_dp - p))
         x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6)/((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
      else
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q/(((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_dp)
      end if
   end function normal_quantile
end module ordinal_profile
