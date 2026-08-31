! Cumulative-link model kernels based on ordinal/R/clm*.R and src/get_fitted.c.
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_clm
   use ordinal_kinds, only : dp
   use ordinal_links, only : link_cdf, link_pdf, link_pdf_gradient, link_logit, link_aranda_ordaz, link_log_gamma
   use ordinal_thresholds, only : threshold_flexible, threshold_parameter_count, threshold_start, thresholds_from_alpha, &
                                  threshold_jacobian
   use ordinal_numerics, only : objective_type, bfgs_minimize, numerical_gradient, numerical_hessian, invert_matrix, &
                                cholesky_factor, solve_cholesky, symmetric_eigenvalues, hessian_diagnostics
   implicit none
   private
   real(dp), parameter :: bad_objective = huge(1.0_dp)/1000.0_dp
   type, extends(objective_type), public :: clm_problem
      integer, allocatable :: y(:)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: nominal_x(:, :)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: offset(:)
      real(dp), allocatable :: scale_x(:, :)
      real(dp), allocatable :: scale_offset(:)
      integer :: nclass = 0
      integer :: link = link_logit
      integer :: threshold = threshold_flexible
      real(dp) :: lambda = 1.0_dp
      logical :: estimate_lambda = .false.
   contains
      procedure :: value => clm_value
   end type clm_problem
   type, public :: clm_fit_result
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: nominal_alpha(:, :)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: zeta(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: vcov(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp) :: lambda = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: max_gradient = huge(1.0_dp)
      real(dp) :: hessian_min_eigenvalue = 0.0_dp
      real(dp) :: hessian_max_eigenvalue = 0.0_dp
      real(dp) :: hessian_condition = huge(1.0_dp)
      integer :: hessian_rank = 0
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      logical :: hessian_positive_definite = .false.
   end type clm_fit_result
   public :: init_clm_problem, clm_nll, fit_clm, clm_predict_proba, clm_predict_proba_nominal
   public :: clm_analytic_gradient, clm_analytic_hessian, clm_nominal_thresholds
contains
   subroutine init_clm_problem(problem, y, x, nclass, link, threshold, weights, offset, scale_x, scale_offset, lambda, &
                               estimate_lambda, nominal_x, status)
      type(clm_problem), intent(out) :: problem !! Initialized cumulative-link problem returned to the caller.
      integer, intent(in) :: y(:) !! Ordered response codes in 1:nclass, one per row of x.
      real(dp), intent(in) :: x(:, :) !! Location design matrix without an intercept; rows are observations.
      integer, intent(in) :: nclass !! Number of ordered response categories; must be at least two.
      integer, intent(in), optional :: link !! Link identifier; defaults to logistic.
      integer, intent(in), optional :: threshold !! Threshold structure identifier; defaults to flexible.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights; defaults to one.
      real(dp), intent(in), optional :: offset(:) !! Location offsets, one per observation; defaults to zero.
      real(dp), intent(in), optional :: scale_x(:, :) !! Log-scale design matrix without an intercept; rows match y.
      real(dp), intent(in), optional :: scale_offset(:) !! Fixed log-scale offsets, one per observation; defaults to zero.
      real(dp), intent(in), optional :: lambda !! Fixed flexible-link shape, or starting value when it is estimated.
      logical, intent(in), optional :: estimate_lambda !! Estimate the Aranda-Ordaz or log-gamma shape when true.
      real(dp), intent(in), optional :: nominal_x(:, :) !! Nominal predictors without an intercept; each shifts all cut points.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, codes, weights, or options.
      integer :: n
      n = size(y)
      status = 0
      if (nclass < 2 .or. size(x, 1) /= n .or. any(y < 1) .or. any(y > nclass)) then
         status = 1
         return
      end if
      problem%nclass = nclass
      if (present(link)) problem%link = link
      if (present(threshold)) problem%threshold = threshold
      if (threshold_parameter_count(nclass, problem%threshold) <= 0) then
         status = 2
         return
      end if
      if (present(lambda)) problem%lambda = lambda
      if (present(estimate_lambda)) problem%estimate_lambda = estimate_lambda
      if (problem%estimate_lambda .and. problem%link /= link_aranda_ordaz .and. problem%link /= link_log_gamma) then
         status = 3
         return
      end if
      problem%y = y
      problem%x = x
      allocate(problem%weights(n), problem%offset(n), problem%scale_offset(n))
      problem%weights = 1.0_dp
      problem%offset = 0.0_dp
      problem%scale_offset = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = 4
            return
         end if
         problem%weights = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            status = 5
            return
         end if
         problem%offset = offset
      end if
      if (present(scale_offset)) then
         if (size(scale_offset) /= n) then
            status = 6
            return
         end if
         problem%scale_offset = scale_offset
      end if
      if (present(scale_x)) then
         if (size(scale_x, 1) /= n) then
            status = 7
            return
         end if
         problem%scale_x = scale_x
      else
         allocate(problem%scale_x(n, 0))
      end if
      if (present(nominal_x)) then
         if (size(nominal_x, 1) /= n) then
            status = 8
            return
         end if
         problem%nominal_x = nominal_x
      else
         allocate(problem%nominal_x(n, 0))
      end if
   end subroutine init_clm_problem

   function clm_nll(par, problem) result(nll)
      real(dp), intent(in) :: par(:) !! Packed threshold, location, scale, and optional flexible-link parameters.
      type(clm_problem), intent(in) :: problem !! Response, design matrices, link, and threshold specification.
      real(dp) :: nll
      nll = clm_nll_problem(par, problem)
   end function clm_nll

   function clm_value(self, x) result(nll)
      class(clm_problem), intent(in) :: self !! Cumulative-link problem evaluated as a numerical objective.
      real(dp), intent(in) :: x(:) !! Packed threshold, location, scale, and optional flexible-link parameters.
      real(dp) :: nll
      nll = clm_nll_problem(x, self)
   end function clm_value

   pure subroutine clm_parameter_counts(problem, nbase, nalpha, nbeta, nzeta, npar)
      type(clm_problem), intent(in) :: problem !! Model specification whose packed parameter dimensions are requested.
      integer, intent(out) :: nbase !! Number of free threshold parameters for one nominal reference pattern.
      integer, intent(out) :: nalpha !! Total threshold and nominal-threshold parameter count.
      integer, intent(out) :: nbeta !! Number of location coefficients.
      integer, intent(out) :: nzeta !! Number of log-scale coefficients.
      integer, intent(out) :: npar !! Total packed parameter count, including an estimated link shape when present.
      nbase = threshold_parameter_count(problem%nclass, problem%threshold)
      nalpha = nbase*(1 + size(problem%nominal_x, 2))
      nbeta = size(problem%x, 2)
      nzeta = size(problem%scale_x, 2)
      npar = nalpha + nbeta + nzeta + merge(1, 0, problem%estimate_lambda)
   end subroutine clm_parameter_counts

   pure subroutine nominal_thresholds_from_row(alpha, nominal_row, tjac, theta, status)
      real(dp), intent(in) :: alpha(:) !! Packed baseline and nominal threshold parameters in predictor-block order.
      real(dp), intent(in) :: nominal_row(:) !! Nominal predictor row without the implicit intercept.
      real(dp), intent(in) :: tjac(:, :) !! Threshold-structure Jacobian mapping one alpha block to actual cut points.
      real(dp), intent(out) :: theta(:) !! Observation-specific actual cut points.
      integer, intent(out) :: status !! Zero on success; nonzero for inconsistent dimensions.
      integer :: nbase, nblock, j, lo, hi
      nbase = size(tjac, 2)
      nblock = 1 + size(nominal_row)
      if (size(alpha) /= nbase*nblock .or. size(theta) /= size(tjac, 1)) then
         status = 1
         theta = 0.0_dp
         return
      end if
      theta = matmul(tjac, alpha(:nbase))
      do j = 1, size(nominal_row)
         lo = j*nbase + 1
         hi = (j + 1)*nbase
         theta = theta + nominal_row(j)*matmul(tjac, alpha(lo:hi))
      end do
      status = 0
   end subroutine nominal_thresholds_from_row

   subroutine clm_nominal_thresholds(alpha, nominal_row, nclass, threshold, theta, status)
      real(dp), intent(in) :: alpha(:) !! Packed baseline and nominal threshold parameters in predictor-block order.
      real(dp), intent(in) :: nominal_row(:) !! Nominal predictor values without the implicit intercept.
      integer, intent(in) :: nclass !! Number of ordered categories used to define the threshold Jacobian.
      integer, intent(in) :: threshold !! Threshold-structure identifier.
      real(dp), intent(out) :: theta(:) !! Cut points for the supplied nominal predictor pattern.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid threshold or parameter dimensions.
      real(dp), allocatable :: tjac(:, :)
      call threshold_jacobian(nclass, threshold, tjac, status)
      if (status /= 0) return
      call nominal_thresholds_from_row(alpha, nominal_row, tjac, theta, status)
   end subroutine clm_nominal_thresholds

   function clm_nll_problem(par, problem) result(nll)
      real(dp), intent(in) :: par(:) !! Packed model parameter vector.
      type(clm_problem), intent(in) :: problem !! Cumulative-link data and model specification.
      real(dp) :: nll
      real(dp), allocatable :: theta(:), tjac(:, :)
      real(dp) :: eta, sigma, prob, lambda
      integer :: nbase, nalpha, nbeta, nzeta, npar, i, pos, st
      call clm_parameter_counts(problem, nbase, nalpha, nbeta, nzeta, npar)
      if (size(par) /= npar) then
         nll = bad_objective
         return
      end if
      call threshold_jacobian(problem%nclass, problem%threshold, tjac, st)
      if (st /= 0) then
         nll = bad_objective
         return
      end if
      allocate(theta(problem%nclass - 1))
      pos = nalpha + nbeta + nzeta
      lambda = problem%lambda
      if (problem%estimate_lambda) then
         if (problem%link == link_aranda_ordaz) then
            lambda = exp(par(pos + 1))
         else
            lambda = par(pos + 1)
         end if
      end if
      nll = 0.0_dp
      do i = 1, size(problem%y)
         if (problem%weights(i) <= 0.0_dp) cycle
         call nominal_thresholds_from_row(par(:nalpha), problem%nominal_x(i, :), tjac, theta, st)
         if (st /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) then
            nll = bad_objective
            return
         end if
         eta = problem%offset(i)
         if (nbeta > 0) eta = eta + dot_product(problem%x(i, :), par(nalpha + 1:nalpha + nbeta))
         sigma = exp(problem%scale_offset(i))
         if (nzeta > 0) sigma = exp(problem%scale_offset(i) + &
              dot_product(problem%scale_x(i, :), par(nalpha + nbeta + 1:nalpha + nbeta + nzeta)))
         prob = category_probability(problem%y(i), theta, eta, sigma, problem%link, lambda)
         if (prob <= tiny(1.0_dp) .or. prob > 1.0_dp) then
            nll = bad_objective
            return
         end if
         nll = nll - problem%weights(i)*log(prob)
      end do
   end function clm_nll_problem

   pure real(dp) function category_probability(y, theta, eta, sigma, link, lambda) result(prob)
      integer, intent(in) :: y !! Ordered category code for one observation.
      real(dp), intent(in) :: theta(:) !! Increasing cut points separating adjacent response categories.
      real(dp), intent(in) :: eta !! Location linear predictor including offset.
      real(dp), intent(in) :: sigma !! Positive latent scale for the observation.
      integer, intent(in) :: link !! Link identifier controlling the latent error distribution.
      real(dp), intent(in) :: lambda !! Flexible-link shape parameter, ignored for fixed-shape links.
      real(dp) :: upper, lower
      integer :: k
      k = size(theta) + 1
      if (y == 1) then
         prob = link_cdf((theta(1) - eta)/sigma, link, lambda, .true.)
         return
      end if
      if (y == k) then
         prob = link_cdf((theta(k - 1) - eta)/sigma, link, lambda, .false.)
         return
      end if
      if ((theta(y - 1) - eta)/sigma > 0.0_dp) then
         upper = link_cdf((theta(y - 1) - eta)/sigma, link, lambda, .false.)
         lower = link_cdf((theta(y) - eta)/sigma, link, lambda, .false.)
         prob = upper - lower
      else
         upper = link_cdf((theta(y) - eta)/sigma, link, lambda, .true.)
         lower = link_cdf((theta(y - 1) - eta)/sigma, link, lambda, .true.)
         prob = upper - lower
      end if
   end function category_probability

   pure subroutine fill_threshold_basis(tjac, threshold_index, nominal_row, basis)
      real(dp), intent(in) :: tjac(:, :) !! Threshold Jacobian for one baseline threshold parameter block.
      integer, intent(in) :: threshold_index !! Actual cut-point row; zero means an infinite boundary.
      real(dp), intent(in) :: nominal_row(:) !! Nominal predictor row without its implicit intercept.
      real(dp), intent(out) :: basis(:) !! Derivative of the selected cut point with respect to packed alpha parameters.
      integer :: nbase, j, lo, hi
      basis = 0.0_dp
      if (threshold_index <= 0 .or. threshold_index > size(tjac, 1)) return
      nbase = size(tjac, 2)
      basis(:nbase) = tjac(threshold_index, :)
      do j = 1, size(nominal_row)
         lo = j*nbase + 1
         hi = (j + 1)*nbase
         basis(lo:hi) = nominal_row(j)*tjac(threshold_index, :)
      end do
   end subroutine fill_threshold_basis

   subroutine clm_analytic_gradient(par, problem, gradient, status)
      real(dp), intent(in) :: par(:) !! Packed CLM parameters; link-shape estimation is not supported analytically.
      type(clm_problem), intent(in) :: problem !! Cumulative-link model data and fixed link specification.
      real(dp), intent(out) :: gradient(:) !! Exact negative-log-likelihood gradient for threshold, location, and scale terms.
      integer, intent(out) :: status !! Zero on success; nonzero for dimensions, shape estimation, or invalid probabilities.
      real(dp), allocatable :: hessian(:, :)
      real(dp) :: nll
      allocate(hessian(size(par), size(par)))
      call clm_nll_derivatives(par, problem, nll, gradient, hessian, .false., status)
   end subroutine clm_analytic_gradient

   subroutine clm_analytic_hessian(par, problem, hessian, status)
      real(dp), intent(in) :: par(:) !! Packed CLM parameters; link-shape estimation is not supported analytically.
      type(clm_problem), intent(in) :: problem !! Cumulative-link model data and fixed link specification.
      real(dp), intent(out) :: hessian(:, :) !! Exact negative-log-likelihood Hessian for threshold, location, and scale terms.
      integer, intent(out) :: status !! Zero on success; nonzero for dimensions, shape estimation, or invalid probabilities.
      real(dp), allocatable :: gradient(:)
      real(dp) :: nll
      allocate(gradient(size(par)))
      call clm_nll_derivatives(par, problem, nll, gradient, hessian, .true., status)
   end subroutine clm_analytic_hessian

   subroutine clm_nll_derivatives(par, problem, nll, gradient, hessian, need_hessian, status)
      real(dp), intent(in) :: par(:) !! Packed threshold, location, and scale parameters for a fixed link shape.
      type(clm_problem), intent(in) :: problem !! Cumulative-link data and design matrices.
      real(dp), intent(out) :: nll !! Negative log-likelihood at the supplied parameters.
      real(dp), intent(out) :: gradient(:) !! Exact gradient of nll, with the same length as par.
      real(dp), intent(out) :: hessian(:, :) !! Exact Hessian of nll when need_hessian is true; otherwise zeroed.
      logical, intent(in) :: need_hessian !! Compute second derivatives when true.
      integer, intent(out) :: status !! Zero on success; nonzero for unsupported or invalid parameter states.
      real(dp), allocatable :: tjac(:, :), theta(:), b1(:), b2(:), c2(:), c3(:)
      real(dp) :: eta, sigma, eta1, eta2, p1, p2, g1, g2, prob, wpi, wpi2, aterm, bterm, epg1, epg2
      integer :: nbase, nalpha, nbeta, nzeta, npar, npsi, i, j, k, st, upper_index, lower_index
      call clm_parameter_counts(problem, nbase, nalpha, nbeta, nzeta, npar)
      status = 0
      gradient = 0.0_dp
      hessian = 0.0_dp
      nll = bad_objective
      if (problem%estimate_lambda) then
         status = 1
         return
      end if
      if (size(par) /= npar .or. size(gradient) /= npar .or. size(hessian, 1) /= npar .or. size(hessian, 2) /= npar) then
         status = 2
         return
      end if
      npsi = nalpha + nbeta
      call threshold_jacobian(problem%nclass, problem%threshold, tjac, st)
      if (st /= 0) then
         status = 3
         return
      end if
      allocate(theta(problem%nclass - 1), b1(npsi), b2(npsi), c2(npsi), c3(nzeta))
      nll = 0.0_dp
      do i = 1, size(problem%y)
         if (problem%weights(i) <= 0.0_dp) cycle
         call nominal_thresholds_from_row(par(:nalpha), problem%nominal_x(i, :), tjac, theta, st)
         if (st /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) then
            status = 4
            nll = bad_objective
            return
         end if
         eta = problem%offset(i)
         if (nbeta > 0) eta = eta + dot_product(problem%x(i, :), par(nalpha + 1:nalpha + nbeta))
         sigma = exp(problem%scale_offset(i))
         if (nzeta > 0) sigma = exp(problem%scale_offset(i) + &
              dot_product(problem%scale_x(i, :), par(nalpha + nbeta + 1:nalpha + nbeta + nzeta)))
         prob = category_probability(problem%y(i), theta, eta, sigma, problem%link, problem%lambda)
         if (prob <= tiny(1.0_dp) .or. prob > 1.0_dp) then
            status = 5
            nll = bad_objective
            return
         end if
         nll = nll - problem%weights(i)*log(prob)
         upper_index = problem%y(i)
         lower_index = problem%y(i) - 1
         b1 = 0.0_dp
         b2 = 0.0_dp
         call fill_threshold_basis(tjac, upper_index, problem%nominal_x(i, :), b1(:nalpha))
         call fill_threshold_basis(tjac, lower_index, problem%nominal_x(i, :), b2(:nalpha))
         if (nbeta > 0) then
            b1(nalpha + 1:npsi) = -problem%x(i, :)
            b2(nalpha + 1:npsi) = -problem%x(i, :)
         end if
         if (upper_index <= size(theta)) then
            eta1 = (theta(upper_index) - eta)/sigma
            p1 = link_pdf(eta1, problem%link, problem%lambda)
            g1 = link_pdf_gradient(eta1, problem%link, problem%lambda)
         else
            eta1 = 0.0_dp
            p1 = 0.0_dp
            g1 = 0.0_dp
         end if
         if (lower_index >= 1) then
            eta2 = (theta(lower_index) - eta)/sigma
            p2 = link_pdf(eta2, problem%link, problem%lambda)
            g2 = link_pdf_gradient(eta2, problem%link, problem%lambda)
         else
            eta2 = 0.0_dp
            p2 = 0.0_dp
            g2 = 0.0_dp
         end if
         c2 = (b1*p1 - b2*p2)/sigma
         if (nzeta > 0) c3 = -(eta1*p1 - eta2*p2)*problem%scale_x(i, :)
         wpi = problem%weights(i)/prob
         wpi2 = problem%weights(i)/(prob*prob)
         gradient(:npsi) = gradient(:npsi) - wpi*c2
         if (nzeta > 0) gradient(npsi + 1:npar) = gradient(npsi + 1:npar) - wpi*c3
         if (.not. need_hessian) cycle
         do j = 1, npsi
            do k = 1, npsi
               hessian(j, k) = hessian(j, k) + wpi2*c2(j)*c2(k) &
                    - wpi*g1*b1(j)*b1(k)/(sigma*sigma) + wpi*g2*b2(j)*b2(k)/(sigma*sigma)
            end do
         end do
         if (nzeta <= 0) cycle
         epg1 = p1 + g1*eta1
         epg2 = p2 + g2*eta2
         do j = 1, npsi
            do k = 1, nzeta
               hessian(j, npsi + k) = hessian(j, npsi + k) &
                    + wpi*b1(j)*epg1*problem%scale_x(i, k)/sigma &
                    - wpi*b2(j)*epg2*problem%scale_x(i, k)/sigma + wpi2*c2(j)*c3(k)
               hessian(npsi + k, j) = hessian(j, npsi + k)
            end do
         end do
         aterm = eta1*p1 - eta2*p2
         bterm = eta1*epg1 - eta2*epg2
         do j = 1, nzeta
            do k = 1, nzeta
               hessian(npsi + j, npsi + k) = hessian(npsi + j, npsi + k) &
                    + wpi*(aterm*aterm/prob - bterm)*problem%scale_x(i, j)*problem%scale_x(i, k)
            end do
         end do
      end do
   end subroutine clm_nll_derivatives

   subroutine newton_fit_clm(problem, par, fval, iterations, status, max_iter, grad_tol)
      type(clm_problem), intent(in) :: problem !! Fixed-shape cumulative-link problem with analytic derivatives.
      real(dp), intent(inout) :: par(:) !! Starting parameter vector on input and fitted parameter vector on output.
      real(dp), intent(out) :: fval !! Negative log-likelihood at the returned parameter vector.
      integer, intent(out) :: iterations !! Number of accepted Newton iterations.
      integer, intent(out) :: status !! Zero on convergence; one for limit, two for line search, three for derivative failure.
      integer, intent(in) :: max_iter !! Maximum Newton iterations.
      real(dp), intent(in) :: grad_tol !! Infinity-norm gradient convergence tolerance.
      real(dp), allocatable :: grad(:), hess(:, :), hwork(:, :), chol(:, :), step(:), trial(:), eig(:)
      real(dp) :: trial_f, alpha_step, mineig, inflation, scale
      integer :: n, iter, ls, st, chst, solst, i
      n = size(par)
      allocate(grad(n), hess(n, n), hwork(n, n), chol(n, n), step(n), trial(n), eig(n))
      status = 1
      iterations = 0
      do iter = 1, max_iter
         call clm_nll_derivatives(par, problem, fval, grad, hess, .true., st)
         if (st /= 0) then
            status = 3
            return
         end if
         if (maxval(abs(grad)) <= grad_tol) then
            status = 0
            return
         end if
         hwork = hess
         call cholesky_factor(hwork, chol, chst)
         if (chst /= 0) then
            call symmetric_eigenvalues(hess, eig, st)
            if (st /= 0) then
               status = 3
               return
            end if
            mineig = minval(eig)
            scale = max(1.0_dp, maxval(abs(eig)))
            inflation = max(1.0e-8_dp*scale, -mineig + 1.0e-8_dp*scale)
            do i = 1, n
               hwork(i, i) = hwork(i, i) + inflation
            end do
            call cholesky_factor(hwork, chol, chst)
            if (chst /= 0) then
               status = 3
               return
            end if
         end if
         call solve_cholesky(chol, grad, step, solst)
         if (solst /= 0) then
            status = 3
            return
         end if
         alpha_step = 1.0_dp
         do ls = 1, 40
            trial = par - alpha_step*step
            trial_f = clm_nll_problem(trial, problem)
            if (trial_f < fval) exit
            alpha_step = 0.5_dp*alpha_step
         end do
         if (ls > 40) then
            status = 2
            return
         end if
         par = trial
         fval = trial_f
         iterations = iter
         if (maxval(abs(alpha_step*step)) <= 1.0e-10_dp*(1.0_dp + maxval(abs(par)))) then
            call clm_nll_derivatives(par, problem, fval, grad, hess, .false., st)
            if (st == 0 .and. maxval(abs(grad)) <= 10.0_dp*grad_tol) then
               status = 0
               return
            end if
         end if
      end do
   end subroutine newton_fit_clm

   subroutine fit_clm(problem, fit, start, max_iter, grad_tol)
      type(clm_problem), intent(in) :: problem !! Fully initialized cumulative-link fitting problem.
      type(clm_fit_result), intent(out) :: fit !! Fitted parameters, likelihood, derivatives, diagnostics, and covariance.
      real(dp), intent(in), optional :: start(:) !! Optional packed starting vector; otherwise ordinal-style starts are generated.
      integer, intent(in), optional :: max_iter !! Maximum optimizer iterations; defaults to 300.
      real(dp), intent(in), optional :: grad_tol !! Infinity-norm gradient tolerance; defaults to 1e-7.
      real(dp), allocatable :: par(:), alpha0(:)
      integer :: nbase, nalpha, nbeta, nzeta, npar, st, invst, limit, diagst, j, lo, hi
      real(dp) :: tol, fval
      call clm_parameter_counts(problem, nbase, nalpha, nbeta, nzeta, npar)
      allocate(par(npar))
      if (present(start)) then
         if (size(start) /= npar) then
            fit%status = 10
            return
         end if
         par = start
      else
         call threshold_start(problem%nclass, problem%threshold, alpha0, st)
         if (st /= 0) then
            fit%status = 11
            return
         end if
         par = 0.0_dp
         par(:nbase) = alpha0
         if (problem%estimate_lambda) then
            if (problem%link == link_aranda_ordaz) then
               par(npar) = log(max(problem%lambda, 1.0e-3_dp))
            else
               par(npar) = problem%lambda
            end if
         end if
      end if
      limit = 300
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-7_dp
      if (present(grad_tol)) tol = grad_tol
      if (problem%estimate_lambda) then
         call bfgs_minimize(problem, par, fval, fit%iterations, fit%status, limit, tol)
      else
         call newton_fit_clm(problem, par, fval, fit%iterations, fit%status, limit, tol)
      end if
      fit%converged = fit%status == 0
      fit%loglik = -fval
      allocate(fit%alpha(nalpha), fit%nominal_alpha(nbase, 1 + size(problem%nominal_x, 2)), fit%beta(nbeta), &
               fit%zeta(nzeta), fit%gradient(npar), fit%hessian(npar, npar), fit%vcov(npar, npar), &
               fit%fitted(size(problem%y)))
      fit%alpha = par(:nalpha)
      do j = 0, size(problem%nominal_x, 2)
         lo = j*nbase + 1
         hi = (j + 1)*nbase
         fit%nominal_alpha(:, j + 1) = fit%alpha(lo:hi)
      end do
      if (nbeta > 0) fit%beta = par(nalpha + 1:nalpha + nbeta)
      if (nzeta > 0) fit%zeta = par(nalpha + nbeta + 1:nalpha + nbeta + nzeta)
      fit%lambda = problem%lambda
      if (problem%estimate_lambda) then
         if (problem%link == link_aranda_ordaz) then
            fit%lambda = exp(par(npar))
         else
            fit%lambda = par(npar)
         end if
      end if
      call thresholds_from_alpha(fit%alpha(:nbase), problem%nclass, problem%threshold, fit%theta, st)
      if (problem%estimate_lambda) then
         call numerical_gradient(problem, par, fit%gradient, 1.0e-6_dp)
         call numerical_hessian(problem, par, fit%hessian, 1.0e-4_dp)
      else
         call clm_nll_derivatives(par, problem, fval, fit%gradient, fit%hessian, .true., st)
         if (st /= 0) then
            call numerical_gradient(problem, par, fit%gradient, 1.0e-6_dp)
            call numerical_hessian(problem, par, fit%hessian, 1.0e-4_dp)
         end if
      end if
      fit%max_gradient = maxval(abs(fit%gradient))
      call invert_matrix(fit%hessian, fit%vcov, invst)
      if (invst /= 0) fit%vcov = 0.0_dp
      call hessian_diagnostics(fit%hessian, fit%hessian_min_eigenvalue, fit%hessian_max_eigenvalue, &
           fit%hessian_condition, fit%hessian_rank, fit%hessian_positive_definite, diagst)
      if (diagst /= 0) fit%hessian_positive_definite = .false.
      call observed_fitted(problem, fit%alpha, fit%beta, fit%zeta, fit%lambda, fit%fitted)
   end subroutine fit_clm

   subroutine observed_fitted(problem, alpha, beta, zeta, lambda, fitted)
      type(clm_problem), intent(in) :: problem !! Cumulative-link data and design information.
      real(dp), intent(in) :: alpha(:) !! Packed baseline and nominal threshold parameters.
      real(dp), intent(in) :: beta(:) !! Fitted location coefficients.
      real(dp), intent(in) :: zeta(:) !! Fitted log-scale coefficients.
      real(dp), intent(in) :: lambda !! Fitted or fixed flexible-link shape parameter.
      real(dp), intent(out) :: fitted(:) !! Probability assigned to each observation's realized category.
      real(dp), allocatable :: tjac(:, :), theta(:)
      real(dp) :: eta, sigma
      integer :: i, st
      call threshold_jacobian(problem%nclass, problem%threshold, tjac, st)
      if (st /= 0) then
         fitted = 0.0_dp
         return
      end if
      allocate(theta(problem%nclass - 1))
      do i = 1, size(problem%y)
         call nominal_thresholds_from_row(alpha, problem%nominal_x(i, :), tjac, theta, st)
         if (st /= 0) then
            fitted(i) = 0.0_dp
            cycle
         end if
         eta = problem%offset(i)
         if (size(beta) > 0) eta = eta + dot_product(problem%x(i, :), beta)
         sigma = exp(problem%scale_offset(i))
         if (size(zeta) > 0) sigma = exp(problem%scale_offset(i) + dot_product(problem%scale_x(i, :), zeta))
         fitted(i) = category_probability(problem%y(i), theta, eta, sigma, problem%link, lambda)
      end do
   end subroutine observed_fitted

   subroutine clm_predict_proba(theta, beta, x, link, probabilities, offset, sigma, lambda, status)
      real(dp), intent(in) :: theta(:) !! Increasing fitted cut points; number of categories is size(theta)+1.
      real(dp), intent(in) :: beta(:) !! Location coefficients matching the columns of x.
      real(dp), intent(in) :: x(:, :) !! New location design matrix without an intercept.
      integer, intent(in) :: link !! Link identifier controlling the latent error distribution.
      real(dp), intent(out) :: probabilities(:, :) !! Category probabilities with rows matching x and columns ordered categories.
      real(dp), intent(in), optional :: offset(:) !! Optional location offsets for prediction rows.
      real(dp), intent(in), optional :: sigma(:) !! Optional positive latent scales for prediction rows; defaults to one.
      real(dp), intent(in), optional :: lambda !! Flexible-link shape parameter; defaults to one.
      integer, intent(out) :: status !! Zero on success; nonzero for dimension or positivity errors.
      real(dp) :: eta, sc, lam
      integer :: i, j, n, k
      n = size(x, 1)
      k = size(theta) + 1
      status = 0
      if (size(x, 2) /= size(beta) .or. size(probabilities, 1) /= n .or. size(probabilities, 2) /= k) then
         status = 1
         return
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            status = 2
            return
         end if
      end if
      if (present(sigma)) then
         if (size(sigma) /= n .or. any(sigma <= 0.0_dp)) then
            status = 3
            return
         end if
      end if
      lam = 1.0_dp
      if (present(lambda)) lam = lambda
      do i = 1, n
         eta = dot_product(x(i, :), beta)
         if (present(offset)) eta = eta + offset(i)
         sc = 1.0_dp
         if (present(sigma)) sc = sigma(i)
         do j = 1, k
            probabilities(i, j) = category_probability(j, theta, eta, sc, link, lam)
         end do
      end do
   end subroutine clm_predict_proba

   subroutine clm_predict_proba_nominal(alpha, beta, x, nominal_x, nclass, threshold, link, probabilities, offset, &
                                        sigma, lambda, status)
      real(dp), intent(in) :: alpha(:) !! Packed baseline and nominal threshold coefficients from a fitted nominal CLM.
      real(dp), intent(in) :: beta(:) !! Location coefficients matching columns of x.
      real(dp), intent(in) :: x(:, :) !! New location design matrix without an intercept.
      real(dp), intent(in) :: nominal_x(:, :) !! New nominal design matrix without an intercept.
      integer, intent(in) :: nclass !! Number of ordered response categories.
      integer, intent(in) :: threshold !! Threshold-structure identifier used for fitting.
      integer, intent(in) :: link !! Link identifier controlling the latent error distribution.
      real(dp), intent(out) :: probabilities(:, :) !! Row-wise probabilities for all ordered categories.
      real(dp), intent(in), optional :: offset(:) !! Optional location offsets for prediction rows.
      real(dp), intent(in), optional :: sigma(:) !! Optional positive latent scales for prediction rows; defaults to one.
      real(dp), intent(in), optional :: lambda !! Flexible-link shape parameter; defaults to one.
      integer, intent(out) :: status !! Zero on success; nonzero for inconsistent dimensions or unordered cut points.
      real(dp), allocatable :: tjac(:, :), theta(:)
      real(dp) :: eta, sc, lam
      integer :: n, i, j, st
      n = size(x, 1)
      status = 0
      if (size(nominal_x, 1) /= n .or. size(x, 2) /= size(beta) .or. &
          size(probabilities, 1) /= n .or. size(probabilities, 2) /= nclass) then
         status = 1
         return
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            status = 2
            return
         end if
      end if
      if (present(sigma)) then
         if (size(sigma) /= n .or. any(sigma <= 0.0_dp)) then
            status = 3
            return
         end if
      end if
      call threshold_jacobian(nclass, threshold, tjac, st)
      if (st /= 0 .or. size(alpha) /= size(tjac, 2)*(1 + size(nominal_x, 2))) then
         status = 4
         return
      end if
      allocate(theta(nclass - 1))
      lam = 1.0_dp
      if (present(lambda)) lam = lambda
      do i = 1, n
         call nominal_thresholds_from_row(alpha, nominal_x(i, :), tjac, theta, st)
         if (st /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) then
            status = 5
            return
         end if
         eta = dot_product(x(i, :), beta)
         if (present(offset)) eta = eta + offset(i)
         sc = 1.0_dp
         if (present(sigma)) sc = sigma(i)
         do j = 1, nclass
            probabilities(i, j) = category_probability(j, theta, eta, sc, link, lam)
         end do
      end do
   end subroutine clm_predict_proba_nominal
end module ordinal_clm
