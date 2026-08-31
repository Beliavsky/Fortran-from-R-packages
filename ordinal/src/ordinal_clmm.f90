! Mixed cumulative-link model kernels based on ordinal/R/clmm*.R and clmm.ssr.R.
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_clmm
   use ordinal_kinds, only : dp
   use ordinal_links, only : link_cdf, link_pdf, link_pdf_gradient, link_logit
   use ordinal_thresholds, only : threshold_flexible, threshold_parameter_count, threshold_start, thresholds_from_alpha
   use ordinal_numerics, only : objective_type, bfgs_minimize, numerical_gradient, numerical_hessian, invert_matrix, &
                                hessian_diagnostics
   use ordinal_quadrature, only : gauss_hermite_rule
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: bad_objective = huge(1.0_dp)/1000.0_dp
   type, extends(objective_type), public :: clmm_problem
      integer, allocatable :: y(:)
      integer, allocatable :: group(:)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: offset(:)
      integer :: nclass = 0
      integer :: ngroup = 0
      integer :: link = link_logit
      integer :: threshold = threshold_flexible
      integer :: nAGQ = 1
      real(dp) :: lambda = 1.0_dp
   contains
      procedure :: value => clmm_value
   end type clmm_problem
   type, public :: clmm_fit_result
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: vcov(:, :)
      real(dp), allocatable :: ranef_mode(:)
      real(dp), allocatable :: ranef_condvar(:)
      real(dp) :: tau = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: max_gradient = huge(1.0_dp)
      real(dp) :: hessian_min_eigenvalue = 0.0_dp
      real(dp) :: hessian_max_eigenvalue = 0.0_dp
      real(dp) :: hessian_condition = huge(1.0_dp)
      integer :: hessian_rank = 0
      integer :: nAGQ = 1
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      logical :: hessian_positive_definite = .false.
   end type clmm_fit_result
   public :: init_clmm_problem, clmm_nll, fit_clmm, clmm_conditional_modes
contains
   subroutine init_clmm_problem(problem, y, x, group, nclass, link, threshold, weights, offset, lambda, nAGQ, status)
      type(clmm_problem), intent(out) :: problem !! Initialized single random-intercept cumulative-link mixed-model problem.
      integer, intent(in) :: y(:) !! Ordered response codes in 1:nclass, one per row of x.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect location design matrix without an intercept.
      integer, intent(in) :: group(:) !! Positive group codes; arbitrary labels are remapped internally to 1:ngroup.
      integer, intent(in) :: nclass !! Number of ordered response categories; must be at least two.
      integer, intent(in), optional :: link !! Link identifier; defaults to logistic.
      integer, intent(in), optional :: threshold !! Threshold structure identifier; defaults to flexible.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights; defaults to one.
      real(dp), intent(in), optional :: offset(:) !! Fixed location offsets, one per observation; defaults to zero.
      real(dp), intent(in), optional :: lambda !! Fixed flexible-link shape parameter; defaults to one.
      integer, intent(in), optional :: nAGQ !! Quadrature control: 0/1 Laplace, >1 adaptive GHQ, <0 nonadaptive GHQ.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, codes, weights, or options.
      integer, allocatable :: labels(:)
      integer :: n, i, j, ng
      n = size(y)
      status = 0
      if (size(x, 1) /= n .or. size(group) /= n .or. nclass < 2 .or. any(y < 1) .or. any(y > nclass) .or. any(group <= 0)) then
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
      if (present(nAGQ)) problem%nAGQ = nAGQ
      if (problem%nAGQ < 0 .and. abs(problem%nAGQ) < 1) then
         status = 3
         return
      end if
      problem%y = y
      problem%x = x
      allocate(problem%group(n), labels(n), problem%weights(n), problem%offset(n))
      problem%weights = 1.0_dp
      problem%offset = 0.0_dp
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
      ng = 0
      labels = 0
      do i = 1, n
         j = find_label(group(i), labels, ng)
         if (j == 0) then
            ng = ng + 1
            labels(ng) = group(i)
            j = ng
         end if
         problem%group(i) = j
      end do
      problem%ngroup = ng
   contains
      pure integer function find_label(value, values, nused) result(pos)
         integer, intent(in) :: value !! Original positive grouping label to locate.
         integer, intent(in) :: values(:) !! Previously encountered original group labels.
         integer, intent(in) :: nused !! Number of valid entries at the front of values.
         integer :: k
         pos = 0
         do k = 1, nused
            if (values(k) == value) then
               pos = k
               return
            end if
         end do
      end function find_label
   end subroutine init_clmm_problem

   function clmm_nll(par, problem) result(nll)
      real(dp), intent(in) :: par(:) !! Packed alpha, fixed-effect beta, and log random-intercept standard deviation.
      type(clmm_problem), intent(in) :: problem !! Response, grouping, design, link, threshold, and quadrature specification.
      real(dp) :: nll
      nll = clmm_nll_problem(par, problem)
   end function clmm_nll

   function clmm_value(self, x) result(nll)
      class(clmm_problem), intent(in) :: self !! Mixed cumulative-link problem evaluated as a numerical objective.
      real(dp), intent(in) :: x(:) !! Packed alpha, fixed-effect beta, and log random-intercept standard deviation.
      real(dp) :: nll
      nll = clmm_nll_problem(x, self)
   end function clmm_value

   function clmm_nll_problem(par, problem) result(nll)
      real(dp), intent(in) :: par(:) !! Packed mixed-model parameter vector.
      type(clmm_problem), intent(in) :: problem !! Single random-intercept mixed-model data and quadrature options.
      real(dp) :: nll
      real(dp), allocatable :: theta(:), nodes(:), qweights(:), terms(:), beta(:)
      real(dp) :: tau, mode, curvature, qmode, scale, b, loglik, maxterm, log_integral
      integer :: nalpha, nbeta, g, q, st, nq
      nalpha = threshold_parameter_count(problem%nclass, problem%threshold)
      nbeta = size(problem%x, 2)
      if (size(par) /= nalpha + nbeta + 1) then
         nll = bad_objective
         return
      end if
      call thresholds_from_alpha(par(:nalpha), problem%nclass, problem%threshold, theta, st)
      if (st /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) then
         nll = bad_objective
         return
      end if
      allocate(beta(nbeta))
      if (nbeta > 0) beta = par(nalpha + 1:nalpha + nbeta)
      tau = exp(par(size(par)))
      if (tau < 1.0e-10_dp .or. tau > 1.0e4_dp) then
         nll = bad_objective
         return
      end if
      nll = 0.0_dp
      if (problem%nAGQ < 0) then
         nq = abs(problem%nAGQ)
         call gauss_hermite_rule(nq, nodes, qweights, st)
         if (st /= 0) then
            nll = bad_objective
            return
         end if
         allocate(terms(nq))
         do g = 1, problem%ngroup
            do q = 1, nq
               b = sqrt(2.0_dp)*tau*nodes(q)
               call group_data_loglik(problem, theta, beta, g, b, loglik, st)
               if (st /= 0) then
                  nll = bad_objective
                  return
               end if
               terms(q) = log(qweights(q)) + loglik
            end do
            maxterm = maxval(terms)
            nll = nll - (maxterm + log(sum(exp(terms - maxterm))) - 0.5_dp*log(pi))
         end do
         return
      end if
      do g = 1, problem%ngroup
         call group_mode_curvature(problem, theta, beta, tau, g, mode, curvature, qmode, st)
         if (st /= 0 .or. curvature <= 0.0_dp) then
            nll = bad_objective
            return
         end if
         if (problem%nAGQ <= 1) then
            nll = nll + qmode + log(tau) + 0.5_dp*log(curvature)
         else
            nq = problem%nAGQ
            if (.not. allocated(nodes)) then
               call gauss_hermite_rule(nq, nodes, qweights, st)
               if (st /= 0) then
                  nll = bad_objective
                  return
               end if
               allocate(terms(nq))
            end if
            scale = 1.0_dp/sqrt(curvature)
            do q = 1, nq
               b = mode + sqrt(2.0_dp)*scale*nodes(q)
               call group_data_loglik(problem, theta, beta, g, b, loglik, st)
               if (st /= 0) then
                  nll = bad_objective
                  return
               end if
               terms(q) = log(qweights(q)) + loglik - 0.5_dp*(b/tau)**2 + nodes(q)*nodes(q)
            end do
            maxterm = maxval(terms)
            log_integral = -0.5_dp*log(pi) - log(tau) + log(scale) &
                 + maxterm + log(sum(exp(terms - maxterm)))
            nll = nll - log_integral
         end if
      end do
   end function clmm_nll_problem

   pure subroutine category_probability_derivatives(y, theta, eta, link, lambda, prob, dpeta, d2peta)
      integer, intent(in) :: y !! Ordered category code for one observation.
      real(dp), intent(in) :: theta(:) !! Increasing cut points separating adjacent categories.
      real(dp), intent(in) :: eta !! Fixed plus random location predictor.
      integer, intent(in) :: link !! Link identifier controlling the latent error distribution.
      real(dp), intent(in) :: lambda !! Fixed flexible-link shape parameter.
      real(dp), intent(out) :: prob !! Probability of category y.
      real(dp), intent(out) :: dpeta !! First derivative of category probability with respect to eta.
      real(dp), intent(out) :: d2peta !! Second derivative of category probability with respect to eta.
      real(dp) :: p1, p2, g1, g2
      integer :: k
      k = size(theta) + 1
      if (y == 1) then
         prob = link_cdf(theta(1) - eta, link, lambda, .true.)
         p1 = link_pdf(theta(1) - eta, link, lambda)
         g1 = link_pdf_gradient(theta(1) - eta, link, lambda)
         dpeta = -p1
         d2peta = g1
      else if (y == k) then
         prob = link_cdf(theta(k - 1) - eta, link, lambda, .false.)
         p2 = link_pdf(theta(k - 1) - eta, link, lambda)
         g2 = link_pdf_gradient(theta(k - 1) - eta, link, lambda)
         dpeta = p2
         d2peta = -g2
      else
         if (theta(y - 1) - eta > 0.0_dp) then
            prob = link_cdf(theta(y - 1) - eta, link, lambda, .false.) &
                 - link_cdf(theta(y) - eta, link, lambda, .false.)
         else
            prob = link_cdf(theta(y) - eta, link, lambda, .true.) &
                 - link_cdf(theta(y - 1) - eta, link, lambda, .true.)
         end if
         p1 = link_pdf(theta(y) - eta, link, lambda)
         p2 = link_pdf(theta(y - 1) - eta, link, lambda)
         g1 = link_pdf_gradient(theta(y) - eta, link, lambda)
         g2 = link_pdf_gradient(theta(y - 1) - eta, link, lambda)
         dpeta = p2 - p1
         d2peta = g1 - g2
      end if
   end subroutine category_probability_derivatives

   subroutine group_data_loglik(problem, theta, beta, g, b, loglik, status)
      type(clmm_problem), intent(in) :: problem !! Single random-intercept mixed-model data and design.
      real(dp), intent(in) :: theta(:) !! Fitted or candidate increasing cut points.
      real(dp), intent(in) :: beta(:) !! Candidate fixed-effect location coefficients.
      integer, intent(in) :: g !! Internally remapped group index in 1:problem%ngroup.
      real(dp), intent(in) :: b !! Candidate random intercept on the response latent-location scale.
      real(dp), intent(out) :: loglik !! Conditional data log-likelihood for group g at random intercept b.
      integer, intent(out) :: status !! Zero on success; nonzero if a category probability is invalid.
      real(dp) :: eta, prob, dpeta, d2peta
      integer :: i
      loglik = 0.0_dp
      status = 0
      do i = 1, size(problem%y)
         if (problem%group(i) /= g .or. problem%weights(i) <= 0.0_dp) cycle
         eta = problem%offset(i) + b
         if (size(beta) > 0) eta = eta + dot_product(problem%x(i, :), beta)
         call category_probability_derivatives(problem%y(i), theta, eta, problem%link, problem%lambda, prob, dpeta, d2peta)
         if (prob <= tiny(1.0_dp) .or. prob > 1.0_dp) then
            status = 1
            loglik = -bad_objective
            return
         end if
         loglik = loglik + problem%weights(i)*log(prob)
      end do
   end subroutine group_data_loglik

   subroutine group_penalized_derivatives(problem, theta, beta, tau, g, b, qvalue, gradient, curvature, status)
      type(clmm_problem), intent(in) :: problem !! Single random-intercept mixed-model data and design.
      real(dp), intent(in) :: theta(:) !! Candidate increasing cut points.
      real(dp), intent(in) :: beta(:) !! Candidate fixed-effect location coefficients.
      real(dp), intent(in) :: tau !! Positive Gaussian random-intercept standard deviation.
      integer, intent(in) :: g !! Internally remapped group index in 1:problem%ngroup.
      real(dp), intent(in) :: b !! Candidate random intercept.
      real(dp), intent(out) :: qvalue !! Penalized negative conditional log-likelihood, excluding normalizing constants.
      real(dp), intent(out) :: gradient !! First derivative of qvalue with respect to b.
      real(dp), intent(out) :: curvature !! Second derivative of qvalue with respect to b.
      integer, intent(out) :: status !! Zero on success; nonzero if a category probability is invalid.
      real(dp) :: eta, prob, dpeta, d2peta, w
      integer :: i
      qvalue = 0.5_dp*(b/tau)**2
      gradient = b/(tau*tau)
      curvature = 1.0_dp/(tau*tau)
      status = 0
      do i = 1, size(problem%y)
         if (problem%group(i) /= g .or. problem%weights(i) <= 0.0_dp) cycle
         eta = problem%offset(i) + b
         if (size(beta) > 0) eta = eta + dot_product(problem%x(i, :), beta)
         call category_probability_derivatives(problem%y(i), theta, eta, problem%link, problem%lambda, prob, dpeta, d2peta)
         if (prob <= tiny(1.0_dp) .or. prob > 1.0_dp) then
            status = 1
            qvalue = bad_objective
            return
         end if
         w = problem%weights(i)
         qvalue = qvalue - w*log(prob)
         gradient = gradient - w*dpeta/prob
         curvature = curvature + w*((dpeta/prob)**2 - d2peta/prob)
      end do
   end subroutine group_penalized_derivatives

   subroutine group_mode_curvature(problem, theta, beta, tau, g, mode, curvature, qmode, status)
      type(clmm_problem), intent(in) :: problem !! Single random-intercept mixed-model data and design.
      real(dp), intent(in) :: theta(:) !! Candidate increasing cut points.
      real(dp), intent(in) :: beta(:) !! Candidate fixed-effect location coefficients.
      real(dp), intent(in) :: tau !! Positive Gaussian random-intercept standard deviation.
      integer, intent(in) :: g !! Internally remapped group index in 1:problem%ngroup.
      real(dp), intent(out) :: mode !! Conditional random-intercept mode for group g.
      real(dp), intent(out) :: curvature !! Penalized nll curvature at the returned mode.
      real(dp), intent(out) :: qmode !! Penalized nll at the returned mode, excluding normalizing constants.
      integer, intent(out) :: status !! Zero on convergence; nonzero for invalid curvature, probability, or line search.
      real(dp) :: b, grad, hess, step, trial, trial_q, trial_g, trial_h, alpha_step
      integer :: iter, ls, st
      b = 0.0_dp
      status = 0
      do iter = 1, 80
         call group_penalized_derivatives(problem, theta, beta, tau, g, b, qmode, grad, hess, st)
         if (st /= 0 .or. hess <= 0.0_dp) then
            status = 1
            return
         end if
         if (abs(grad) <= 1.0e-9_dp*(1.0_dp + abs(qmode))) exit
         step = grad/hess
         step = max(-5.0_dp*tau, min(5.0_dp*tau, step))
         alpha_step = 1.0_dp
         do ls = 1, 40
            trial = b - alpha_step*step
            call group_penalized_derivatives(problem, theta, beta, tau, g, trial, trial_q, trial_g, trial_h, st)
            if (st == 0 .and. trial_q < qmode) exit
            alpha_step = 0.5_dp*alpha_step
         end do
         if (ls > 40) then
            if (abs(grad) <= 1.0e-7_dp) exit
            status = 2
            return
         end if
         b = trial
         if (abs(alpha_step*step) <= 1.0e-10_dp*(1.0_dp + abs(b))) exit
      end do
      call group_penalized_derivatives(problem, theta, beta, tau, g, b, qmode, grad, curvature, st)
      if (st /= 0 .or. curvature <= 0.0_dp) then
         status = 3
         return
      end if
      mode = b
   end subroutine group_mode_curvature

   subroutine fit_clmm(problem, fit, start, max_iter, grad_tol)
      type(clmm_problem), intent(in) :: problem !! Fully initialized single random-intercept cumulative-link mixed-model problem.
      type(clmm_fit_result), intent(out) :: fit !! Fitted fixed effects, thresholds, random-effect scale, modes, and diagnostics.
      real(dp), intent(in), optional :: start(:) !! Optional packed starting vector; otherwise threshold starts and tau=1 are used.
      integer, intent(in), optional :: max_iter !! Maximum BFGS iterations; defaults to 250.
      real(dp), intent(in), optional :: grad_tol !! Infinity-norm numerical-gradient tolerance; defaults to 1e-6.
      real(dp), allocatable :: par(:), alpha0(:)
      real(dp) :: fval, tol
      integer :: nalpha, nbeta, npar, st, invst, limit, diagst
      nalpha = threshold_parameter_count(problem%nclass, problem%threshold)
      nbeta = size(problem%x, 2)
      npar = nalpha + nbeta + 1
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
         par(:nalpha) = alpha0
         par(npar) = log(1.0_dp)
      end if
      limit = 250
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-6_dp
      if (present(grad_tol)) tol = grad_tol
      call bfgs_minimize(problem, par, fval, fit%iterations, fit%status, limit, tol)
      fit%converged = fit%status == 0
      fit%loglik = -fval
      fit%tau = exp(par(npar))
      fit%nAGQ = problem%nAGQ
      allocate(fit%alpha(nalpha), fit%beta(nbeta), fit%gradient(npar), fit%hessian(npar, npar), fit%vcov(npar, npar), &
               fit%ranef_mode(problem%ngroup), fit%ranef_condvar(problem%ngroup))
      fit%alpha = par(:nalpha)
      if (nbeta > 0) fit%beta = par(nalpha + 1:nalpha + nbeta)
      call thresholds_from_alpha(fit%alpha, problem%nclass, problem%threshold, fit%theta, st)
      call numerical_gradient(problem, par, fit%gradient, 1.0e-5_dp)
      call numerical_hessian(problem, par, fit%hessian, 5.0e-4_dp)
      fit%max_gradient = maxval(abs(fit%gradient))
      call invert_matrix(fit%hessian, fit%vcov, invst)
      if (invst /= 0) fit%vcov = 0.0_dp
      call hessian_diagnostics(fit%hessian, fit%hessian_min_eigenvalue, fit%hessian_max_eigenvalue, &
           fit%hessian_condition, fit%hessian_rank, fit%hessian_positive_definite, diagst)
      if (diagst /= 0) fit%hessian_positive_definite = .false.
      call clmm_conditional_modes(problem, fit%theta, fit%beta, fit%tau, fit%ranef_mode, fit%ranef_condvar)
   end subroutine fit_clmm

   subroutine clmm_conditional_modes(problem, theta, beta, tau, modes, condvar)
      type(clmm_problem), intent(in) :: problem !! Random-intercept cumulative-link data and model specification.
      real(dp), intent(in) :: theta(:) !! Fitted increasing cut points.
      real(dp), intent(in) :: beta(:) !! Fitted fixed-effect location coefficients.
      real(dp), intent(in) :: tau !! Positive Gaussian random-intercept standard deviation.
      real(dp), intent(out) :: modes(:) !! Posterior modes of random intercepts, one per internally remapped group.
      real(dp), intent(out), optional :: condvar(:) !! Approximate conditional variances, inverse penalized curvature by group.
      real(dp) :: curvature, qmode
      integer :: g, st
      if (size(modes) /= problem%ngroup) return
      if (present(condvar)) then
         if (size(condvar) /= problem%ngroup) return
      end if
      do g = 1, problem%ngroup
         call group_mode_curvature(problem, theta, beta, tau, g, modes(g), curvature, qmode, st)
         if (st /= 0) then
            modes(g) = 0.0_dp
            if (present(condvar)) condvar(g) = huge(1.0_dp)
         else
            if (present(condvar)) condvar(g) = 1.0_dp/curvature
         end if
      end do
   end subroutine clmm_conditional_modes
end module ordinal_clmm
