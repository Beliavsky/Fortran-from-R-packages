! General Laplace cumulative-link mixed-model kernels based on ordinal/R/clmm.R.
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_clmm_general
   use ordinal_kinds, only : dp
   use ordinal_links, only : link_cdf, link_pdf, link_pdf_gradient, link_logit
   use ordinal_thresholds, only : threshold_flexible, threshold_parameter_count, threshold_start, thresholds_from_alpha
   use ordinal_numerics, only : objective_type, bfgs_minimize, numerical_gradient, numerical_hessian, invert_matrix, &
                                cholesky_factor, solve_cholesky, logdet_spd, symmetric_eigenvalues, hessian_diagnostics
   implicit none
   private
   real(dp), parameter :: bad_objective = huge(1.0_dp)/1000.0_dp
   type, extends(objective_type), public :: clmm_laplace_problem
      integer, allocatable :: y(:)
      integer, allocatable :: re_group(:, :)
      integer, allocatable :: term_q(:)
      integer, allocatable :: term_ngroup(:)
      integer, allocatable :: term_col_start(:)
      integer, allocatable :: term_re_start(:)
      integer, allocatable :: term_cov_start(:)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: re_z(:, :)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: offset(:)
      integer :: nclass = 0
      integer :: nterms = 0
      integer :: nrandom = 0
      integer :: ncovpar = 0
      integer :: link = link_logit
      integer :: threshold = threshold_flexible
      real(dp) :: lambda = 1.0_dp
   contains
      procedure :: value => clmm_laplace_value
   end type clmm_laplace_problem
   type, public :: clmm_laplace_fit_result
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: covariance_parameters(:)
      real(dp), allocatable :: re_cov(:, :)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: vcov(:, :)
      real(dp), allocatable :: ranef_mode(:)
      real(dp), allocatable :: ranef_condvar(:, :)
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
   end type clmm_laplace_fit_result
   public :: init_clmm_laplace_problem, clmm_laplace_nll, fit_clmm_laplace, clmm_laplace_modes
contains
   subroutine init_clmm_laplace_problem(problem, y, x, re_group, re_z, term_q, nclass, link, threshold, weights, &
                                        offset, lambda, status)
      type(clmm_laplace_problem), intent(out) :: problem !! Initialized dense general-Laplace cumulative-link mixed model.
      integer, intent(in) :: y(:) !! Ordered response codes in 1:nclass, one per observation.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect location design matrix without an intercept.
      integer, intent(in) :: re_group(:, :) !! Positive grouping labels; columns correspond to random-effect terms.
      real(dp), intent(in) :: re_z(:, :) !! Random-effect covariates; term column blocks are concatenated left to right.
      integer, intent(in) :: term_q(:) !! Number of random coefficients per group level for each random-effect term.
      integer, intent(in) :: nclass !! Number of ordered response categories; must be at least two.
      integer, intent(in), optional :: link !! Link identifier; defaults to logistic.
      integer, intent(in), optional :: threshold !! Threshold structure identifier; defaults to flexible.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights; defaults to one.
      real(dp), intent(in), optional :: offset(:) !! Fixed location offsets, one per observation; defaults to zero.
      real(dp), intent(in), optional :: lambda !! Fixed flexible-link shape parameter; defaults to one.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, group labels, or options.
      integer, allocatable :: labels(:)
      integer :: n, nt, t, i, g, ng, qsum, re_start, col_start, cov_start
      n = size(y)
      nt = size(term_q)
      status = 0
      if (nclass < 2 .or. nt < 1 .or. size(x, 1) /= n .or. size(re_group, 1) /= n .or. &
          size(re_group, 2) /= nt .or. size(re_z, 1) /= n .or. any(term_q < 1) .or. &
          any(y < 1) .or. any(y > nclass) .or. any(re_group <= 0)) then
         status = 1
         return
      end if
      qsum = sum(term_q)
      if (size(re_z, 2) /= qsum) then
         status = 2
         return
      end if
      problem%nclass = nclass
      problem%nterms = nt
      if (present(link)) problem%link = link
      if (present(threshold)) problem%threshold = threshold
      if (threshold_parameter_count(nclass, problem%threshold) <= 0) then
         status = 3
         return
      end if
      if (present(lambda)) problem%lambda = lambda
      problem%y = y
      problem%x = x
      problem%re_z = re_z
      problem%term_q = term_q
      allocate(problem%re_group(n, nt), problem%term_ngroup(nt), problem%term_col_start(nt), &
               problem%term_re_start(nt), problem%term_cov_start(nt), problem%weights(n), problem%offset(n), labels(n))
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
      col_start = 1
      re_start = 1
      cov_start = 1
      problem%nrandom = 0
      problem%ncovpar = 0
      do t = 1, nt
         problem%term_col_start(t) = col_start
         problem%term_re_start(t) = re_start
         problem%term_cov_start(t) = cov_start
         labels = 0
         ng = 0
         do i = 1, n
            g = find_label(re_group(i, t), labels, ng)
            if (g == 0) then
               ng = ng + 1
               labels(ng) = re_group(i, t)
               g = ng
            end if
            problem%re_group(i, t) = g
         end do
         problem%term_ngroup(t) = ng
         problem%nrandom = problem%nrandom + ng*term_q(t)
         problem%ncovpar = problem%ncovpar + term_q(t)*(term_q(t) + 1)/2
         col_start = col_start + term_q(t)
         re_start = re_start + ng*term_q(t)
         cov_start = cov_start + term_q(t)*(term_q(t) + 1)/2
      end do
   contains
      pure integer function find_label(value, values, nused) result(pos)
         integer, intent(in) :: value !! Original positive grouping label to locate.
         integer, intent(in) :: values(:) !! Previously encountered original labels for one random-effect term.
         integer, intent(in) :: nused !! Number of valid entries at the front of values.
         integer :: j
         pos = 0
         do j = 1, nused
            if (values(j) == value) then
               pos = j
               return
            end if
         end do
      end function find_label
   end subroutine init_clmm_laplace_problem

   function clmm_laplace_nll(par, problem) result(nll)
      real(dp), intent(in) :: par(:) !! Packed thresholds, fixed effects, and random covariance-factor parameters.
      type(clmm_laplace_problem), intent(in) :: problem !! General cumulative-link mixed-model data and structure.
      real(dp) :: nll
      nll = clmm_laplace_nll_problem(par, problem)
   end function clmm_laplace_nll

   function clmm_laplace_value(self, x) result(nll)
      class(clmm_laplace_problem), intent(in) :: self !! General mixed-model objective and associated data.
      real(dp), intent(in) :: x(:) !! Packed thresholds, fixed effects, and covariance-factor parameters.
      real(dp) :: nll
      nll = clmm_laplace_nll_problem(x, self)
   end function clmm_laplace_value

   function clmm_laplace_nll_problem(par, problem) result(nll)
      real(dp), intent(in) :: par(:) !! Packed general-Laplace mixed-model parameter vector.
      type(clmm_laplace_problem), intent(in) :: problem !! General cumulative-link mixed-model data and structure.
      real(dp) :: nll
      real(dp), allocatable :: theta(:), beta(:), covpar(:), mode(:), mode_hessian(:, :), covsmall(:, :), invsmall(:, :)
      real(dp) :: qmode, logdet_g, logdet_h
      integer :: nalpha, nbeta, npar, st
      nalpha = threshold_parameter_count(problem%nclass, problem%threshold)
      nbeta = size(problem%x, 2)
      npar = nalpha + nbeta + problem%ncovpar
      if (size(par) /= npar) then
         nll = bad_objective
         return
      end if
      call thresholds_from_alpha(par(:nalpha), problem%nclass, problem%threshold, theta, st)
      if (st /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) then
         nll = bad_objective
         return
      end if
      allocate(beta(nbeta), covpar(problem%ncovpar), mode(problem%nrandom), mode_hessian(problem%nrandom, problem%nrandom))
      if (nbeta > 0) beta = par(nalpha + 1:nalpha + nbeta)
      covpar = par(nalpha + nbeta + 1:)
      call build_covariance_info(problem, covpar, covsmall, invsmall, logdet_g, st)
      if (st /= 0) then
         nll = bad_objective
         return
      end if
      call find_general_mode(problem, theta, beta, invsmall, mode, qmode, mode_hessian, st)
      if (st /= 0) then
         nll = bad_objective
         return
      end if
      call logdet_spd(mode_hessian, logdet_h, st)
      if (st /= 0) then
         nll = bad_objective
         return
      end if
      nll = qmode + 0.5_dp*logdet_g + 0.5_dp*logdet_h
   end function clmm_laplace_nll_problem

   subroutine build_covariance_info(problem, covpar, covariance, inverse_covariance, logdet_g, status)
      type(clmm_laplace_problem), intent(in) :: problem !! Random-effect term dimensions and covariance parameter layout.
      real(dp), intent(in) :: covpar(:) !! Packed lower-Cholesky parameters with log diagonals for every random term.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Block-diagonal term covariance matrices over concatenated q columns.
      real(dp), allocatable, intent(out) :: inverse_covariance(:, :) !! Block-diagonal inverses matching covariance.
      real(dp), intent(out) :: logdet_g !! Log determinant of the full replicated random-effects covariance matrix.
      integer, intent(out) :: status !! Zero on success; nonzero for parameter dimensions or numerical inversion failure.
      real(dp), allocatable :: l(:, :), gmat(:, :), ginv(:, :)
      integer :: qsum, t, q, i, j, idx, c0, p0, invst
      qsum = sum(problem%term_q)
      allocate(covariance(qsum, qsum), inverse_covariance(qsum, qsum))
      covariance = 0.0_dp
      inverse_covariance = 0.0_dp
      logdet_g = 0.0_dp
      status = 0
      if (size(covpar) /= problem%ncovpar) then
         status = 1
         return
      end if
      do t = 1, problem%nterms
         q = problem%term_q(t)
         c0 = problem%term_col_start(t)
         p0 = problem%term_cov_start(t)
         allocate(l(q, q), gmat(q, q), ginv(q, q))
         l = 0.0_dp
         idx = p0
         do i = 1, q
            do j = 1, i
               if (i == j) then
                  if (covpar(idx) < log(1.0e-8_dp) .or. covpar(idx) > log(1.0e4_dp)) then
                     status = 2
                     return
                  end if
                  l(i, j) = exp(covpar(idx))
               else
                  if (abs(covpar(idx)) > 1.0e4_dp) then
                     status = 2
                     return
                  end if
                  l(i, j) = covpar(idx)
               end if
               idx = idx + 1
            end do
         end do
         gmat = matmul(l, transpose(l))
         call invert_matrix(gmat, ginv, invst)
         if (invst /= 0) then
            status = 3
            return
         end if
         covariance(c0:c0 + q - 1, c0:c0 + q - 1) = gmat
         inverse_covariance(c0:c0 + q - 1, c0:c0 + q - 1) = ginv
         do i = 1, q
            logdet_g = logdet_g + 2.0_dp*real(problem%term_ngroup(t), dp)*log(l(i, i))
         end do
         deallocate(l, gmat, ginv)
      end do
   end subroutine build_covariance_info

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

   subroutine general_penalized_derivatives(problem, theta, beta, inverse_covariance, b, qvalue, gradient, hessian, status)
      type(clmm_laplace_problem), intent(in) :: problem !! General cumulative-link mixed-model data and random structure.
      real(dp), intent(in) :: theta(:) !! Candidate increasing cut points.
      real(dp), intent(in) :: beta(:) !! Candidate fixed-effect location coefficients.
      real(dp), intent(in) :: inverse_covariance(:, :) !! Block-diagonal inverse covariance for each random term's q coefficients.
      real(dp), intent(in) :: b(:) !! Candidate direct random coefficients over all term levels.
      real(dp), intent(out) :: qvalue !! Penalized negative conditional log-likelihood excluding Gaussian constants.
      real(dp), intent(out) :: gradient(:) !! Gradient of qvalue with respect to all random coefficients.
      real(dp), intent(out) :: hessian(:, :) !! Hessian of qvalue with respect to all random coefficients.
      integer, intent(out) :: status !! Zero on success; nonzero for dimensions or invalid category probabilities.
      integer, allocatable :: active_index(:)
      real(dp), allocatable :: active_value(:), bv(:), gv(:)
      real(dp) :: eta, prob, dpeta, d2peta, grad_eta, hess_eta, w
      integer :: t, q, g, i, j, k, c0, r0, nactive, idx
      qvalue = 0.0_dp
      gradient = 0.0_dp
      hessian = 0.0_dp
      status = 0
      if (size(b) /= problem%nrandom .or. size(gradient) /= problem%nrandom .or. &
          size(hessian, 1) /= problem%nrandom .or. size(hessian, 2) /= problem%nrandom) then
         status = 1
         return
      end if
      do t = 1, problem%nterms
         q = problem%term_q(t)
         c0 = problem%term_col_start(t)
         allocate(bv(q), gv(q))
         do g = 1, problem%term_ngroup(t)
            r0 = problem%term_re_start(t) + (g - 1)*q
            bv = b(r0:r0 + q - 1)
            gv = matmul(inverse_covariance(c0:c0 + q - 1, c0:c0 + q - 1), bv)
            qvalue = qvalue + 0.5_dp*dot_product(bv, gv)
            gradient(r0:r0 + q - 1) = gradient(r0:r0 + q - 1) + gv
            hessian(r0:r0 + q - 1, r0:r0 + q - 1) = hessian(r0:r0 + q - 1, r0:r0 + q - 1) &
                 + inverse_covariance(c0:c0 + q - 1, c0:c0 + q - 1)
         end do
         deallocate(bv, gv)
      end do
      allocate(active_index(sum(problem%term_q)), active_value(sum(problem%term_q)))
      do i = 1, size(problem%y)
         if (problem%weights(i) <= 0.0_dp) cycle
         eta = problem%offset(i)
         if (size(beta) > 0) eta = eta + dot_product(problem%x(i, :), beta)
         nactive = 0
         do t = 1, problem%nterms
            q = problem%term_q(t)
            c0 = problem%term_col_start(t)
            g = problem%re_group(i, t)
            r0 = problem%term_re_start(t) + (g - 1)*q
            eta = eta + dot_product(problem%re_z(i, c0:c0 + q - 1), b(r0:r0 + q - 1))
            do j = 1, q
               nactive = nactive + 1
               active_index(nactive) = r0 + j - 1
               active_value(nactive) = problem%re_z(i, c0 + j - 1)
            end do
         end do
         call category_probability_derivatives(problem%y(i), theta, eta, problem%link, problem%lambda, prob, dpeta, d2peta)
         if (prob <= tiny(1.0_dp) .or. prob > 1.0_dp) then
            status = 2
            qvalue = bad_objective
            return
         end if
         w = problem%weights(i)
         qvalue = qvalue - w*log(prob)
         grad_eta = -w*dpeta/prob
         hess_eta = w*((dpeta/prob)**2 - d2peta/prob)
         do j = 1, nactive
            idx = active_index(j)
            gradient(idx) = gradient(idx) + grad_eta*active_value(j)
         end do
         do j = 1, nactive
            do k = 1, nactive
               hessian(active_index(j), active_index(k)) = hessian(active_index(j), active_index(k)) &
                    + hess_eta*active_value(j)*active_value(k)
            end do
         end do
      end do
   end subroutine general_penalized_derivatives

   subroutine find_general_mode(problem, theta, beta, inverse_covariance, mode, qmode, mode_hessian, status)
      type(clmm_laplace_problem), intent(in) :: problem !! General cumulative-link mixed-model data and random structure.
      real(dp), intent(in) :: theta(:) !! Candidate increasing cut points.
      real(dp), intent(in) :: beta(:) !! Candidate fixed-effect location coefficients.
      real(dp), intent(in) :: inverse_covariance(:, :) !! Per-term inverse covariance blocks over concatenated q columns.
      real(dp), intent(out) :: mode(:) !! Conditional modes of all direct random coefficients.
      real(dp), intent(out) :: qmode !! Penalized negative conditional log-likelihood at mode.
      real(dp), intent(out) :: mode_hessian(:, :) !! Random-effect Hessian at the returned conditional mode.
      integer, intent(out) :: status !! Zero on convergence; nonzero for derivative, factorization, or line-search failure.
      real(dp), allocatable :: grad(:), hwork(:, :), chol(:, :), step(:), trial(:), eig(:), trial_hess(:, :), trial_grad(:)
      real(dp) :: trial_q, alpha_step, mineig, inflation, scale
      integer :: n, iter, ls, st, chst, solst, i
      n = problem%nrandom
      if (size(mode) /= n .or. size(mode_hessian, 1) /= n .or. size(mode_hessian, 2) /= n) then
         status = 1
         return
      end if
      allocate(grad(n), hwork(n, n), chol(n, n), step(n), trial(n), eig(n), trial_hess(n, n), trial_grad(n))
      mode = 0.0_dp
      status = 0
      do iter = 1, 100
         call general_penalized_derivatives(problem, theta, beta, inverse_covariance, mode, qmode, grad, mode_hessian, st)
         if (st /= 0) then
            status = 2
            return
         end if
         if (maxval(abs(grad)) <= 1.0e-8_dp*(1.0_dp + abs(qmode))) exit
         hwork = mode_hessian
         call cholesky_factor(hwork, chol, chst)
         if (chst /= 0) then
            call symmetric_eigenvalues(mode_hessian, eig, st)
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
            trial = mode - alpha_step*step
            call general_penalized_derivatives(problem, theta, beta, inverse_covariance, trial, trial_q, trial_grad, &
                                               trial_hess, st)
            if (st == 0 .and. trial_q < qmode) exit
            alpha_step = 0.5_dp*alpha_step
         end do
         if (ls > 40) then
            if (maxval(abs(grad)) <= 1.0e-6_dp) exit
            status = 4
            return
         end if
         mode = trial
         if (maxval(abs(alpha_step*step)) <= 1.0e-10_dp*(1.0_dp + maxval(abs(mode)))) exit
      end do
      call general_penalized_derivatives(problem, theta, beta, inverse_covariance, mode, qmode, grad, mode_hessian, st)
      if (st /= 0) status = 5
   end subroutine find_general_mode

   subroutine clmm_laplace_modes(problem, alpha, beta, covariance_parameters, modes, condvar, status)
      type(clmm_laplace_problem), intent(in) :: problem !! General cumulative-link mixed-model data and random structure.
      real(dp), intent(in) :: alpha(:) !! Free threshold parameters in the selected threshold structure.
      real(dp), intent(in) :: beta(:) !! Fixed-effect location coefficients.
      real(dp), intent(in) :: covariance_parameters(:) !! Packed random-effect Cholesky parameters with log diagonals.
      real(dp), intent(out) :: modes(:) !! Conditional modes of all direct random coefficients.
      real(dp), intent(out) :: condvar(:, :) !! Approximate conditional covariance, inverse random-effect Hessian at the modes.
      integer, intent(out) :: status !! Zero on success; nonzero for parameter dimensions or numerical failure.
      real(dp), allocatable :: theta(:), covariance(:, :), inverse_covariance(:, :), mode_hessian(:, :)
      real(dp) :: qmode, logdet_g
      integer :: st, invst
      status = 0
      if (size(alpha) /= threshold_parameter_count(problem%nclass, problem%threshold) .or. &
          size(beta) /= size(problem%x, 2) .or. size(covariance_parameters) /= problem%ncovpar .or. &
          size(modes) /= problem%nrandom .or. size(condvar, 1) /= problem%nrandom .or. &
          size(condvar, 2) /= problem%nrandom) then
         status = 1
         return
      end if
      call thresholds_from_alpha(alpha, problem%nclass, problem%threshold, theta, st)
      if (st /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) then
         status = 2
         return
      end if
      call build_covariance_info(problem, covariance_parameters, covariance, inverse_covariance, logdet_g, st)
      if (st /= 0) then
         status = 3
         return
      end if
      allocate(mode_hessian(problem%nrandom, problem%nrandom))
      call find_general_mode(problem, theta, beta, inverse_covariance, modes, qmode, mode_hessian, st)
      if (st /= 0) then
         status = 4
         return
      end if
      call invert_matrix(mode_hessian, condvar, invst)
      if (invst /= 0) then
         status = 5
         condvar = 0.0_dp
      end if
   end subroutine clmm_laplace_modes

   subroutine fit_clmm_laplace(problem, fit, start, max_iter, grad_tol)
      type(clmm_laplace_problem), intent(in) :: problem !! Fully initialized general-Laplace cumulative-link mixed model.
      type(clmm_laplace_fit_result), intent(out) :: fit !! Fitted parameters, modes, covariance estimates, and diagnostics.
      real(dp), intent(in), optional :: start(:) !! Optional packed start; defaults use threshold starts and identity covariance.
      integer, intent(in), optional :: max_iter !! Maximum BFGS iterations; defaults to 200.
      real(dp), intent(in), optional :: grad_tol !! Infinity-norm numerical-gradient tolerance; defaults to 2e-6.
      real(dp), allocatable :: par(:), alpha0(:), covariance(:, :), inverse_covariance(:, :)
      real(dp) :: fval, tol, logdet_g
      integer :: nalpha, nbeta, npar, limit, st, invst, diagst
      nalpha = threshold_parameter_count(problem%nclass, problem%threshold)
      nbeta = size(problem%x, 2)
      npar = nalpha + nbeta + problem%ncovpar
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
      end if
      limit = 200
      if (present(max_iter)) limit = max_iter
      tol = 2.0e-6_dp
      if (present(grad_tol)) tol = grad_tol
      call bfgs_minimize(problem, par, fval, fit%iterations, fit%status, limit, tol)
      fit%converged = fit%status == 0
      fit%loglik = -fval
      allocate(fit%alpha(nalpha), fit%beta(nbeta), fit%covariance_parameters(problem%ncovpar), fit%gradient(npar), &
               fit%hessian(npar, npar), fit%vcov(npar, npar), fit%ranef_mode(problem%nrandom), &
               fit%ranef_condvar(problem%nrandom, problem%nrandom))
      fit%alpha = par(:nalpha)
      if (nbeta > 0) fit%beta = par(nalpha + 1:nalpha + nbeta)
      fit%covariance_parameters = par(nalpha + nbeta + 1:)
      call thresholds_from_alpha(fit%alpha, problem%nclass, problem%threshold, fit%theta, st)
      call build_covariance_info(problem, fit%covariance_parameters, covariance, inverse_covariance, logdet_g, st)
      if (st == 0) fit%re_cov = covariance
      call numerical_gradient(problem, par, fit%gradient, 2.0e-5_dp)
      call numerical_hessian(problem, par, fit%hessian, 8.0e-4_dp)
      fit%max_gradient = maxval(abs(fit%gradient))
      call invert_matrix(fit%hessian, fit%vcov, invst)
      if (invst /= 0) fit%vcov = 0.0_dp
      call hessian_diagnostics(fit%hessian, fit%hessian_min_eigenvalue, fit%hessian_max_eigenvalue, &
           fit%hessian_condition, fit%hessian_rank, fit%hessian_positive_definite, diagst)
      if (diagst /= 0) fit%hessian_positive_definite = .false.
      call clmm_laplace_modes(problem, fit%alpha, fit%beta, fit%covariance_parameters, fit%ranef_mode, &
                              fit%ranef_condvar, st)
      if (st /= 0 .and. fit%status == 0) then
         fit%status = 20 + st
         fit%converged = .false.
      end if
   end subroutine fit_clmm_laplace
end module ordinal_clmm_general
