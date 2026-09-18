! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_refit
   use flexmix_kinds, only : dp
   use flexmix_types
   use flexmix_components, only : gaussian_regression_log_density, poisson_regression_log_density
   use flexmix_components, only : binomial_regression_log_density, gamma_regression_log_density
   use flexmix_numeric, only : normalize_log_probabilities
   implicit none
   private

   public :: flexmix_get_design
   public :: flexmix_parameter_vector
   public :: flexmix_replace_parameters
   public :: flexmix_exist_gradient
   public :: flexmix_loglik_regression_at
   public :: flexmix_gradient_regression
   public :: flexmix_refit_optim_regression

contains

   subroutine flexmix_get_design(result, design, info)
      type(flexmix_result), intent(in) :: result !! Fitted model whose component-to-parameter incidence matrix is requested.
      logical, allocatable, intent(out) :: design(:,:) !! Component incidence matrix; rows are retained components and columns.
      integer, intent(out) :: info !! Zero on success; nonzero when the fitted family is unsupported or structurally incomplete.
      integer, allocatable :: counts(:)
      integer :: j, first, last

      if (allocated(result%parameter_design)) then
         call fixed_parameter_design(result, design, info)
         return
      end if
      call component_parameter_counts(result, counts, info)
      if (info /= 0) then
         allocate(design(0,0))
         return
      end if
      allocate(design(result%k,sum(counts)))
      design = .false.
      first = 1
      do j = 1, result%k
         last = first + counts(j) - 1
         if (counts(j) > 0) design(j,first:last) = .true.
         first = last + 1
      end do
   end subroutine flexmix_get_design

   subroutine flexmix_parameter_vector(result, parameters, info)
      type(flexmix_result), intent(in) :: result !! Fitted model supplying component/concomitant parameters for refitting.
      real(dp), allocatable, intent(out) :: parameters(:) !! Unconstrained component blocks followed by prior coefficients.
      integer, intent(out) :: info !! Zero on success; nonzero when required family-specific arrays are absent or unsupported.
      integer, allocatable :: counts(:)
      integer :: component_total, j, pos, q

      if (allocated(result%parameter_design)) then
         call pack_fixed_parameters(result, parameters, info)
         return
      end if
      call component_parameter_counts(result, counts, info)
      if (info /= 0) then
         allocate(parameters(0))
         return
      end if
      component_total = sum(counts)
      q = prior_parameter_count(result)
      allocate(parameters(component_total + q))
      pos = 1
      do j = 1, result%k
         call pack_component_parameters(result, j, parameters, pos, info)
         if (info /= 0) return
      end do
      call pack_prior_parameters(result, parameters, pos, info)
   end subroutine flexmix_parameter_vector

   subroutine flexmix_replace_parameters(result, parameters, updated, info)
      type(flexmix_result), intent(in) :: result !! Template defining family, dimensions, component roles, and parameter layout.
      real(dp), intent(in) :: parameters(:) !! Replacement parameter vector in the `flexmix_parameter_vector` layout.
      type(flexmix_result), intent(out) :: updated !! Copy of `result` with numerical component and concomitant parameters.
      integer, intent(out) :: info !! Zero on success; nonzero for a length mismatch, unsupported family, or malformed stored.
      integer, allocatable :: counts(:)
      integer :: expected, j, pos

      if (allocated(result%parameter_design)) then
         call unpack_fixed_parameters(result, parameters, updated, info)
         return
      end if
      call component_parameter_counts(result, counts, info)
      if (info /= 0) then
         updated = result
         return
      end if
      expected = sum(counts) + prior_parameter_count(result)
      if (size(parameters) /= expected) then
         updated = result
         info = 2
         return
      end if
      updated = result
      pos = 1
      do j = 1, result%k
         call unpack_component_parameters(updated, j, parameters, pos, info)
         if (info /= 0) return
      end do
      call unpack_prior_parameters(updated, parameters, pos, info)
   end subroutine flexmix_replace_parameters

   pure logical function flexmix_exist_gradient(result) result(available)
      type(flexmix_result), intent(in) :: result !! Fitted model whose analytic observed-data score availability is queried.
      if (allocated(result%parameter_design)) then
         available = .false.
         return
      end if
      select case (result%model_kind)
      case (flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial)
         available = .true.
      case default
         available = .false.
      end select
   end function flexmix_exist_gradient

   subroutine flexmix_loglik_regression_at(template, x, y, parameters, value, info, trials, offset, case_weights, concomitant_x)
      type(flexmix_result), intent(in) :: template !! Fitted regression-mixture template defining component/prior parameter layout.
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n,p)`, used to evaluate each component likelihood.
      real(dp), intent(in) :: y(:) !! Response vector, size `n`; for binomial models this contains successes.
      real(dp), intent(in) :: parameters(:) !! Candidate unconstrained parameter vector returned by `flexmix_parameter_vector`.
      real(dp), intent(out) :: value !! Observed-data log likelihood at `parameters` for the supplied rows.
      integer, intent(out) :: info !! Zero on success; nonzero for dimensions, unsupported family, or parameter layout.
      real(dp), intent(in), optional :: trials(:) !! Optional binomial trial counts, size `n`; omitted values imply Bernoulli.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case weights, size `n`; defaults to one.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant design matrix, shape `(n,q)`, matching stored.
      type(flexmix_result) :: work
      real(dp), allocatable :: log_density(:,:), prior_rows(:,:), logp(:,:), posterior(:,:), logsum(:), weights(:), total(:)
      integer :: j, n

      value = -huge(1.0_dp)
      n = size(y)
      if (size(x,1) /= n .or. size(x,2) /= template%p) then
         info = 10
         return
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            info = 11
            return
         end if
      end if
      if (present(trials)) then
         if (size(trials) /= n) then
            info = 12
            return
         end if
      end if
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            info = 13
            return
         end if
      end if
      call flexmix_replace_parameters(template, parameters, work, info)
      if (info /= 0) return
      allocate(log_density(n,work%k), prior_rows(n,work%k), logp(n,work%k), posterior(n,work%k), logsum(n), weights(n), total(n))
      weights = 1.0_dp
      if (present(case_weights)) weights = case_weights
      do j = 1, work%k
         select case (work%model_kind)
         case (flexmix_model_gaussian)
            call gaussian_regression_log_density(x, y, work%beta(:,j), work%sigma(j), log_density(:,j), offset)
         case (flexmix_model_poisson)
            call poisson_regression_log_density(x, y, work%beta(:,j), log_density(:,j), offset)
         case (flexmix_model_binomial)
            if (present(trials)) then
               call binomial_regression_log_density(x, y, trials, work%beta(:,j), log_density(:,j), offset)
            else
               total = 1.0_dp
               call binomial_regression_log_density(x, y, total, work%beta(:,j), log_density(:,j), offset)
            end if
         case (flexmix_model_gamma_regression)
            call gamma_regression_log_density(x, y, work%beta(:,j), work%shape(j), log_density(:,j), offset)
         case default
            info = 14
            return
         end select
      end do
      call prior_matrix(work, n, prior_rows, info, concomitant_x)
      if (info /= 0) return
      logp = log_density + log(max(prior_rows,tiny(1.0_dp)))
      call normalize_log_probabilities(logp, posterior, logsum)
      value = sum(weights * logsum)
      info = 0
   end subroutine flexmix_loglik_regression_at

   subroutine flexmix_gradient_regression(template, x, y, parameters, gradient, info, trials, offset, case_weights, concomitant_x)
      type(flexmix_result), intent(in) :: template !! Fitted Gaussian, Poisson, or binomial mixture template defining score.
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n,p)`, used for component scores.
      real(dp), intent(in) :: y(:) !! Response vector, size `n`; for binomial models this contains successes.
      real(dp), intent(in) :: parameters(:) !! Candidate unconstrained parameter vector returned by `flexmix_parameter_vector`.
      real(dp), allocatable, intent(out) :: gradient(:) !! Analytic observed-data score in the same order as `parameters`.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported models, invalid dimensions, or malformed.
      real(dp), intent(in), optional :: trials(:) !! Optional binomial trial counts, size `n`; omitted values imply Bernoulli.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant design matrix, shape `(n,q)`, matching stored.
      type(flexmix_result) :: work
      real(dp), allocatable :: log_density(:,:), prior_rows(:,:), logp(:,:), posterior(:,:), logsum(:), weights(:)
      real(dp), allocatable :: eta(:), mu(:), total(:), residual(:)
      integer, allocatable :: counts(:)
      integer :: j, n, p, pos, q, c

      if (.not. flexmix_exist_gradient(template)) then
         allocate(gradient(0))
         info = 20
         return
      end if
      n = size(y)
      p = size(x,2)
      if (size(x,1) /= n .or. p /= template%p) then
         allocate(gradient(0))
         info = 21
         return
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            allocate(gradient(0))
            info = 22
            return
         end if
      end if
      if (present(trials)) then
         if (size(trials) /= n) then
            allocate(gradient(0))
            info = 23
            return
         end if
      end if
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            allocate(gradient(0))
            info = 24
            return
         end if
      end if
      call flexmix_replace_parameters(template, parameters, work, info)
      if (info /= 0) then
         allocate(gradient(0))
         return
      end if
      call component_parameter_counts(work, counts, info)
      if (info /= 0) then
         allocate(gradient(0))
         return
      end if
      allocate(gradient(size(parameters)))
      gradient = 0.0_dp
      allocate(log_density(n,work%k), prior_rows(n,work%k), logp(n,work%k), posterior(n,work%k), logsum(n), weights(n))
      allocate(eta(n), mu(n), total(n), residual(n))
      weights = 1.0_dp
      if (present(case_weights)) weights = case_weights
      if (present(trials)) then
         total = trials
      else
         total = 1.0_dp
      end if
      do j = 1, work%k
         select case (work%model_kind)
         case (flexmix_model_gaussian)
            call gaussian_regression_log_density(x, y, work%beta(:,j), work%sigma(j), log_density(:,j), offset)
         case (flexmix_model_poisson)
            call poisson_regression_log_density(x, y, work%beta(:,j), log_density(:,j), offset)
         case (flexmix_model_binomial)
            call binomial_regression_log_density(x, y, total, work%beta(:,j), log_density(:,j), offset)
         end select
      end do
      call prior_matrix(work, n, prior_rows, info, concomitant_x)
      if (info /= 0) return
      logp = log_density + log(max(prior_rows,tiny(1.0_dp)))
      call normalize_log_probabilities(logp, posterior, logsum)

      pos = 1
      do j = 1, work%k
         if (counts(j) == 0) cycle
         eta = matmul(x, work%beta(:,j))
         if (present(offset)) eta = eta + offset
         select case (work%model_kind)
         case (flexmix_model_gaussian)
            residual = y - eta
            do c = 1, p
               gradient(pos+c-1) = sum(weights * posterior(:,j) * residual * x(:,c)) / (work%sigma(j)**2)
            end do
            gradient(pos+p) = sum(weights * posterior(:,j) * (-1.0_dp + (residual/work%sigma(j))**2))
         case (flexmix_model_poisson)
            mu = exp(max(-30.0_dp,min(30.0_dp,eta)))
            residual = y - mu
            do c = 1, p
               gradient(pos+c-1) = sum(weights * posterior(:,j) * residual * x(:,c))
            end do
         case (flexmix_model_binomial)
            mu = 1.0_dp / (1.0_dp + exp(-max(-30.0_dp,min(30.0_dp,eta))))
            residual = y - total * mu
            do c = 1, p
               gradient(pos+c-1) = sum(weights * posterior(:,j) * residual * x(:,c))
            end do
         end select
         pos = pos + counts(j)
      end do

      if (work%k > 1) then
         if (allocated(work%concomitant_coef)) then
            q = size(work%concomitant_coef,1)
            if (.not. present(concomitant_x)) then
               info = 25
               return
            end if
            if (size(concomitant_x,1) /= n .or. size(concomitant_x,2) /= q) then
               info = 26
               return
            end if
            do j = 2, work%k
               do c = 1, q
                  gradient(pos) = sum(weights * (posterior(:,j) - prior_rows(:,j)) * concomitant_x(:,c))
                  pos = pos + 1
               end do
            end do
         else
            do j = 2, work%k
               gradient(pos) = sum(weights * (posterior(:,j) - prior_rows(:,j)))
               pos = pos + 1
            end do
         end if
      end if
      info = 0
   end subroutine flexmix_gradient_regression

   subroutine flexmix_refit_optim_regression(result, x, y, estimate, standard_error, covariance, info, trials, offset, &
                                               case_weights, concomitant_x, difference_step)
      type(flexmix_result), intent(in) :: result !! Fitted Gaussian, Poisson, or binomial mixture evaluated at its current.
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n,p)`, used for score and Hessian evaluation.
      real(dp), intent(in) :: y(:) !! Response vector, size `n`; for binomial models this contains successes.
      real(dp), allocatable, intent(out) :: estimate(:) !! Packed current parameter estimates in `flexmix_parameter_vector`
      real(dp), allocatable, intent(out) :: standard_error(:) !! Standard errors from the inverse observed information, same.
      real(dp), allocatable, intent(out) :: covariance(:,:) !! Approximate covariance matrix from the negative inverse numerical.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported gradients, score failures, or a singular.
      real(dp), intent(in), optional :: trials(:) !! Optional binomial trial counts, size `n`; omitted values imply Bernoulli.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant design matrix, shape `(n,q)`, matching stored.
      real(dp), intent(in), optional :: difference_step !! Positive score-difference step; default is `1e-5`.
      real(dp), allocatable :: hessian(:,:), observed(:,:), plus(:), minus(:), gplus(:), gminus(:), rhs(:), solution(:)
      real(dp) :: h, ridge
      integer :: j, m, solve_info, attempt

      if (.not. flexmix_exist_gradient(result)) then
         allocate(estimate(0), standard_error(0), covariance(0,0))
         info = 30
         return
      end if
      h = 1.0e-5_dp
      if (present(difference_step)) h = difference_step
      if (h <= 0.0_dp) then
         allocate(estimate(0), standard_error(0), covariance(0,0))
         info = 31
         return
      end if
      call flexmix_parameter_vector(result, estimate, info)
      if (info /= 0) then
         allocate(standard_error(0), covariance(0,0))
         return
      end if
      m = size(estimate)
      allocate(hessian(m,m), observed(m,m), covariance(m,m), standard_error(m), plus(m), minus(m), rhs(m), solution(m))
      hessian = 0.0_dp
      do j = 1, m
         plus = estimate
         minus = estimate
         plus(j) = plus(j) + h
         minus(j) = minus(j) - h
         call flexmix_gradient_regression(result, x, y, plus, gplus, info, trials, offset, case_weights, concomitant_x)
         if (info /= 0) return
         call flexmix_gradient_regression(result, x, y, minus, gminus, info, trials, offset, case_weights, concomitant_x)
         if (info /= 0) return
         hessian(:,j) = (gplus - gminus) / (2.0_dp * h)
      end do
      hessian = 0.5_dp * (hessian + transpose(hessian))
      observed = -hessian
      covariance = 0.0_dp
      ridge = 0.0_dp
      do attempt = 1, 8
         solve_info = 0
         do j = 1, m
            rhs = 0.0_dp
            rhs(j) = 1.0_dp
            call solve_information_system(observed, ridge, rhs, solution, solve_info)
            if (solve_info /= 0) exit
            covariance(:,j) = solution
         end do
         if (solve_info == 0) exit
         if (attempt == 1) then
            ridge = 1.0e-10_dp
         else
            ridge = 10.0_dp * ridge
         end if
      end do
      if (solve_info /= 0) then
         standard_error = huge(1.0_dp)
         info = 32
         return
      end if
      covariance = 0.5_dp * (covariance + transpose(covariance))
      do j = 1, m
         standard_error(j) = sqrt(max(0.0_dp,covariance(j,j)))
      end do
      info = 0
   end subroutine flexmix_refit_optim_regression

   subroutine solve_information_system(observed, ridge, rhs, solution, info)
      use flexmix_numeric, only : solve_linear_system
      real(dp), intent(in) :: observed(:,:) !! Symmetric observed-information matrix before optional diagonal regularization.
      real(dp), intent(in) :: ridge !! Nonnegative diagonal regularization added to stabilize a nearly singular information.
      real(dp), intent(in) :: rhs(:) !! Right-hand side vector for one inverse-information column.
      real(dp), intent(out) :: solution(:) !! Solution of `(observed + ridge*I) * solution = rhs`.
      integer, intent(out) :: info !! Zero when the linear solve succeeds; nonzero when the regularized matrix is singular.
      real(dp), allocatable :: matrix(:,:)
      integer :: i

      matrix = observed
      do i = 1, size(matrix,1)
         matrix(i,i) = matrix(i,i) + ridge
      end do
      call solve_linear_system(matrix, rhs, solution, info)
   end subroutine solve_information_system

   subroutine fixed_parameter_design(result, design, info)
      type(flexmix_result), intent(in) :: result !! Fixed/shared-coefficient result carrying its coefficient-incidence matrix.
      logical, allocatable, intent(out) :: design(:,:) !! Incidence matrix including shared scale/shape columns when applicable.
      integer, intent(out) :: info !! Zero on success; nonzero when fixed-design metadata are incomplete.
      integer :: q, extra, g, ng

      if (.not. allocated(result%parameter_design)) then
         allocate(design(0,0))
         info = 1
         return
      end if
      q = size(result%parameter_design,2)
      extra = 0
      if (result%model_kind == flexmix_model_gaussian) then
         if (.not. allocated(result%variance_group)) then
            allocate(design(0,0))
            info = 2
            return
         end if
         extra = maxval(result%variance_group)
      else if (result%model_kind == flexmix_model_gamma_regression) then
         extra = 1
      end if
      allocate(design(result%k,q+extra))
      design = .false.
      design(:,1:q) = result%parameter_design
      if (result%model_kind == flexmix_model_gaussian) then
         ng = extra
         do g = 1, ng
            design(:,q+g) = result%variance_group == g
         end do
      else if (result%model_kind == flexmix_model_gamma_regression) then
         design(:,q+1) = .true.
      end if
      info = 0
   end subroutine fixed_parameter_design

   subroutine pack_fixed_parameters(result, parameters, info)
      type(flexmix_result), intent(in) :: result !! Fixed/shared result whose unique numerical parameters are packed.
      real(dp), allocatable, intent(out) :: parameters(:) !! Unique coefficient, scale/shape, and prior parameters in refit.
      integer, intent(out) :: info !! Zero on success; nonzero when fixed-design metadata or family arrays are incomplete.
      integer :: q, extra, prior_count, c, j, g, pos, ng

      q = size(result%parameter_design,2)
      extra = 0
      if (result%model_kind == flexmix_model_gaussian) then
         if (.not. allocated(result%variance_group) .or. .not. allocated(result%sigma)) then
            allocate(parameters(0))
            info = 1
            return
         end if
         extra = maxval(result%variance_group)
      else if (result%model_kind == flexmix_model_gamma_regression) then
         if (.not. allocated(result%shape)) then
            allocate(parameters(0))
            info = 1
            return
         end if
         extra = 1
      end if
      prior_count = prior_parameter_count(result)
      allocate(parameters(q+extra+prior_count))
      pos = 1
      do c = 1, q
         j = findloc(result%parameter_design(:,c),.true.,dim=1)
         if (j < 1) then
            info = 2
            return
         end if
         parameters(pos) = result%beta(c,j)
         pos = pos + 1
      end do
      if (result%model_kind == flexmix_model_gaussian) then
         ng = maxval(result%variance_group)
         do g = 1, ng
            j = findloc(result%variance_group,g,dim=1)
            parameters(pos) = log(max(result%sigma(j),tiny(1.0_dp)))
            pos = pos + 1
         end do
      else if (result%model_kind == flexmix_model_gamma_regression) then
         parameters(pos) = result%shape(1)
         pos = pos + 1
      end if
      call pack_prior_parameters(result, parameters, pos, info)
   end subroutine pack_fixed_parameters

   subroutine unpack_fixed_parameters(result, parameters, updated, info)
      type(flexmix_result), intent(in) :: result !! Fixed/shared template defining unique parameter incidence and scale sharing.
      real(dp), intent(in) :: parameters(:) !! Replacement vector in the order returned by `flexmix_parameter_vector`.
      type(flexmix_result), intent(out) :: updated !! Copy with shared parameters propagated to every affected component.
      integer, intent(out) :: info !! Zero on success; nonzero for length mismatch or incomplete fixed-design metadata.
      integer :: q, extra, expected, c, j, g, pos, ng

      q = size(result%parameter_design,2)
      extra = 0
      if (result%model_kind == flexmix_model_gaussian) then
         if (.not. allocated(result%variance_group)) then
            updated = result
            info = 1
            return
         end if
         extra = maxval(result%variance_group)
      else if (result%model_kind == flexmix_model_gamma_regression) then
         extra = 1
      end if
      expected = q + extra + prior_parameter_count(result)
      updated = result
      if (size(parameters) /= expected) then
         info = 2
         return
      end if
      pos = 1
      updated%beta = 0.0_dp
      do c = 1, q
         do j = 1, result%k
            if (result%parameter_design(j,c)) updated%beta(c,j) = parameters(pos)
         end do
         pos = pos + 1
      end do
      if (result%model_kind == flexmix_model_gaussian) then
         ng = maxval(result%variance_group)
         do g = 1, ng
            do j = 1, result%k
               if (result%variance_group(j) == g) updated%sigma(j) = exp(parameters(pos))
            end do
            pos = pos + 1
         end do
      else if (result%model_kind == flexmix_model_gamma_regression) then
         updated%shape = max(parameters(pos),tiny(1.0_dp))
         pos = pos + 1
      end if
      call unpack_prior_parameters(updated, parameters, pos, info)
   end subroutine unpack_fixed_parameters

   subroutine component_parameter_counts(result, counts, info)
      type(flexmix_result), intent(in) :: result !! Fitted result whose optimizable component parameter counts are required.
      integer, allocatable, intent(out) :: counts(:) !! Number of packed component parameters for each retained component.
      integer, intent(out) :: info !! Zero on success; nonzero when the family-specific result arrays are unavailable.
      integer :: j, d

      allocate(counts(max(0,result%k)))
      counts = 0
      info = 0
      do j = 1, result%k
         if (allocated(result%component_role)) then
            if (result%component_role(j) /= flexmix_role_regular) cycle
         end if
         select case (result%model_kind)
         case (flexmix_model_gaussian, flexmix_model_robust_gaussian)
            if (.not. allocated(result%beta) .or. .not. allocated(result%sigma)) then
               info = 1
               return
            end if
            counts(j) = size(result%beta,1) + 1
         case (flexmix_model_poisson, flexmix_model_binomial, flexmix_model_ziglm_poisson, &
               flexmix_model_ziglm_binomial, flexmix_model_robust_poisson)
            if (.not. allocated(result%beta)) then
               info = 1
               return
            end if
            counts(j) = size(result%beta,1)
         case (flexmix_model_gamma_regression)
            if (.not. allocated(result%beta) .or. .not. allocated(result%shape)) then
               info = 1
               return
            end if
            counts(j) = size(result%beta,1) + 1
         case (flexmix_model_multinomial)
            if (.not. allocated(result%multinomial_coef)) then
               info = 1
               return
            end if
            counts(j) = size(result%multinomial_coef,1) * size(result%multinomial_coef,2)
         case (flexmix_model_mvnorm)
            if (.not. allocated(result%center) .or. .not. allocated(result%covariance)) then
               info = 1
               return
            end if
            d = result%d
            if (result%diagonal_covariance) then
               counts(j) = 2*d
            else
               counts(j) = d + d*(d+1)/2
            end if
         case (flexmix_model_mvbinary, flexmix_model_mvpois)
            counts(j) = result%d
         case (flexmix_model_mvcombi)
            if (.not. allocated(result%binary_mask)) then
               info = 1
               return
            end if
            counts(j) = result%d + count(.not. result%binary_mask)
         case (flexmix_model_lognormal, flexmix_model_inverse_gaussian, flexmix_model_gamma, flexmix_model_weibull)
            counts(j) = 2
         case (flexmix_model_exponential)
            counts(j) = 1
         case default
            info = 1
            return
         end select
      end do
   end subroutine component_parameter_counts

   pure integer function prior_parameter_count(result) result(count_value)
      type(flexmix_result), intent(in) :: result !! Fitted result whose independent prior/concomitant coefficients are counted.
      if (result%k <= 1) then
         count_value = 0
      else if (allocated(result%concomitant_coef)) then
         count_value = size(result%concomitant_coef,1) * (result%k - 1)
      else
         count_value = result%k - 1
      end if
   end function prior_parameter_count

   subroutine pack_component_parameters(result, component, parameters, pos, info)
      type(flexmix_result), intent(in) :: result !! Fitted result supplying one component's numerical parameters.
      integer, intent(in) :: component !! One-based retained component index to pack.
      real(dp), intent(inout) :: parameters(:) !! Destination parameter vector being filled in family-specific order.
      integer, intent(inout) :: pos !! One-based insertion position, advanced past the packed component block.
      integer, intent(out) :: info !! Zero on success; nonzero for an unsupported family or malformed stored arrays.
      integer :: a, b, d, p

      info = 0
      if (allocated(result%component_role)) then
         if (result%component_role(component) /= flexmix_role_regular) return
      end if
      select case (result%model_kind)
      case (flexmix_model_gaussian, flexmix_model_robust_gaussian)
         p = size(result%beta,1)
         parameters(pos:pos+p-1) = result%beta(:,component)
         parameters(pos+p) = log(max(result%sigma(component),tiny(1.0_dp)))
         pos = pos + p + 1
      case (flexmix_model_poisson, flexmix_model_binomial, flexmix_model_ziglm_poisson, &
            flexmix_model_ziglm_binomial, flexmix_model_robust_poisson)
         p = size(result%beta,1)
         parameters(pos:pos+p-1) = result%beta(:,component)
         pos = pos + p
      case (flexmix_model_gamma_regression)
         p = size(result%beta,1)
         parameters(pos:pos+p-1) = result%beta(:,component)
         parameters(pos+p) = result%shape(component)
         pos = pos + p + 1
      case (flexmix_model_multinomial)
         p = size(result%multinomial_coef,1) * size(result%multinomial_coef,2)
         parameters(pos:pos+p-1) = reshape(result%multinomial_coef(:,:,component), [p])
         pos = pos + p
      case (flexmix_model_mvnorm)
         d = result%d
         parameters(pos:pos+d-1) = result%center(:,component)
         pos = pos + d
         if (result%diagonal_covariance) then
            do a = 1, d
               parameters(pos) = result%covariance(a,a,component)
               pos = pos + 1
            end do
         else
            do b = 1, d
               do a = b, d
                  parameters(pos) = result%covariance(a,b,component)
                  pos = pos + 1
               end do
            end do
         end if
      case (flexmix_model_mvbinary)
         parameters(pos:pos+result%d-1) = result%probability(:,component)
         pos = pos + result%d
      case (flexmix_model_mvpois)
         parameters(pos:pos+result%d-1) = result%lambda(:,component)
         pos = pos + result%d
      case (flexmix_model_mvcombi)
         parameters(pos:pos+result%d-1) = result%center(:,component)
         pos = pos + result%d
         do a = 1, result%d
            if (.not. result%binary_mask(a)) then
               parameters(pos) = result%covariance(a,a,component)
               pos = pos + 1
            end if
         end do
      case (flexmix_model_lognormal)
         parameters(pos:pos+1) = [result%center(1,component), result%sigma(component)]
         pos = pos + 2
      case (flexmix_model_exponential)
         parameters(pos) = result%rate(component)
         pos = pos + 1
      case (flexmix_model_inverse_gaussian)
         parameters(pos:pos+1) = [result%center(1,component), result%lambda(1,component)]
         pos = pos + 2
      case (flexmix_model_gamma)
         parameters(pos:pos+1) = [result%shape(component), result%rate(component)]
         pos = pos + 2
      case (flexmix_model_weibull)
         parameters(pos:pos+1) = [result%shape(component), result%scale(component)]
         pos = pos + 2
      case default
         info = 1
      end select
   end subroutine pack_component_parameters

   subroutine unpack_component_parameters(result, component, parameters, pos, info)
      type(flexmix_result), intent(inout) :: result !! Result whose one component is replaced from the packed parameter vector.
      integer, intent(in) :: component !! One-based retained component index to replace.
      real(dp), intent(in) :: parameters(:) !! Source parameter vector in `flexmix_parameter_vector` order.
      integer, intent(inout) :: pos !! One-based source position, advanced past the unpacked component block.
      integer, intent(out) :: info !! Zero on success; nonzero for an unsupported family.
      integer :: a, b, d, p
      real(dp) :: value

      info = 0
      if (allocated(result%component_role)) then
         if (result%component_role(component) /= flexmix_role_regular) return
      end if
      select case (result%model_kind)
      case (flexmix_model_gaussian, flexmix_model_robust_gaussian)
         p = size(result%beta,1)
         result%beta(:,component) = parameters(pos:pos+p-1)
         result%sigma(component) = exp(parameters(pos+p))
         pos = pos + p + 1
      case (flexmix_model_poisson, flexmix_model_binomial, flexmix_model_ziglm_poisson, &
            flexmix_model_ziglm_binomial, flexmix_model_robust_poisson)
         p = size(result%beta,1)
         result%beta(:,component) = parameters(pos:pos+p-1)
         pos = pos + p
      case (flexmix_model_gamma_regression)
         p = size(result%beta,1)
         result%beta(:,component) = parameters(pos:pos+p-1)
         result%shape(component) = max(parameters(pos+p),tiny(1.0_dp))
         pos = pos + p + 1
      case (flexmix_model_multinomial)
         p = size(result%multinomial_coef,1) * size(result%multinomial_coef,2)
         result%multinomial_coef(:,:,component) = reshape(parameters(pos:pos+p-1), shape(result%multinomial_coef(:,:,component)))
         pos = pos + p
      case (flexmix_model_mvnorm)
         d = result%d
         result%center(:,component) = parameters(pos:pos+d-1)
         pos = pos + d
         result%covariance(:,:,component) = 0.0_dp
         if (result%diagonal_covariance) then
            do a = 1, d
               result%covariance(a,a,component) = max(parameters(pos),tiny(1.0_dp))
               pos = pos + 1
            end do
         else
            do b = 1, d
               do a = b, d
                  value = parameters(pos)
                  result%covariance(a,b,component) = value
                  result%covariance(b,a,component) = value
                  pos = pos + 1
               end do
            end do
         end if
      case (flexmix_model_mvbinary)
         result%probability(:,component) = min(1.0_dp-tiny(1.0_dp),max(tiny(1.0_dp),parameters(pos:pos+result%d-1)))
         pos = pos + result%d
      case (flexmix_model_mvpois)
         result%lambda(:,component) = max(tiny(1.0_dp),parameters(pos:pos+result%d-1))
         pos = pos + result%d
      case (flexmix_model_mvcombi)
         result%center(:,component) = parameters(pos:pos+result%d-1)
         pos = pos + result%d
         do a = 1, result%d
            if (result%binary_mask(a)) then
               result%center(a,component) = min(1.0_dp-tiny(1.0_dp),max(tiny(1.0_dp),result%center(a,component)))
            else
               result%covariance(a,a,component) = max(parameters(pos),tiny(1.0_dp))
               pos = pos + 1
            end if
         end do
      case (flexmix_model_lognormal)
         result%center(1,component) = parameters(pos)
         result%sigma(component) = max(parameters(pos+1),tiny(1.0_dp))
         pos = pos + 2
      case (flexmix_model_exponential)
         result%rate(component) = max(parameters(pos),tiny(1.0_dp))
         pos = pos + 1
      case (flexmix_model_inverse_gaussian)
         result%center(1,component) = max(parameters(pos),tiny(1.0_dp))
         result%lambda(1,component) = max(parameters(pos+1),tiny(1.0_dp))
         pos = pos + 2
      case (flexmix_model_gamma)
         result%shape(component) = max(parameters(pos),tiny(1.0_dp))
         result%rate(component) = max(parameters(pos+1),tiny(1.0_dp))
         pos = pos + 2
      case (flexmix_model_weibull)
         result%shape(component) = max(parameters(pos),tiny(1.0_dp))
         result%scale(component) = max(parameters(pos+1),tiny(1.0_dp))
         pos = pos + 2
      case default
         info = 1
      end select
   end subroutine unpack_component_parameters

   subroutine pack_prior_parameters(result, parameters, pos, info)
      type(flexmix_result), intent(in) :: result !! Fitted result supplying constant or multinomial-concomitant prior.
      real(dp), intent(inout) :: parameters(:) !! Destination full parameter vector being filled after component blocks.
      integer, intent(inout) :: pos !! One-based insertion position, advanced past prior parameters.
      integer, intent(out) :: info !! Zero on success; nonzero when stored priors cannot be transformed safely.
      integer :: j, q

      info = 0
      if (result%k <= 1) return
      if (allocated(result%concomitant_coef)) then
         q = size(result%concomitant_coef,1)
         do j = 2, result%k
            parameters(pos:pos+q-1) = result%concomitant_coef(:,j) - result%concomitant_coef(:,1)
            pos = pos + q
         end do
      else
         if (.not. allocated(result%prior) .or. any(result%prior <= 0.0_dp)) then
            info = 1
            return
         end if
         do j = 2, result%k
            parameters(pos) = log(result%prior(j)) - log(result%prior(1))
            pos = pos + 1
         end do
      end if
   end subroutine pack_prior_parameters

   subroutine unpack_prior_parameters(result, parameters, pos, info)
      type(flexmix_result), intent(inout) :: result !! Result whose constant or multinomial-concomitant priors are replaced.
      real(dp), intent(in) :: parameters(:) !! Source full parameter vector containing prior parameters after component blocks.
      integer, intent(inout) :: pos !! One-based source position, advanced past prior parameters.
      integer, intent(out) :: info !! Zero on success; nonzero when prior storage is malformed.
      real(dp), allocatable :: logits(:)
      real(dp) :: maximum
      integer :: j, q

      info = 0
      if (result%k <= 1) return
      if (allocated(result%concomitant_coef)) then
         q = size(result%concomitant_coef,1)
         result%concomitant_coef(:,1) = 0.0_dp
         do j = 2, result%k
            result%concomitant_coef(:,j) = parameters(pos:pos+q-1)
            pos = pos + q
         end do
      else
         allocate(logits(result%k))
         logits(1) = 0.0_dp
         do j = 2, result%k
            logits(j) = parameters(pos)
            pos = pos + 1
         end do
         maximum = maxval(logits)
         logits = exp(logits - maximum)
         if (.not. allocated(result%prior)) allocate(result%prior(result%k))
         result%prior = logits / sum(logits)
      end if
   end subroutine unpack_prior_parameters

   subroutine prior_matrix(result, n, prior_rows, info, concomitant_x)
      type(flexmix_result), intent(in) :: result !! Fitted result supplying constant priors or stored multinomial-concomitant.
      integer, intent(in) :: n !! Number of observations for which prior probabilities are requested.
      real(dp), intent(out) :: prior_rows(:,:) !! Prior probability matrix, shape `(n,k)`.
      integer, intent(out) :: info !! Zero on success; nonzero when the supplied concomitant design is missing or dimensionally.
      real(dp), allocatable :: eta(:,:)
      real(dp) :: maximum, denominator
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant design `(n,q)`, required for softmax priors.
      integer :: i

      if (size(prior_rows,1) /= n .or. size(prior_rows,2) /= result%k) then
         info = 1
         return
      end if
      if (allocated(result%concomitant_coef)) then
         if (.not. present(concomitant_x)) then
            info = 2
            return
         end if
         if (size(concomitant_x,1) /= n .or. size(concomitant_x,2) /= size(result%concomitant_coef,1)) then
            info = 3
            return
         end if
         allocate(eta(n,result%k))
         eta = matmul(concomitant_x,result%concomitant_coef)
         do i = 1, n
            maximum = maxval(eta(i,:))
            prior_rows(i,:) = exp(eta(i,:) - maximum)
            denominator = sum(prior_rows(i,:))
            prior_rows(i,:) = prior_rows(i,:) / denominator
         end do
      else
         if (.not. allocated(result%prior) .or. size(result%prior) /= result%k) then
            info = 4
            return
         end if
         do i = 1, n
            prior_rows(i,:) = result%prior
         end do
      end if
      info = 0
   end subroutine prior_matrix

end module flexmix_refit
