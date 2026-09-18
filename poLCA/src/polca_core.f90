module polca_core
   use polca_kinds, only : dp
   use polca_linalg, only : solve_linear, symmetric_pinv
   use polca_types, only : polca_model
   implicit none
   private
   public :: polca_fit
   public :: polca_item_likelihood
   public :: polca_postclass
   public :: polca_probhat
   public :: polca_update_prior
   public :: polca_beta_derivatives
   public :: polca_standard_errors

contains

   subroutine polca_item_likelihood(probs, n_choices, y, likelihood)
      real(dp), intent(in) :: probs(:, :, :) !! Class-by-category-by-item response probabilities.
      integer, intent(in) :: n_choices(:) !! Number of valid response categories for each manifest item.
      integer, intent(in) :: y(:, :) !! Observation-by-item integer responses; zero denotes a missing response.
      real(dp), intent(out) :: likelihood(:, :) !! Observation-by-class conditional response likelihoods.
      real(dp) :: logp
      integer :: i, j, r, n, nclass

      n = size(y, 1)
      nclass = size(probs, 1)
      likelihood = 0.0_dp
      do i = 1, n
         do r = 1, nclass
            logp = 0.0_dp
            do j = 1, size(y, 2)
               if (y(i, j) > 0) then
                  if (y(i, j) <= n_choices(j)) then
                     logp = logp + log(max(probs(r, y(i, j), j), tiny(1.0_dp)))
                  else
                     logp = log(tiny(1.0_dp))
                     exit
                  end if
               end if
            end do
            likelihood(i, r) = exp(logp)
         end do
      end do
   end subroutine polca_item_likelihood

   subroutine polca_postclass(prior, probs, n_choices, y, posterior, loglik)
      real(dp), intent(in) :: prior(:, :) !! Observation-by-class prior membership probabilities.
      real(dp), intent(in) :: probs(:, :, :) !! Class-by-category-by-item response probabilities.
      integer, intent(in) :: n_choices(:) !! Number of response categories for each manifest item.
      integer, intent(in) :: y(:, :) !! Observation-by-item responses with zero used for missing values.
      real(dp), intent(out) :: posterior(:, :) !! Observation-by-class posterior membership probabilities.
      real(dp), intent(out), optional :: loglik !! Observed-data log likelihood summed over observations.
      real(dp), allocatable :: z(:)
      real(dp) :: amax, denom, ll
      integer :: i, j, r, nclass

      nclass = size(prior, 2)
      allocate(z(nclass))
      ll = 0.0_dp
      do i = 1, size(y, 1)
         do r = 1, nclass
            z(r) = log(max(prior(i, r), tiny(1.0_dp)))
            do j = 1, size(y, 2)
               if (y(i, j) > 0) then
                  if (y(i, j) <= n_choices(j)) then
                     z(r) = z(r) + log(max(probs(r, y(i, j), j), tiny(1.0_dp)))
                  else
                     z(r) = log(tiny(1.0_dp))
                     exit
                  end if
               end if
            end do
         end do
         amax = maxval(z)
         z = exp(z - amax)
         denom = sum(z)
         if (denom <= tiny(1.0_dp)) then
            posterior(i, :) = 1.0_dp / real(nclass, dp)
            ll = ll + log(tiny(1.0_dp))
         else
            posterior(i, :) = z / denom
            ll = ll + amax + log(denom)
         end if
      end do
      if (present(loglik)) loglik = ll
   end subroutine polca_postclass

   subroutine polca_probhat(posterior, y, n_choices, probs)
      real(dp), intent(in) :: posterior(:, :) !! Observation-by-class posterior probabilities.
      integer, intent(in) :: y(:, :) !! Observation-by-item responses; zero denotes missing data.
      integer, intent(in) :: n_choices(:) !! Number of response categories for each item.
      real(dp), intent(out) :: probs(:, :, :) !! Updated class-by-category-by-item response probabilities.
      real(dp) :: denom
      integer :: j, k, r

      probs = 0.0_dp
      do r = 1, size(posterior, 2)
         do j = 1, size(y, 2)
            denom = sum(posterior(:, r), mask=y(:, j) > 0)
            if (denom <= tiny(1.0_dp)) then
               probs(r, 1:n_choices(j), j) = 1.0_dp / real(n_choices(j), dp)
            else
               do k = 1, n_choices(j)
                  probs(r, k, j) = sum(posterior(:, r), mask=y(:, j) == k) / denom
               end do
            end if
         end do
      end do
   end subroutine polca_probhat

   subroutine polca_update_prior(coeff, x, prior)
      real(dp), intent(in) :: coeff(:, :) !! Predictor coefficients for classes 2..R relative to reference class 1.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix including an intercept column when desired.
      real(dp), intent(out) :: prior(:, :) !! Observation-by-class multinomial-logit prior probabilities.
      real(dp), allocatable :: eta(:)
      real(dp) :: amax, denom
      integer :: i, r, nclass

      nclass = size(coeff, 2) + 1
      allocate(eta(nclass))
      do i = 1, size(x, 1)
         eta(1) = 0.0_dp
         do r = 2, nclass
            eta(r) = dot_product(x(i, :), coeff(:, r - 1))
         end do
         amax = maxval(eta)
         eta = exp(eta - amax)
         denom = sum(eta)
         prior(i, :) = eta / denom
      end do
   end subroutine polca_update_prior

   subroutine polca_beta_derivatives(posterior, prior, x, grad, neg_hess)
      real(dp), intent(in) :: posterior(:, :) !! Observation-by-class posterior probabilities.
      real(dp), intent(in) :: prior(:, :) !! Observation-by-class prior probabilities at current coefficients.
      real(dp), intent(in) :: x(:, :) !! Predictor design matrix including any intercept column.
      real(dp), intent(out) :: grad(:) !! Score vector for non-reference multinomial-logit coefficients.
      real(dp), intent(out) :: neg_hess(:, :) !! Negative observed Hessian matching poLCA's native update.
      real(dp) :: w
      integer :: a, b, i, j, k, r, s, nclass, nx

      nclass = size(prior, 2)
      nx = size(x, 2)
      grad = 0.0_dp
      neg_hess = 0.0_dp

      do i = 1, size(x, 1)
         do r = 2, nclass
            do j = 1, nx
               a = (r - 2) * nx + j
               grad(a) = grad(a) + x(i, j) * (posterior(i, r) - prior(i, r))
            end do
         end do
         do r = 2, nclass
            do s = 2, nclass
               if (r == s) then
                  w = prior(i, r) * (1.0_dp - prior(i, r)) &
                      - posterior(i, r) * (1.0_dp - posterior(i, r))
               else
                  w = posterior(i, r) * posterior(i, s) - prior(i, r) * prior(i, s)
               end if
               do j = 1, nx
                  a = (r - 2) * nx + j
                  do k = 1, nx
                     b = (s - 2) * nx + k
                     neg_hess(a, b) = neg_hess(a, b) + x(i, j) * x(i, k) * w
                  end do
               end do
            end do
         end do
      end do
   end subroutine polca_beta_derivatives

   subroutine polca_fit(y, n_choices, nclass, model, x, probs_start, maxiter, tol, nrep, seed, calc_se)
      integer, intent(in) :: y(:, :) !! Observation-by-item categorical responses; use zero for missing values.
      integer, intent(in) :: n_choices(:) !! Number of categories for each manifest item.
      integer, intent(in) :: nclass !! Number of latent classes; must be at least one.
      type(polca_model), intent(out) :: model !! Fitted latent-class or latent-class-regression model.
      real(dp), intent(in), optional :: x(:, :) !! Optional predictor matrix; include an intercept as column one.
      real(dp), intent(in), optional :: probs_start(:, :, :) !! Optional class-by-category-by-item starting probabilities.
      integer, intent(in), optional :: maxiter !! Maximum EM iterations per restart; default is 1000.
      real(dp), intent(in), optional :: tol !! Absolute log-likelihood improvement tolerance; default is 1e-10.
      integer, intent(in), optional :: nrep !! Number of independent starts; default is one.
      integer, intent(in), optional :: seed !! Deterministic seed used for generated starts.
      logical, intent(in), optional :: calc_se !! Whether to compute outer-product standard errors; default is true.
      type(polca_model) :: candidate
      real(dp), allocatable :: xx(:, :)
      real(dp) :: use_tol
      integer :: irep, nr, mi, use_seed
      logical :: do_se

      mi = 1000
      if (present(maxiter)) mi = maxiter
      use_tol = 1.0e-10_dp
      if (present(tol)) use_tol = tol
      nr = 1
      if (present(nrep)) nr = max(1, nrep)
      use_seed = 12345
      if (present(seed)) use_seed = seed
      do_se = .true.
      if (present(calc_se)) do_se = calc_se

      call validate_inputs(y, n_choices, nclass)
      if (present(x)) then
         if (size(x, 1) /= size(y, 1)) error stop "polca_fit: x and y row counts differ"
         allocate(xx(size(x, 1), size(x, 2)))
         xx = x
      else
         allocate(xx(size(y, 1), 1))
         xx = 1.0_dp
      end if

      call set_random_seed(use_seed)
      model%loglik = -huge(1.0_dp)
      do irep = 1, nr
         if (present(probs_start) .and. irep == 1) then
            call fit_once(y, n_choices, nclass, xx, candidate, mi, use_tol, probs_start)
         else
            call fit_once(y, n_choices, nclass, xx, candidate, mi, use_tol)
         end if
         if (candidate%loglik > model%loglik) model = candidate
      end do

      if (do_se) call polca_standard_errors(y, xx, model)
      call compute_fit_statistics(y, model)
   end subroutine polca_fit

   subroutine validate_inputs(y, n_choices, nclass)
      integer, intent(in) :: y(:, :) !! Observation-by-item categorical responses to validate.
      integer, intent(in) :: n_choices(:) !! Number of categories declared for each item.
      integer, intent(in) :: nclass !! Requested number of latent classes.
      integer :: j

      if (nclass < 1) error stop "polca_fit: nclass must be positive"
      if (size(y, 2) /= size(n_choices)) error stop "polca_fit: n_choices has wrong length"
      if (any(n_choices < 1)) error stop "polca_fit: each item must have at least one category"
      if (any(y < 0)) error stop "polca_fit: responses must be nonnegative integers"
      do j = 1, size(y, 2)
         if (any(y(:, j) > n_choices(j))) error stop "polca_fit: response exceeds declared category count"
      end do
   end subroutine validate_inputs

   subroutine fit_once(y, n_choices, nclass, x, model, maxiter, tol, probs_start)
      integer, intent(in) :: y(:, :) !! Observation-by-item categorical responses with zero for missing data.
      integer, intent(in) :: n_choices(:) !! Number of categories for each manifest item.
      integer, intent(in) :: nclass !! Number of latent classes.
      real(dp), intent(in) :: x(:, :) !! Predictor design matrix including the intercept if desired.
      type(polca_model), intent(out) :: model !! Model fitted from one set of initial response probabilities.
      integer, intent(in) :: maxiter !! Maximum number of EM updates.
      real(dp), intent(in) :: tol !! Absolute log-likelihood improvement tolerance.
      real(dp), intent(in), optional :: probs_start(:, :, :) !! Optional response-probability starting array.
      real(dp), allocatable :: probs(:, :, :), prior(:, :), posterior(:, :), coeff(:, :)
      real(dp), allocatable :: grad(:), neg_hess(:, :), delta(:)
      real(dp) :: ll, ll_old, dll
      integer :: i, j, k, r, iter, maxk, nx, q
      logical :: ok, use_covariates

      maxk = maxval(n_choices)
      nx = size(x, 2)
      use_covariates = nx > 1 .and. nclass > 1
      allocate(probs(nclass, maxk, size(y, 2)))
      allocate(prior(size(y, 1), nclass), posterior(size(y, 1), nclass))
      allocate(coeff(nx, max(0, nclass - 1)))
      coeff = 0.0_dp

      if (present(probs_start)) then
         if (size(probs_start, 1) /= nclass .or. size(probs_start, 2) < maxk &
             .or. size(probs_start, 3) /= size(y, 2)) then
            error stop "polca_fit: probs_start has incompatible shape"
         end if
         probs = probs_start(:, 1:maxk, :)
         call normalize_probs(probs, n_choices)
      else
         call random_start(probs, n_choices)
      end if

      if (nclass == 1) then
         prior = 1.0_dp
      else
         call polca_update_prior(coeff, x, prior)
      end if
      ll_old = -huge(1.0_dp)
      ll = ll_old

      if (use_covariates) then
         q = nx * (nclass - 1)
         allocate(grad(q), neg_hess(q, q), delta(q))
      end if

      do iter = 1, maxiter
         call polca_postclass(prior, probs, n_choices, y, posterior, ll)
         call polca_probhat(posterior, y, n_choices, probs)
         if (nclass > 1) then
            if (use_covariates) then
               call polca_beta_derivatives(posterior, prior, x, grad, neg_hess)
               call solve_linear(neg_hess, grad, delta, ok)
               if (ok) then
                  do r = 2, nclass
                     coeff(:, r - 1) = coeff(:, r - 1) + delta((r - 2) * nx + 1:(r - 1) * nx)
                  end do
               end if
               call polca_update_prior(coeff, x, prior)
            else
               do r = 1, nclass
                  prior(:, r) = sum(posterior(:, r)) / real(size(y, 1), dp)
               end do
            end if
         end if
         call polca_postclass(prior, probs, n_choices, y, posterior, ll)
         if (iter > 1) then
            dll = ll - ll_old
            if (abs(dll) <= tol) exit
         end if
         ll_old = ll
      end do

      model%n = size(y, 1)
      model%n_items = size(y, 2)
      model%n_classes = nclass
      model%n_predictors = nx
      model%numiter = min(iter, maxiter)
      model%converged = iter <= maxiter
      model%has_covariates = use_covariates
      model%loglik = ll
      allocate(model%n_choices(size(n_choices)))
      model%n_choices = n_choices
      allocate(model%probs(nclass, maxk, size(y, 2)))
      model%probs = probs
      allocate(model%class_share(nclass))
      model%class_share = sum(posterior, dim=1) / real(size(y, 1), dp)
      allocate(model%posterior(size(y, 1), nclass))
      model%posterior = posterior
      allocate(model%predclass(size(y, 1)))
      do i = 1, size(y, 1)
         model%predclass(i) = maxloc(posterior(i, :), dim=1)
      end do
      if (nclass > 1 .and. use_covariates) then
         allocate(model%coeff(nx, nclass - 1))
         model%coeff = coeff
      else
         allocate(model%coeff(0, 0))
      end if

      model%npar = nclass * sum(n_choices - 1) + (nclass - 1)
      if (use_covariates) model%npar = model%npar + nx * (nclass - 1) - (nclass - 1)
      model%aic = -2.0_dp * ll + 2.0_dp * real(model%npar, dp)
      model%bic = -2.0_dp * ll + log(real(model%n, dp)) * real(model%npar, dp)

      do j = 1, size(y, 2)
         do r = 1, nclass
            do k = n_choices(j) + 1, maxk
               model%probs(r, k, j) = 0.0_dp
            end do
         end do
      end do
   end subroutine fit_once

   subroutine normalize_probs(probs, n_choices)
      real(dp), intent(inout) :: probs(:, :, :) !! Class-by-category-by-item probabilities normalized in place.
      integer, intent(in) :: n_choices(:) !! Number of active categories for each item.
      real(dp) :: s
      integer :: j, r

      do j = 1, size(n_choices)
         do r = 1, size(probs, 1)
            probs(r, 1:n_choices(j), j) = max(probs(r, 1:n_choices(j), j), tiny(1.0_dp))
            s = sum(probs(r, 1:n_choices(j), j))
            probs(r, 1:n_choices(j), j) = probs(r, 1:n_choices(j), j) / s
            if (n_choices(j) < size(probs, 2)) probs(r, n_choices(j) + 1:, j) = 0.0_dp
         end do
      end do
   end subroutine normalize_probs

   subroutine random_start(probs, n_choices)
      real(dp), intent(out) :: probs(:, :, :) !! Generated class-by-category-by-item starting probabilities.
      integer, intent(in) :: n_choices(:) !! Number of active categories for each item.
      integer :: j, r

      probs = 0.0_dp
      do j = 1, size(n_choices)
         do r = 1, size(probs, 1)
            call random_number(probs(r, 1:n_choices(j), j))
         end do
      end do
      call normalize_probs(probs, n_choices)
   end subroutine random_start

   subroutine set_random_seed(seed)
      integer, intent(in) :: seed !! Scalar seed expanded deterministically to the compiler's random-seed vector.
      integer, allocatable :: put(:)
      integer :: i, n

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine set_random_seed

   subroutine polca_standard_errors(y, x, model)
      integer, intent(in) :: y(:, :) !! Fitted observation-by-item responses with zero for missing values.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix used by the fitted model.
      type(polca_model), intent(inout) :: model !! Fitted model receiving probability, share, and coefficient SEs.
      real(dp), allocatable :: prior(:, :), score(:, :), info(:, :), vce(:, :)
      real(dp), allocatable :: jac(:, :), vprob(:, :), jmix(:, :), vbeta(:, :), vmix(:, :)
      real(dp) :: indicator
      integer :: a, b, cpos, i, j, k, l, q, qbeta, qprob, r, rpos, nx, maxk

      nx = size(x, 2)
      maxk = size(model%probs, 2)
      qprob = model%n_classes * sum(model%n_choices - 1)
      qbeta = 0
      if (model%n_classes > 1) qbeta = nx * (model%n_classes - 1)
      q = qprob + qbeta
      if (q == 0) return

      allocate(prior(model%n, model%n_classes))
      if (model%has_covariates) then
         call polca_update_prior(model%coeff, x, prior)
      else
         do r = 1, model%n_classes
            prior(:, r) = model%class_share(r)
         end do
      end if

      allocate(score(model%n, q))
      score = 0.0_dp
      cpos = 0
      do r = 1, model%n_classes
         do j = 1, model%n_items
            do k = 2, model%n_choices(j)
               cpos = cpos + 1
               do i = 1, model%n
                  if (y(i, j) > 0) then
                     indicator = merge(1.0_dp, 0.0_dp, y(i, j) == k)
                     score(i, cpos) = model%posterior(i, r) * (indicator - model%probs(r, k, j))
                  end if
               end do
            end do
         end do
      end do
      if (model%n_classes > 1) then
         do r = 2, model%n_classes
            do l = 1, nx
               cpos = cpos + 1
               score(:, cpos) = x(:, l) * (model%posterior(:, r) - prior(:, r))
            end do
         end do
      end if

      allocate(info(q, q), vce(q, q))
      info = matmul(transpose(score), score)
      call symmetric_pinv(info, vce)

      allocate(jac(model%n_classes * sum(model%n_choices), qprob))
      jac = 0.0_dp
      rpos = 1
      cpos = 1
      do r = 1, model%n_classes
         do j = 1, model%n_items
            do k = 1, model%n_choices(j)
               do l = 2, model%n_choices(j)
                  if (k == l) then
                     indicator = model%probs(r, k, j) * (1.0_dp - model%probs(r, k, j))
                  else
                     indicator = -model%probs(r, k, j) * model%probs(r, l, j)
                  end if
                  jac(rpos + k - 1, cpos + l - 2) = indicator
               end do
            end do
            rpos = rpos + model%n_choices(j)
            cpos = cpos + model%n_choices(j) - 1
         end do
      end do
      allocate(vprob(size(jac, 1), size(jac, 1)))
      vprob = matmul(jac, matmul(vce(1:qprob, 1:qprob), transpose(jac)))
      allocate(model%probs_se(model%n_classes, maxk, model%n_items))
      model%probs_se = 0.0_dp
      rpos = 1
      do r = 1, model%n_classes
         do j = 1, model%n_items
            do k = 1, model%n_choices(j)
               model%probs_se(r, k, j) = sqrt(max(0.0_dp, vprob(rpos, rpos)))
               rpos = rpos + 1
            end do
         end do
      end do

      allocate(model%class_share_se(model%n_classes))
      model%class_share_se = 0.0_dp
      if (qbeta > 0) then
         allocate(vbeta(qbeta, qbeta))
         vbeta = vce(qprob + 1:q, qprob + 1:q)
         if (model%has_covariates) then
            allocate(model%coeff_se(nx, model%n_classes - 1))
            allocate(model%coeff_v(qbeta, qbeta))
            model%coeff_v = vbeta
            do r = 1, model%n_classes - 1
               do l = 1, nx
                  a = (r - 1) * nx + l
                  model%coeff_se(l, r) = sqrt(max(0.0_dp, vbeta(a, a)))
               end do
            end do
         else
            allocate(model%coeff_se(0, 0), model%coeff_v(0, 0))
         end if

         allocate(jmix(model%n_classes, qbeta))
         jmix = 0.0_dp
         do r = 2, model%n_classes
            do l = 1, nx
               b = (r - 2) * nx + l
               do i = 1, model%n
                  do k = 1, model%n_classes
                     if (k == r) then
                        indicator = prior(i, k) * (1.0_dp - prior(i, r))
                     else
                        indicator = -prior(i, k) * prior(i, r)
                     end if
                     jmix(k, b) = jmix(k, b) + indicator * x(i, l) / real(model%n, dp)
                  end do
               end do
            end do
         end do
         allocate(vmix(model%n_classes, model%n_classes))
         vmix = matmul(jmix, matmul(vbeta, transpose(jmix)))
         do r = 1, model%n_classes
            model%class_share_se(r) = sqrt(max(0.0_dp, vmix(r, r)))
         end do
      else
         allocate(model%coeff_se(0, 0), model%coeff_v(0, 0))
      end if
      model%has_se = .true.
   end subroutine polca_standard_errors

   subroutine compute_fit_statistics(y, model)
      integer, intent(in) :: y(:, :) !! Fitted categorical response matrix with zero indicating missing data.
      type(polca_model), intent(inout) :: model !! Fitted model receiving Pearson and deviance statistics.
      integer, allocatable :: cells(:, :), freq(:)
      real(dp), allocatable :: expected(:)
      integer :: i

      model%nobs_complete = count(all(y > 0, dim=2))
      if (model%nobs_complete == 0) then
         model%chisq = 0.0_dp
         model%gsq = 0.0_dp
         return
      end if
      call compress_complete(y, cells, freq)
      allocate(expected(size(freq)))
      call expected_cell_counts(model, cells, model%nobs_complete, expected)
      model%chisq = 0.0_dp
      model%gsq = 0.0_dp
      do i = 1, size(freq)
         if (expected(i) > tiny(1.0_dp)) then
            model%chisq = model%chisq + (real(freq(i), dp) - expected(i))**2 / expected(i)
            if (freq(i) > 0) then
               model%gsq = model%gsq + 2.0_dp * real(freq(i), dp) &
                    * log(real(freq(i), dp) / expected(i))
            end if
         end if
      end do
      model%chisq = model%chisq + real(model%nobs_complete, dp) - sum(expected)
   end subroutine compute_fit_statistics

   subroutine compress_complete(y, cells, freq)
      integer, intent(in) :: y(:, :) !! Response matrix whose fully observed rows are compressed into unique cells.
      integer, allocatable, intent(out) :: cells(:, :) !! Unique complete response patterns.
      integer, allocatable, intent(out) :: freq(:) !! Frequency for each returned response pattern.
      integer, allocatable :: tmp(:, :), tf(:)
      integer :: i, j, nc
      logical :: found

      allocate(tmp(size(y, 1), size(y, 2)), tf(size(y, 1)))
      nc = 0
      tf = 0
      do i = 1, size(y, 1)
         if (any(y(i, :) == 0)) cycle
         found = .false.
         do j = 1, nc
            if (all(tmp(j, :) == y(i, :))) then
               tf(j) = tf(j) + 1
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            nc = nc + 1
            tmp(nc, :) = y(i, :)
            tf(nc) = 1
         end if
      end do
      allocate(cells(nc, size(y, 2)), freq(nc))
      cells = tmp(1:nc, :)
      freq = tf(1:nc)
   end subroutine compress_complete

   subroutine expected_cell_counts(model, cells, nobs, expected)
      type(polca_model), intent(in) :: model !! Fitted latent-class model providing response probabilities and shares.
      integer, intent(in) :: cells(:, :) !! Complete response patterns to score.
      integer, intent(in) :: nobs !! Number of observations scaling the predicted cell probabilities.
      real(dp), intent(out) :: expected(:) !! Expected counts for the supplied complete response patterns.
      real(dp) :: p
      integer :: i, j, r

      do i = 1, size(cells, 1)
         expected(i) = 0.0_dp
         do r = 1, model%n_classes
            p = model%class_share(r)
            do j = 1, model%n_items
               p = p * model%probs(r, cells(i, j), j)
            end do
            expected(i) = expected(i) + p
         end do
         expected(i) = real(nobs, dp) * expected(i)
      end do
   end subroutine expected_cell_counts

end module polca_core
