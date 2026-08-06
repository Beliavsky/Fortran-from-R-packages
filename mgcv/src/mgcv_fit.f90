! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_fit
   use mgcv_kinds, only : dp
   use mgcv_linalg, only : spd_solve, spd_inverse, trace_product, logdet_spd, jacobi_eigen
   use mgcv_families, only : family_t, family_gaussian, family_binomial, &
      link_inverse, link_function, mu_eta, variance_function, initialize_mu, &
      deviance_sum, family_name
   implicit none
   private

   integer, parameter, public :: method_fixed = 0
   integer, parameter, public :: method_gcv = 1
   integer, parameter, public :: method_reml = 2
   integer, parameter, public :: method_ubre = 3

   type, public :: gam_control_t
      integer :: max_irls = 100
      integer :: max_outer = 40
      real(dp) :: irls_tolerance = 1.0e-8_dp
      real(dp) :: outer_tolerance = 1.0e-4_dp
      real(dp) :: initial_log_lambda = 0.0_dp
      real(dp) :: initial_step = 3.0_dp
      logical :: trace = .false.
   end type gam_control_t

   type, public :: gam_model_t
      type(family_t) :: family
      integer :: method = method_fixed
      integer :: iterations = 0
      integer :: outer_iterations = 0
      logical :: converged = .false.
      integer :: rank = 0
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: scale = 1.0_dp
      real(dp) :: edf = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: linear_predictor(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: leverage(:)
      real(dp), allocatable :: lambda(:)
   contains
      procedure :: summary => print_gam_summary
   end type gam_model_t

   public :: gam_fit, magic_fit, predict_gam, fitted_values

contains

   subroutine gam_fit(x, y, penalties, model, status, family, lambda, method, weights, offset, control)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), intent(in), optional :: penalties(:, :, :)
      type(gam_model_t), intent(out) :: model
      integer, intent(out) :: status
      type(family_t), intent(in), optional :: family
      real(dp), intent(in), optional :: lambda(:), weights(:), offset(:)
      integer, intent(in), optional :: method
      type(gam_control_t), intent(in), optional :: control
      type(family_t) :: fam
      type(gam_control_t) :: ctrl
      real(dp), allocatable :: w(:), off(:), lam(:), loglam(:), candidate(:)
      real(dp), allocatable :: s(:, :, :)
      type(gam_model_t) :: trial, best
      integer :: meth, m, i, outer
      real(dp) :: step, best_obj, trial_obj, improvement
      logical :: changed

      status = 0
      fam = family_t(); if (present(family)) fam = family
      ctrl = gam_control_t(); if (present(control)) ctrl = control
      meth = method_gcv; if (present(method)) meth = method
      if (size(x, 1) /= size(y) .or. size(y) < 2 .or. size(x, 2) < 1) then
         status = 1; return
      end if
      allocate(w(size(y)), off(size(y)))
      w = 1.0_dp; if (present(weights)) then
         if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) then; status = 2; return; end if
         w = weights
      end if
      off = 0.0_dp; if (present(offset)) then
         if (size(offset) /= size(y)) then; status = 3; return; end if
         off = offset
      end if
      if (present(penalties)) then
         if (size(penalties, 1) /= size(x, 2) .or. size(penalties, 2) /= size(x, 2)) then
            status = 4; return
         end if
         allocate(s(size(penalties, 1), size(penalties, 2), size(penalties, 3)))
         s = penalties; m = size(s, 3)
      else
         allocate(s(size(x, 2), size(x, 2), 0)); m = 0
      end if
      allocate(lam(m))
      if (present(lambda)) then
         if (size(lambda) /= m .or. any(lambda < 0.0_dp)) then; status = 5; return; end if
         lam = lambda
      else
         lam = exp(ctrl%initial_log_lambda)
      end if

      if (meth == method_fixed .or. m == 0 .or. present(lambda)) then
         call fit_fixed_lambda(x, y, s, lam, fam, w, off, ctrl, model, status)
         if (status == 0) then
            model%method = meth
            model%objective = model_objective(model, x, s, lam, w, meth)
         end if
         return
      end if

      allocate(loglam(m), candidate(m))
      loglam = ctrl%initial_log_lambda
      lam = exp(loglam)
      call fit_fixed_lambda(x, y, s, lam, fam, w, off, ctrl, best, status)
      if (status /= 0) return
      best_obj = model_objective(best, x, s, lam, w, meth)
      step = ctrl%initial_step
      outer = 0
      do while (outer < ctrl%max_outer .and. step > ctrl%outer_tolerance)
         outer = outer + 1; changed = .false.; improvement = 0.0_dp
         do i = 1, m
            candidate = loglam; candidate(i) = min(25.0_dp, loglam(i) + step)
            lam = exp(candidate)
            call fit_fixed_lambda(x, y, s, lam, fam, w, off, ctrl, trial, status)
            if (status == 0) then
               trial_obj = model_objective(trial, x, s, lam, w, meth)
               if (trial_obj < best_obj) then
                  improvement = max(improvement, best_obj - trial_obj)
                  best = trial; best_obj = trial_obj; loglam = candidate; changed = .true.
               end if
            end if
            candidate = loglam; candidate(i) = max(-25.0_dp, loglam(i) - step)
            lam = exp(candidate)
            call fit_fixed_lambda(x, y, s, lam, fam, w, off, ctrl, trial, status)
            if (status == 0) then
               trial_obj = model_objective(trial, x, s, lam, w, meth)
               if (trial_obj < best_obj) then
                  improvement = max(improvement, best_obj - trial_obj)
                  best = trial; best_obj = trial_obj; loglam = candidate; changed = .true.
               end if
            end if
         end do
         if (.not. changed) step = 0.5_dp * step
         if (ctrl%trace) write(*,'(a,i0,a,es12.4,a,es12.4)') &
            'outer=', outer, ' objective=', best_obj, ' step=', step
         if (changed .and. improvement <= ctrl%outer_tolerance * max(1.0_dp, abs(best_obj))) step = 0.5_dp * step
      end do
      model = best
      model%method = meth; model%outer_iterations = outer; model%objective = best_obj
      model%lambda = exp(loglam)
      status = 0
   end subroutine gam_fit

   subroutine magic_fit(x, y, penalties, model, status, lambda, method, weights, control)
      real(dp), intent(in) :: x(:, :), y(:), penalties(:, :, :)
      type(gam_model_t), intent(out) :: model
      integer, intent(out) :: status
      real(dp), intent(in), optional :: lambda(:), weights(:)
      integer, intent(in), optional :: method
      type(gam_control_t), intent(in), optional :: control
      type(family_t) :: gaussian
      gaussian%id = family_gaussian
      call gam_fit(x, y, penalties, model, status, gaussian, lambda, method, weights=weights, control=control)
   end subroutine magic_fit

   subroutine fit_fixed_lambda(x, y, penalties, lambda, family, weights, offset, control, model, status)
      real(dp), intent(in) :: x(:, :), y(:), penalties(:, :, :), lambda(:), weights(:), offset(:)
      type(family_t), intent(in) :: family
      type(gam_control_t), intent(in) :: control
      type(gam_model_t), intent(out) :: model
      integer, intent(out) :: status
      real(dp), allocatable :: mu(:), eta(:), z(:), work_w(:), beta(:), beta_old(:)
      real(dp), allocatable :: penalty(:, :), xtwx(:, :), a(:, :), rhs(:, :), solution(:, :)
      real(dp), allocatable :: ainv(:, :), wx(:, :), old_mu(:)
      real(dp) :: dev, old_dev, penalized, old_penalized, alpha, dmu
      integer :: i, j, iter, p, n
      logical :: converged

      status = 0; n = size(y); p = size(x, 2)
      allocate(penalty(p, p)); penalty = 0.0_dp
      do j = 1, size(lambda)
         penalty = penalty + lambda(j) * penalties(:, :, j)
      end do
      allocate(mu(n), eta(n), z(n), work_w(n), beta(p), beta_old(p), old_mu(n))
      call initialize_mu(y, family, mu)
      eta = link_function(mu, family) - offset
      beta = 0.0_dp; old_dev = huge(1.0_dp); old_penalized = huge(1.0_dp); converged = .false.

      do iter = 1, control%max_irls
         beta_old = beta; old_mu = mu
         do i = 1, n
            dmu = mu_eta(eta(i) + offset(i), family)
            work_w(i) = weights(i) * dmu * dmu / max(variance_function(mu(i), family), 1.0e-14_dp)
            z(i) = eta(i) + (y(i) - mu(i)) / max(dmu, 1.0e-12_dp)
         end do
         call weighted_crossproducts(x, z, work_w, xtwx, rhs)
         allocate(a(p, p)); a = xtwx + penalty
         do i = 1, p
            a(i, i) = a(i, i) + 1.0e-10_dp * max(1.0_dp, abs(xtwx(i, i)))
         end do
         call spd_solve(a, rhs, solution, status, 1.0e-12_dp)
         if (status /= 0) return
         beta = solution(:, 1)
         eta = matmul(x, beta)
         mu = link_inverse(eta + offset, family)
         dev = deviance_sum(y, mu, weights, family)
         penalized = dev + dot_product(beta, matmul(penalty, beta))
         if (iter > 1 .and. penalized > old_penalized * (1.0_dp + 1.0e-10_dp)) then
            alpha = 0.5_dp
            do while (alpha > 1.0e-6_dp)
               beta = beta_old + alpha * (beta - beta_old)
               eta = matmul(x, beta); mu = link_inverse(eta + offset, family)
               dev = deviance_sum(y, mu, weights, family)
               penalized = dev + dot_product(beta, matmul(penalty, beta))
               if (penalized <= old_penalized) exit
               alpha = 0.5_dp * alpha
            end do
         end if
         if (abs(old_dev - dev) <= control%irls_tolerance * (1.0_dp + abs(dev))) then
            converged = .true.; exit
         end if
         old_dev = dev; old_penalized = penalized
         deallocate(xtwx, rhs, a, solution)
      end do
      if (allocated(xtwx)) deallocate(xtwx)
      if (allocated(rhs)) deallocate(rhs)
      if (allocated(a)) deallocate(a)
      if (allocated(solution)) deallocate(solution)
      call weighted_crossproducts(x, z, work_w, xtwx, rhs)
      allocate(a(p, p)); a = xtwx + penalty
      do i = 1, p
         a(i, i) = a(i, i) + 1.0e-10_dp * max(1.0_dp, abs(xtwx(i, i)))
      end do
      call spd_inverse(a, ainv, status, 1.0e-12_dp)
      if (status /= 0) return

      model%family = family; model%iterations = iter; model%converged = converged
      model%deviance = deviance_sum(y, mu, weights, family)
      model%edf = trace_product(ainv, xtwx)
      model%rank = count([(a(i, i) > 1.0e-10_dp, i=1,p)])
      if (family%id == family_gaussian) then
         model%scale = model%deviance / max(1.0_dp, real(count(weights > 0.0_dp), dp) - model%edf)
      else
         model%scale = 1.0_dp
      end if
      allocate(model%coefficients(p), model%covariance(p, p), model%fitted(n), &
               model%linear_predictor(n), model%residuals(n), model%leverage(n), model%lambda(size(lambda)))
      model%coefficients = beta; model%covariance = model%scale * ainv
      model%fitted = mu; model%linear_predictor = eta + offset
      model%lambda = lambda
      do i = 1, n
         model%residuals(i) = sign(sqrt(max(0.0_dp, unit_deviance(y(i), mu(i), family))), y(i) - mu(i))
      end do
      allocate(wx(n, p))
      wx = x * spread(sqrt(max(work_w, 0.0_dp)), 2, p)
      do i = 1, n
         model%leverage(i) = work_w(i) * dot_product(x(i, :), matmul(ainv, x(i, :)))
      end do
   end subroutine fit_fixed_lambda

   subroutine weighted_crossproducts(x, z, weights, xtwx, rhs)
      real(dp), intent(in) :: x(:, :), z(:), weights(:)
      real(dp), allocatable, intent(out) :: xtwx(:, :), rhs(:, :)
      real(dp), allocatable :: wx(:, :)
      integer :: p
      p = size(x, 2)
      allocate(wx(size(x, 1), p), xtwx(p, p), rhs(p, 1))
      wx = x * spread(weights, 2, p)
      xtwx = matmul(transpose(x), wx)
      rhs(:, 1) = matmul(transpose(x), weights * z)
   end subroutine weighted_crossproducts

   function model_objective(model, x, penalties, lambda, weights, method) result(value)
      type(gam_model_t), intent(in) :: model
      real(dp), intent(in) :: x(:, :), penalties(:, :, :), lambda(:), weights(:)
      integer, intent(in) :: method
      real(dp) :: value, n_eff, denominator, logdet_a, logdet_p, scale
      real(dp), allocatable :: xtx(:, :), pmat(:, :), a(:, :), vals(:), vecs(:, :)
      integer :: i, j, status, rankp
      n_eff = real(count(weights > 0.0_dp), dp)
      denominator = max(1.0_dp, n_eff - model%edf)
      select case (method)
      case (method_ubre)
         scale = max(model%scale, 1.0e-12_dp)
         value = model%deviance / n_eff - 2.0_dp * scale * denominator / n_eff + scale
      case (method_reml)
         allocate(xtx(size(x, 2), size(x, 2)), pmat(size(x, 2), size(x, 2)))
         xtx = matmul(transpose(x), x); pmat = 0.0_dp
         do j = 1, size(lambda); pmat = pmat + lambda(j) * penalties(:, :, j); end do
         a = xtx + pmat
         do i = 1, size(a, 1); a(i, i) = a(i, i) + 1.0e-10_dp; end do
         logdet_a = logdet_spd(a, status, 1.0e-12_dp)
         call jacobi_eigen(pmat, vals, vecs, status)
         logdet_p = 0.0_dp; rankp = 0
         if (status == 0 .and. size(vals) > 0) then
            do i = 1, size(vals)
               if (vals(i) > 1.0e-10_dp * max(1.0_dp, vals(1))) then
                  logdet_p = logdet_p + log(vals(i)); rankp = rankp + 1
               end if
            end do
         end if
         scale = max(model%deviance / denominator, 1.0e-14_dp)
         value = 0.5_dp * (denominator * log(scale) + logdet_a - logdet_p)
      case default
         value = n_eff * model%deviance / (denominator * denominator)
      end select
      if (.not. (value < huge(1.0_dp))) value = huge(1.0_dp)
   end function model_objective

   elemental function unit_deviance(y, mu, family) result(term)
      real(dp), intent(in) :: y, mu
      type(family_t), intent(in) :: family
      real(dp) :: term, m
      m = max(mu, 1.0e-12_dp)
      select case (family%id)
      case (family_gaussian)
         term = (y - m)**2
      case (family_binomial)
         term = 0.0_dp
         if (y > 0.0_dp) term = term + y * log(y / m)
         if (y < 1.0_dp) term = term + (1.0_dp - y) * log((1.0_dp - y) / (1.0_dp - m))
         term = 2.0_dp * term
      case default
         if (y <= tiny(1.0_dp)) then
            term = 2.0_dp * m
         else
            term = 2.0_dp * (y * log(max(y, 1.0e-12_dp) / m) - (y - m))
         end if
      end select
   end function unit_deviance

   subroutine predict_gam(model, x, eta, response, status, offset, se)
      type(gam_model_t), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: eta(:), response(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: offset(:)
      real(dp), allocatable, intent(out), optional :: se(:)
      real(dp), allocatable :: off(:)
      integer :: i
      if (.not. allocated(model%coefficients) .or. size(x, 2) /= size(model%coefficients)) then
         allocate(eta(0), response(0)); if (present(se)) allocate(se(0)); status = 1; return
      end if
      allocate(off(size(x, 1))); off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= size(x, 1)) then
            allocate(eta(0), response(0)); if (present(se)) allocate(se(0)); status = 2; return
         end if
         off = offset
      end if
      allocate(eta(size(x, 1)), response(size(x, 1)))
      eta = matmul(x, model%coefficients) + off
      response = link_inverse(eta, model%family)
      if (present(se)) then
         allocate(se(size(x, 1)))
         do i = 1, size(x, 1)
            se(i) = sqrt(max(0.0_dp, dot_product(x(i, :), matmul(model%covariance, x(i, :)))))
         end do
      end if
      status = 0
   end subroutine predict_gam

   function fitted_values(model) result(values)
      type(gam_model_t), intent(in) :: model
      real(dp), allocatable :: values(:)
      if (allocated(model%fitted)) then
         allocate(values(size(model%fitted))); values = model%fitted
      else
         allocate(values(0))
      end if
   end function fitted_values

   subroutine print_gam_summary(self)
      class(gam_model_t), intent(in) :: self
      character(len=16) :: method_name
      select case (self%method)
      case (method_fixed); method_name = 'fixed'
      case (method_gcv); method_name = 'GCV'
      case (method_reml); method_name = 'REML-like'
      case (method_ubre); method_name = 'UBRE'
      case default; method_name = 'unknown'
      end select
      write(*,'(a)') 'Generalized additive model'
      write(*,'(a,a)') '  family: ', family_name(self%family)
      write(*,'(a,a)') '  method: ', trim(method_name)
      write(*,'(a,l1)') '  converged: ', self%converged
      write(*,'(a,i0)') '  PIRLS iterations: ', self%iterations
      write(*,'(a,f12.5)') '  deviance: ', self%deviance
      write(*,'(a,f12.5)') '  effective df: ', self%edf
      write(*,'(a,f12.5)') '  scale: ', self%scale
      if (allocated(self%lambda)) write(*,'(a,*(es12.4,1x))') '  smoothing parameters: ', self%lambda
   end subroutine print_gam_summary

end module mgcv_fit
