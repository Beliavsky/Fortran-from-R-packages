module gamm4_fit_mod
   use gamm4_kinds, only : dp
   use gamm4_types, only : gamm4_smooth_t, gamm4_control_t, gamm4_result_t, &
      gamm4_vb_result_t, gamm4_family_gaussian
   use gamm4_reparam, only : smooth_reparam_t, reparameterize_smooth
   use gamm4_covariance, only : gamm4_get_vb
   use lme4, only : random_term_t, family_binomial, family_poisson, family_gamma, &
      family_inverse_gaussian, family_negative_binomial, covariance_unstructured, &
      covariance_diagonal, covariance_compound_symmetry, covariance_ar1, &
      build_random_design, term_covariance_from_eta
   use r_linalg, only : spd_inverse_logdet
   implicit none
   private
   public :: gamm4_fit

contains

   subroutine gamm4_fit(y, x_fixed, smooths, result, status, random_terms, family, &
      weights, offset, reml, dispersion, control)
      real(dp), intent(in) :: y(:) !! Response vector with one value per observation.
      real(dp), intent(in) :: x_fixed(:, :) !! Parametric fixed-effect design; rows correspond to observations.
      type(gamm4_smooth_t), intent(in) :: smooths(:) !! Smooth terms represented by mgcv-compatible basis and penalty metadata.
      type(gamm4_result_t), intent(out) :: result !! Fitted coefficients, covariance, means, smoothing parameters, and diagnostics.
      integer, intent(out) :: status !! Zero on convergence; nonzero for validation, support, or numerical failures.
      type(random_term_t), intent(in), optional :: random_terms(:) !! Ordinary grouped random effects in lme4 numeric-array form.
      integer, intent(in), optional :: family !! Zero selects Gaussian; lme4 family constants select supported GLMM families.
      real(dp), intent(in), optional :: weights(:) !! Positive prior weights; defaults to one for every observation.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset; defaults to zero.
      logical, intent(in), optional :: reml !! For Gaussian fits, use REML when true and ML when false; defaults to true.
      real(dp), intent(in), optional :: dispersion !! Fixed GLMM dispersion for Gamma/inverse-Gaussian families; defaults to one.
      type(gamm4_control_t), intent(in), optional :: control !! Iteration limits, tolerances, and variance-parameter bounds.
      type(gamm4_control_t) :: ctrl
      type(smooth_reparam_t), allocatable :: reps(:)
      type(random_term_t), allocatable :: terms(:)
      real(dp), allocatable :: w(:), off(:), xfit(:, :), zs(:, :), zr(:, :), zall(:, :)
      real(dp), allocatable :: eta(:), beta(:), u(:), work_w(:), gre(:, :), gall(:, :)
      integer, allocatable :: offsets(:), fixed_start(:), random_start(:), eta_s(:), smooth_eta_col(:)
      integer :: fam, i, info, n, ns, p0, pfit, qs, qr, nt_re, nt, evals
      logical :: use_reml, ok
      real(dp) :: disp, objective, scale

      call initialize_result(result)
      status = 0
      n = size(y)
      ns = size(smooths)
      p0 = size(x_fixed, 2)
      if (n < 1 .or. size(x_fixed, 1) /= n) then
         call fail_result(result, status, 1, 'response and fixed-design dimensions are inconsistent')
         return
      end if
      allocate(w(n), off(n))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            call fail_result(result, status, 2, 'weights must be positive and match the response length')
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            call fail_result(result, status, 3, 'offset must match the response length')
            return
         end if
         off = offset
      end if
      fam = gamm4_family_gaussian
      if (present(family)) fam = family
      if (.not. supported_family(fam)) then
         call fail_result(result, status, 4, 'unsupported response-family code')
         return
      end if
      if (fam == family_binomial .and. (any(y < 0.0_dp) .or. any(y > 1.0_dp))) then
         call fail_result(result, status, 5, 'binomial responses must lie in [0,1]')
         return
      end if
      if ((fam == family_poisson .or. fam == family_gamma .or. fam == family_inverse_gaussian .or. &
           fam == family_negative_binomial) .and. any(y < 0.0_dp)) then
         call fail_result(result, status, 6, 'non-Gaussian positive-mean family received a negative response')
         return
      end if
      if ((fam == family_gamma .or. fam == family_inverse_gaussian) .and. any(y <= 0.0_dp)) then
         call fail_result(result, status, 7, 'Gamma and inverse-Gaussian responses must be strictly positive')
         return
      end if
      ctrl = gamm4_control_t()
      if (present(control)) ctrl = control
      use_reml = .true.
      if (present(reml)) use_reml = reml
      disp = 1.0_dp
      if (present(dispersion)) disp = dispersion
      if (disp <= 0.0_dp) then
         call fail_result(result, status, 8, 'dispersion must be positive')
         return
      end if

      allocate(reps(ns), fixed_start(ns), random_start(ns), eta_s(ns))
      pfit = p0
      qs = 0
      nt = 0
      do i = 1, ns
         call reparameterize_smooth(smooths(i), reps(i), ctrl%eigen_tolerance, info)
         if (info /= 0) then
            call fail_result(result, status, 20 + info, &
               'smooth reparameterization failed, multi-penalty smooths are not yet supported')
            return
         end if
         if (size(smooths(i)%basis, 1) /= n) then
            call fail_result(result, status, 25, 'smooth basis row count must match the response length')
            return
         end if
         fixed_start(i) = pfit + 1
         pfit = pfit + reps(i)%nullity
         random_start(i) = qs + 1
         qs = qs + reps(i)%rank
         if (reps(i)%rank > 0) then
            nt = nt + 1
            eta_s(i) = nt
         else
            eta_s(i) = 0
         end if
      end do
      allocate(xfit(n, pfit), zs(n, qs), smooth_eta_col(qs))
      xfit = 0.0_dp
      zs = 0.0_dp
      smooth_eta_col = 0
      if (p0 > 0) xfit(:, 1:p0) = x_fixed
      do i = 1, ns
         if (reps(i)%nullity > 0) then
            xfit(:, fixed_start(i):fixed_start(i) + reps(i)%nullity - 1) = reps(i)%fixed_design
         end if
         if (reps(i)%rank > 0) then
            zs(:, random_start(i):random_start(i) + reps(i)%rank - 1) = reps(i)%random_design
            smooth_eta_col(random_start(i):random_start(i) + reps(i)%rank - 1) = eta_s(i)
         end if
      end do

      if (present(random_terms)) then
         allocate(terms(size(random_terms)))
         terms = random_terms
         do i = 1, size(terms)
            call terms(i)%validate(n, ok, result%message)
            if (.not. ok) then
               call fail_result(result, status, 30, 'invalid lme4 random-effect term')
               return
            end if
         end do
         call build_random_design(terms, zr, offsets)
      else
         allocate(terms(0), zr(n, 0), offsets(1))
         offsets(1) = 1
      end if
      qr = size(zr, 2)
      nt_re = random_parameter_count(terms)
      nt = nt + nt_re
      allocate(zall(n, qs + qr))
      if (qs > 0) zall(:, 1:qs) = zs
      if (qr > 0) zall(:, qs + 1:qs + qr) = zr
      allocate(eta(nt))
      call initialize_eta(eta, eta_s, terms, ns)

      call optimize_variances(y, xfit, zall, smooth_eta_col, terms, eta, w, off, fam, use_reml, disp, ctrl, &
         objective, beta, u, work_w, scale, gre, gall, evals, info)
      if (info /= 0) then
         call fail_result(result, status, 40 + info, 'joint smooth/random-effect variance optimization failed')
         return
      end if
      result%iterations = evals
      call assemble_result(y, x_fixed, smooths, reps, xfit, zs, zr, eta, eta_s, w, off, fam, &
         use_reml, beta, u, work_w, scale, gre, objective, result, info)
      if (info /= 0) then
         call fail_result(result, status, 60 + info, 'fit converged but GAMM covariance/result assembly failed')
         return
      end if
      result%status = 0
      result%converged = .true.
      result%message = 'ok'
      status = 0
   end subroutine gamm4_fit

   subroutine optimize_variances(y, x, z, smooth_eta_col, terms, eta, weights, offset, family, reml, dispersion, &
      control, objective, beta, u, work_w, scale, gre, gall, evaluations, status)
      real(dp), intent(in) :: y(:) !! Response vector used in the joint mixed-model objective.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design after smooth null-space reparameterization.
      real(dp), intent(in) :: z(:, :) !! Combined smooth and ordinary random-effect design matrix.
      integer, intent(in) :: smooth_eta_col(:) !! Variance-parameter index for each leading smooth random-effect column.
      type(random_term_t), intent(in) :: terms(:) !! Ordinary lme4 random-effect terms whose covariance follows the smooth blocks.
      real(dp), intent(inout) :: eta(:) !! Log-standard-deviation/covariance parameters optimized in place.
      real(dp), intent(in) :: weights(:) !! Positive prior weights entering the Gaussian or PIRLS objective.
      real(dp), intent(in) :: offset(:) !! Additive linear-predictor offset.
      integer, intent(in) :: family !! Gaussian or supported lme4 family code.
      logical, intent(in) :: reml !! Gaussian REML flag; ignored for non-Gaussian families.
      real(dp), intent(in) :: dispersion !! Fixed non-Gaussian dispersion where relevant.
      type(gamm4_control_t), intent(in) :: control !! Variance-optimization and PIRLS controls.
      real(dp), intent(out) :: objective !! Final negative twice log-likelihood/deviance criterion.
      real(dp), allocatable, intent(out) :: beta(:) !! Final fixed-effect coefficient vector.
      real(dp), allocatable, intent(out) :: u(:) !! Final random-effect mode vector.
      real(dp), allocatable, intent(out) :: work_w(:) !! Final PIRLS working weights; Gaussian fits return the prior weights.
      real(dp), intent(out) :: scale !! Gaussian residual variance or supplied non-Gaussian dispersion.
      real(dp), allocatable, intent(out) :: gre(:, :) !! Relative covariance matrix for ordinary, non-smooth random effects.
      real(dp), allocatable, intent(out) :: gall(:, :) !! Relative covariance matrix for all smooth and ordinary random effects.
      integer, intent(out) :: evaluations !! Number of objective evaluations made by coordinate optimization.
      integer, intent(out) :: status !! Zero on convergence; nonzero on numerical failure.
      real(dp), allocatable :: trial(:), best_beta(:), best_u(:), best_w(:), best_gre(:, :), best_gall(:, :)
      real(dp) :: best, candidate, old_best, lower, upper, x1, x2, f1, f2, phi, tol
      real(dp) :: candidate_scale
      integer :: i, iter, info

      status = 0
      evaluations = 0
      call evaluate_model(y, x, z, smooth_eta_col, terms, eta, weights, offset, family, reml, dispersion, control, &
         best, best_beta, best_u, best_w, scale, best_gre, best_gall, info)
      evaluations = evaluations + 1
      if (info /= 0) then
         status = 1
         return
      end if
      allocate(trial(size(eta)))
      phi = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)
      tol = max(1.0e-7_dp, control%tolerance)
      do iter = 1, control%max_outer
         old_best = best
         do i = 1, size(eta)
            call eta_bounds(i, size(smooth_eta_col), terms, eta, lower, upper, control)
            x1 = upper - phi * (upper - lower)
            x2 = lower + phi * (upper - lower)
            trial = eta
            trial(i) = x1
            call evaluate_model(y, x, z, smooth_eta_col, terms, trial, weights, offset, family, reml, dispersion, control, &
               f1, beta, u, work_w, candidate_scale, gre, gall, info)
            evaluations = evaluations + 1
            if (info /= 0) f1 = huge(1.0_dp)
            trial(i) = x2
            call evaluate_model(y, x, z, smooth_eta_col, terms, trial, weights, offset, family, reml, dispersion, control, &
               f2, beta, u, work_w, candidate_scale, gre, gall, info)
            evaluations = evaluations + 1
            if (info /= 0) f2 = huge(1.0_dp)
            do while (upper - lower > tol .and. evaluations < 100000)
               if (f1 > f2) then
                  lower = x1
                  x1 = x2
                  f1 = f2
                  x2 = lower + phi * (upper - lower)
                  trial = eta
                  trial(i) = x2
                  call evaluate_model(y, x, z, smooth_eta_col, terms, trial, weights, offset, family, reml, dispersion, control, &
                     f2, beta, u, work_w, candidate_scale, gre, gall, info)
                  evaluations = evaluations + 1
                  if (info /= 0) f2 = huge(1.0_dp)
               else
                  upper = x2
                  x2 = x1
                  f2 = f1
                  x1 = upper - phi * (upper - lower)
                  trial = eta
                  trial(i) = x1
                  call evaluate_model(y, x, z, smooth_eta_col, terms, trial, weights, offset, family, reml, dispersion, control, &
                     f1, beta, u, work_w, candidate_scale, gre, gall, info)
                  evaluations = evaluations + 1
                  if (info /= 0) f1 = huge(1.0_dp)
               end if
               if (abs(upper - lower) <= max(tol, 1.0e-3_dp)) exit
            end do
            trial = eta
            trial(i) = 0.5_dp * (lower + upper)
            call evaluate_model(y, x, z, smooth_eta_col, terms, trial, weights, offset, family, reml, dispersion, control, &
               candidate, beta, u, work_w, candidate_scale, gre, gall, info)
            evaluations = evaluations + 1
            if (info == 0 .and. candidate < best) then
               eta = trial
               best = candidate
               best_beta = beta
               best_u = u
               best_w = work_w
               scale = candidate_scale
               best_gre = gre
               best_gall = gall
            end if
         end do
         if (abs(old_best - best) <= control%tolerance * (1.0_dp + abs(best))) exit
      end do
      objective = best
      beta = best_beta
      u = best_u
      work_w = best_w
      gre = best_gre
      gall = best_gall
   end subroutine optimize_variances

   subroutine evaluate_model(y, x, z, smooth_eta_col, terms, eta, weights, offset, family, reml, dispersion, &
      control, objective, beta, u, work_w, scale, gre, gall, status)
      real(dp), intent(in) :: y(:) !! Response vector for the current objective evaluation.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix.
      real(dp), intent(in) :: z(:, :) !! Combined smooth and ordinary random-effect design matrix.
      integer, intent(in) :: smooth_eta_col(:) !! Variance-parameter index for each leading smooth random-effect column.
      type(random_term_t), intent(in) :: terms(:) !! Ordinary random-effect terms.
      real(dp), intent(in) :: eta(:) !! Current variance/covariance parameter vector.
      real(dp), intent(in) :: weights(:) !! Positive prior observation weights.
      real(dp), intent(in) :: offset(:) !! Additive linear-predictor offset.
      integer, intent(in) :: family !! Gaussian or supported lme4 family code.
      logical, intent(in) :: reml !! Gaussian REML flag.
      real(dp), intent(in) :: dispersion !! Fixed non-Gaussian dispersion parameter.
      type(gamm4_control_t), intent(in) :: control !! PIRLS iteration controls.
      real(dp), intent(out) :: objective !! Current Gaussian REML/ML or Laplace-PIRLS objective.
      real(dp), allocatable, intent(out) :: beta(:) !! Fixed-effect coefficients at this variance parameter.
      real(dp), allocatable, intent(out) :: u(:) !! Random-effect mode at this variance parameter.
      real(dp), allocatable, intent(out) :: work_w(:) !! Final working weights.
      real(dp), intent(out) :: scale !! Gaussian residual variance or supplied dispersion.
      real(dp), allocatable, intent(out) :: gre(:, :) !! Ordinary random-effect relative covariance matrix.
      real(dp), allocatable, intent(out) :: gall(:, :) !! Combined relative covariance matrix.
      integer, intent(out) :: status !! Zero on success; nonzero on invalid covariance or failed linear algebra.
      call build_combined_covariance(smooth_eta_col, terms, eta, gre, gall, status)
      if (status /= 0) return
      if (family == gamm4_family_gaussian) then
         call evaluate_gaussian(y, x, z, gall, weights, offset, reml, objective, beta, u, scale, status)
         work_w = weights
      else
         scale = dispersion
         call evaluate_glmm(y, x, z, gall, weights, offset, family, dispersion, control, objective, beta, u, work_w, status)
      end if
   end subroutine evaluate_model

   subroutine evaluate_gaussian(y, x, z, g, weights, offset, reml, objective, beta, u, scale, status)
      real(dp), intent(in) :: y(:) !! Gaussian response vector.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix.
      real(dp), intent(in) :: z(:, :) !! Random-effect design matrix.
      real(dp), intent(in) :: g(:, :) !! Relative random-effect covariance matrix.
      real(dp), intent(in) :: weights(:) !! Positive prior precision weights.
      real(dp), intent(in) :: offset(:) !! Additive fixed offset removed before GLS fitting.
      logical, intent(in) :: reml !! Use restricted likelihood when true.
      real(dp), intent(out) :: objective !! Negative twice profiled Gaussian log likelihood up to lme4 constants.
      real(dp), allocatable, intent(out) :: beta(:) !! Generalized least-squares fixed-effect coefficients.
      real(dp), allocatable, intent(out) :: u(:) !! Conditional random-effect modes on the relative scale.
      real(dp), intent(out) :: scale !! Profiled residual variance.
      integer, intent(out) :: status !! Zero on success; nonzero when required matrices are singular.
      real(dp), allocatable :: v(:, :), vinv(:, :), xtvx(:, :), xtvx_inv(:, :), r(:), centered(:)
      real(dp) :: logdetv, logdetx, rss, denom, pi
      integer :: i, info, n, p

      status = 0
      n = size(y)
      p = size(x, 2)
      allocate(v(n, n), centered(n), r(n))
      v = matmul(matmul(z, g), transpose(z))
      do i = 1, n
         v(i, i) = v(i, i) + 1.0_dp / weights(i)
      end do
      call spd_inverse_logdet(v, vinv, logdetv, info)
      if (info /= 0) then
         status = 1
         return
      end if
      centered = y - offset
      xtvx = matmul(transpose(x), matmul(vinv, x))
      call spd_inverse_logdet(xtvx, xtvx_inv, logdetx, info)
      if (info /= 0) then
         status = 2
         return
      end if
      beta = matmul(xtvx_inv, matmul(transpose(x), matmul(vinv, centered)))
      r = centered - matmul(x, beta)
      rss = dot_product(r, matmul(vinv, r))
      denom = real(n, dp)
      if (reml) denom = real(n - p, dp)
      if (denom <= 0.0_dp .or. rss <= 0.0_dp) then
         status = 3
         return
      end if
      scale = rss / denom
      pi = acos(-1.0_dp)
      objective = denom * (log(2.0_dp * pi) + 1.0_dp + log(scale)) + logdetv
      if (reml) objective = objective + logdetx
      if (size(z, 2) > 0) then
         u = matmul(g, matmul(transpose(z), matmul(vinv, r)))
      else
         allocate(u(0))
      end if
   end subroutine evaluate_gaussian

   subroutine evaluate_glmm(y, x, z, g, weights, offset, family, dispersion, control, &
      objective, beta, u, work_w, status)
      real(dp), intent(in) :: y(:) !! Non-Gaussian response vector.
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design matrix.
      real(dp), intent(in) :: z(:, :) !! Random-effect design matrix.
      real(dp), intent(in) :: g(:, :) !! Positive-definite random-effect covariance matrix.
      real(dp), intent(in) :: weights(:) !! Prior observation weights.
      real(dp), intent(in) :: offset(:) !! Additive linear-predictor offset.
      integer, intent(in) :: family !! Supported lme4 family code.
      real(dp), intent(in) :: dispersion !! Fixed Gamma/inverse-Gaussian dispersion; one for binomial/Poisson.
      type(gamm4_control_t), intent(in) :: control !! PIRLS iteration limit and tolerance.
      real(dp), intent(out) :: objective !! Laplace/PIRLS objective at the conditional mode.
      real(dp), allocatable, intent(out) :: beta(:) !! Fixed-effect conditional-mode coefficients.
      real(dp), allocatable, intent(out) :: u(:) !! Random-effect conditional modes.
      real(dp), allocatable, intent(out) :: work_w(:) !! Final PIRLS working weights.
      integer, intent(out) :: status !! Zero on convergence; nonzero for failed linear algebra or PIRLS convergence.
      real(dp), allocatable :: ginv(:, :), design(:, :), h(:, :), hinv(:, :), rhs(:), coef(:)
      real(dp), allocatable :: eta_lin(:), mu(:), zwork(:), oldcoef(:), wz(:)
      real(dp) :: logdetg, logdeth, penalty, cond_dev, ignored
      integer :: info, iter, n, p, q

      status = 0
      n = size(y)
      p = size(x, 2)
      q = size(z, 2)
      if (q < 1) then
         status = 1
         return
      end if
      call spd_inverse_logdet(g, ginv, logdetg, info)
      if (info /= 0) then
         status = 2
         return
      end if
      allocate(design(n, p + q), beta(p), u(q), work_w(n), eta_lin(n), mu(n), zwork(n), &
         oldcoef(p + q), coef(p + q), rhs(p + q), wz(n))
      design(:, 1:p) = x
      design(:, p + 1:p + q) = z
      beta = 0.0_dp
      u = 0.0_dp
      coef = 0.0_dp
      do iter = 1, control%max_pirls
         oldcoef = coef
         eta_lin = offset + matmul(x, beta) + matmul(z, u)
         call family_working_values(y, eta_lin, weights, family, dispersion, mu, work_w, zwork, status)
         if (status /= 0) return
         wz = work_w * (zwork - offset)
         h = matmul(transpose(design), spread(work_w, 2, p + q) * design)
         h(p + 1:p + q, p + 1:p + q) = h(p + 1:p + q, p + 1:p + q) + ginv
         rhs = matmul(transpose(design), wz)
         call spd_inverse_logdet(h, hinv, ignored, info)
         if (info /= 0) then
            status = 3
            return
         end if
         coef = matmul(hinv, rhs)
         beta = coef(1:p)
         u = coef(p + 1:p + q)
         if (maxval(abs(coef - oldcoef)) <= control%pirls_tolerance * (1.0_dp + maxval(abs(coef)))) exit
      end do
      if (iter > control%max_pirls) then
         status = 4
         return
      end if
      eta_lin = offset + matmul(x, beta) + matmul(z, u)
      call family_working_values(y, eta_lin, weights, family, dispersion, mu, work_w, zwork, status)
      if (status /= 0) return
      h = matmul(transpose(z), spread(work_w, 2, q) * z) + ginv
      call spd_inverse_logdet(h, hinv, logdeth, info)
      if (info /= 0) then
         status = 5
         return
      end if
      penalty = dot_product(u, matmul(ginv, u))
      cond_dev = conditional_deviance(y, mu, weights, family, dispersion)
      objective = cond_dev + penalty + logdetg + logdeth
   end subroutine evaluate_glmm

   subroutine family_working_values(y, eta, prior_w, family, dispersion, mu, work_w, zwork, status)
      real(dp), intent(in) :: y(:) !! Observed responses for one PIRLS update.
      real(dp), intent(in) :: eta(:) !! Current linear predictor including offset.
      real(dp), intent(in) :: prior_w(:) !! Prior observation weights.
      integer, intent(in) :: family !! Supported lme4 family code.
      real(dp), intent(in) :: dispersion !! Fixed response dispersion.
      real(dp), intent(out) :: mu(:) !! Mean response at the current linear predictor.
      real(dp), intent(out) :: work_w(:) !! PIRLS working weights.
      real(dp), intent(out) :: zwork(:) !! PIRLS working response on the link scale.
      integer, intent(out) :: status !! Zero on success; nonzero for unsupported family or nonfinite working quantities.
      real(dp) :: e, dmu, variance
      integer :: i

      status = 0
      do i = 1, size(y)
         e = max(-35.0_dp, min(35.0_dp, eta(i)))
         select case (family)
         case (family_binomial)
            mu(i) = 1.0_dp / (1.0_dp + exp(-e))
            mu(i) = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, mu(i)))
            dmu = mu(i) * (1.0_dp - mu(i))
            variance = dmu
         case (family_poisson)
            mu(i) = max(1.0e-12_dp, exp(e))
            dmu = mu(i)
            variance = mu(i)
         case (family_gamma)
            mu(i) = max(1.0e-12_dp, exp(e))
            dmu = mu(i)
            variance = dispersion * mu(i) * mu(i)
         case (family_inverse_gaussian)
            mu(i) = max(1.0e-12_dp, exp(e))
            dmu = mu(i)
            variance = dispersion * mu(i) ** 3
         case (family_negative_binomial)
            mu(i) = max(1.0e-12_dp, exp(e))
            dmu = mu(i)
            variance = mu(i) + dispersion * mu(i) * mu(i)
         case default
            status = 1
            return
         end select
         work_w(i) = prior_w(i) * dmu * dmu / max(variance, 1.0e-15_dp)
         zwork(i) = eta(i) + (y(i) - mu(i)) / max(dmu, 1.0e-15_dp)
      end do
   end subroutine family_working_values

   real(dp) function conditional_deviance(y, mu, weights, family, dispersion) result(value)
      real(dp), intent(in) :: y(:) !! Observed responses entering the family conditional deviance.
      real(dp), intent(in) :: mu(:) !! Fitted conditional means.
      real(dp), intent(in) :: weights(:) !! Prior observation weights.
      integer, intent(in) :: family !! Supported lme4 family code.
      real(dp), intent(in) :: dispersion !! Dispersion or negative-binomial overdispersion parameter.
      real(dp) :: yi, mui, term
      integer :: i

      value = 0.0_dp
      do i = 1, size(y)
         yi = y(i)
         mui = max(mu(i), 1.0e-15_dp)
         select case (family)
         case (family_binomial)
            term = 0.0_dp
            if (yi > 0.0_dp) term = term + yi * log(yi / mui)
            if (yi < 1.0_dp) term = term + (1.0_dp - yi) * log((1.0_dp - yi) / max(1.0_dp - mui, 1.0e-15_dp))
            value = value + 2.0_dp * weights(i) * term
         case (family_poisson)
            if (yi > 0.0_dp) then
               term = yi * log(yi / mui) - (yi - mui)
            else
               term = mui
            end if
            value = value + 2.0_dp * weights(i) * term
         case (family_gamma)
            value = value + 2.0_dp * weights(i) * ((yi - mui) / mui - log(yi / mui)) / dispersion
         case (family_inverse_gaussian)
            value = value + weights(i) * (yi - mui) ** 2 / (dispersion * yi * mui * mui)
         case (family_negative_binomial)
            if (yi > 0.0_dp) then
               term = yi * log(yi / mui)
            else
               term = 0.0_dp
            end if
            term = term - (yi + 1.0_dp / dispersion) * &
               log((yi + 1.0_dp / dispersion) / (mui + 1.0_dp / dispersion))
            value = value + 2.0_dp * weights(i) * term
         end select
      end do
   end function conditional_deviance

   subroutine build_combined_covariance(smooth_eta_col, terms, eta, gre, gall, status)
      integer, intent(in) :: smooth_eta_col(:) !! Variance-parameter index for every smooth random-effect design column.
      type(random_term_t), intent(in) :: terms(:) !! Ordinary random-effect terms following the smooth blocks.
      real(dp), intent(in) :: eta(:) !! Smooth log-standard-deviations followed by lme4 covariance parameters.
      real(dp), allocatable, intent(out) :: gre(:, :) !! Expanded relative covariance for ordinary random effects only.
      real(dp), allocatable, intent(out) :: gall(:, :) !! Expanded relative covariance for smooth and ordinary random effects.
      integer, intent(out) :: status !! Zero on success; nonzero when an lme4 covariance parameterization is invalid.
      real(dp), allocatable :: cov(:, :)
      integer :: i, info, j, np, p0, q, qr, r0, lev, qs, nsmooth

      status = 0
      qr = 0
      do i = 1, size(terms)
         qr = qr + terms(i)%n_random_effects()
      end do
      qs = size(smooth_eta_col)
      nsmooth = size(eta) - random_parameter_count(terms)
      allocate(gre(qr, qr), gall(qs + qr, qs + qr))
      gre = 0.0_dp
      gall = 0.0_dp
      do j = 1, qs
         if (smooth_eta_col(j) < 1 .or. smooth_eta_col(j) > nsmooth) then
            status = 2
            return
         end if
         gall(j, j) = exp(2.0_dp * eta(smooth_eta_col(j)))
      end do
      p0 = nsmooth
      r0 = 0
      do i = 1, size(terms)
         q = terms(i)%n_coefficients()
         np = terms(i)%n_parameters()
         call term_covariance_from_eta(terms(i), eta(p0 + 1:p0 + np), cov, info)
         if (info /= 0) then
            status = 1
            return
         end if
         do lev = 1, terms(i)%n_levels
            gre(r0 + (lev - 1) * q + 1:r0 + lev * q, r0 + (lev - 1) * q + 1:r0 + lev * q) = cov
         end do
         p0 = p0 + np
         r0 = r0 + terms(i)%n_random_effects()
      end do
      if (qr > 0) gall(qs + 1:qs + qr, qs + 1:qs + qr) = gre
   end subroutine build_combined_covariance

   subroutine assemble_result(y, x_fixed, smooths, reps, xfit, zs, zr, eta, eta_s, weights, offset, family, &
      reml, beta, u, work_w, scale, gre, objective, result, status)
      real(dp), intent(in) :: y(:) !! Observed response vector.
      real(dp), intent(in) :: x_fixed(:, :) !! Original parametric fixed-effect design.
      type(gamm4_smooth_t), intent(in) :: smooths(:) !! Original smooth basis/penalty objects.
      type(smooth_reparam_t), intent(in) :: reps(:) !! Fixed/random smooth reparameterizations used for fitting.
      real(dp), intent(in) :: xfit(:, :) !! Reparameterized fixed-effect design used by the mixed-model fit.
      real(dp), intent(in) :: zs(:, :) !! Smooth random-effect design matrix.
      real(dp), intent(in) :: zr(:, :) !! Ordinary random-effect design matrix.
      real(dp), intent(in) :: eta(:) !! Final variance/covariance parameter vector.
      integer, intent(in) :: eta_s(:) !! Index of each smooth variance parameter, zero for unpenalized smooths.
      real(dp), intent(in) :: weights(:) !! Prior observation weights.
      real(dp), intent(in) :: offset(:) !! Additive linear-predictor offset.
      integer, intent(in) :: family !! Gaussian or supported lme4 response-family code.
      logical, intent(in) :: reml !! Gaussian restricted-likelihood flag.
      real(dp), intent(in) :: beta(:) !! Reparameterized fixed-effect estimates.
      real(dp), intent(in) :: u(:) !! Combined smooth and ordinary random-effect modes.
      real(dp), intent(in) :: work_w(:) !! Final PIRLS working weights or Gaussian prior weights.
      real(dp), intent(in) :: scale !! Gaussian residual variance or non-Gaussian dispersion.
      real(dp), intent(in) :: gre(:, :) !! Relative covariance of ordinary random effects.
      real(dp), intent(in) :: objective !! Final mixed-model objective value.
      type(gamm4_result_t), intent(inout) :: result !! Result populated with original-basis GAM and conditional-fit quantities.
      integer, intent(out) :: status !! Zero on success; nonzero when covariance reconstruction fails.
      type(gamm4_vb_result_t) :: vbres
      real(dp), allocatable :: xorig(:, :), xfp(:, :), bmap(:, :), theta(:), spcol(:), v(:), phi(:, :)
      real(dp), allocatable :: gam_eta(:), cond_eta(:), mu(:)
      integer :: i, j, n, ns, p0, porig, col_orig, col_fit, us, qs, qr

      status = 0
      n = size(y)
      ns = size(smooths)
      p0 = size(x_fixed, 2)
      porig = p0
      do i = 1, ns
         porig = porig + size(smooths(i)%basis, 2)
      end do
      allocate(xorig(n, porig), xfp(n, porig), bmap(porig, porig), theta(porig), spcol(porig))
      xorig = 0.0_dp
      xfp = 0.0_dp
      bmap = 0.0_dp
      theta = 0.0_dp
      spcol = 0.0_dp
      if (p0 > 0) then
         xorig(:, 1:p0) = x_fixed
         xfp(:, 1:p0) = x_fixed
         do i = 1, p0
            bmap(i, i) = 1.0_dp
            theta(i) = beta(i)
         end do
      end if
      col_orig = p0 + 1
      col_fit = p0 + 1
      us = 1
      qs = size(zs, 2)
      do i = 1, ns
         j = size(smooths(i)%basis, 2)
         xorig(:, col_orig:col_orig + j - 1) = smooths(i)%basis
         xfp(:, col_orig:col_orig + j - 1) = matmul(smooths(i)%basis, reps(i)%transform)
         bmap(col_orig:col_orig + j - 1, col_orig:col_orig + j - 1) = reps(i)%transform
         if (reps(i)%nullity > 0) then
            theta(col_orig:col_orig + reps(i)%nullity - 1) = &
               beta(col_fit:col_fit + reps(i)%nullity - 1)
            col_fit = col_fit + reps(i)%nullity
         end if
         if (reps(i)%rank > 0) then
            theta(col_orig + reps(i)%nullity:col_orig + j - 1) = u(us:us + reps(i)%rank - 1)
            spcol(col_orig + reps(i)%nullity:col_orig + j - 1) = exp(-2.0_dp * eta(eta_s(i)))
            us = us + reps(i)%rank
         end if
         col_orig = col_orig + j
      end do
      result%coefficients = matmul(bmap, theta)
      allocate(result%sp(ns), result%smooth_sd(ns))
      result%sp = 0.0_dp
      result%smooth_sd = 0.0_dp
      do i = 1, ns
         if (eta_s(i) > 0) then
            result%sp(i) = exp(-2.0_dp * eta(eta_s(i)))
            result%smooth_sd(i) = sqrt(scale) * exp(eta(eta_s(i)))
         end if
      end do
      result%variance_parameters = eta
      result%random_effects = u
      result%scale = scale
      result%family = family
      result%reml = reml
      result%deviance = objective
      result%log_likelihood = -0.5_dp * objective
      if (family == gamm4_family_gaussian) then
         if (reml) then
            result%method = 'lmer.REML'
         else
            result%method = 'lmer.ML'
         end if
      else
         result%method = 'glmer.Laplace'
      end if

      allocate(gam_eta(n), cond_eta(n), mu(n))
      gam_eta = offset + matmul(xfit, beta)
      if (qs > 0) gam_eta = gam_eta + matmul(zs, u(1:qs))
      cond_eta = gam_eta
      qr = size(zr, 2)
      if (qr > 0) cond_eta = cond_eta + matmul(zr, u(qs + 1:qs + qr))
      result%linear_predictor = gam_eta
      result%conditional_linear_predictor = cond_eta
      call inverse_link(gam_eta, family, result%fitted)
      call inverse_link(cond_eta, family, result%conditional_fitted)
      result%residuals = y - result%conditional_fitted

      allocate(v(n))
      if (family == gamm4_family_gaussian) then
         v = scale / weights
      else
         v = max(scale / max(work_w, 1.0e-12_dp), 1.0e-12_dp)
      end if
      if (size(zr, 2) > 0) then
         phi = gre
      else
         allocate(phi(0, 0))
      end if
      call gamm4_get_vb(v, zr, phi, scale, xorig, xfp, spcol, bmap, vbres)
      if (vbres%status /= 0) then
         status = vbres%status
         return
      end if
      result%covariance = vbres%vb
      allocate(result%edf(porig))
      do i = 1, porig
         result%edf(i) = sum(result%covariance(i, :) * vbres%xvx(i, :))
      end do
   end subroutine assemble_result

   subroutine inverse_link(eta, family, mu)
      real(dp), intent(in) :: eta(:) !! Linear predictor values to transform to response means.
      integer, intent(in) :: family !! Gaussian or supported lme4 response-family code.
      real(dp), allocatable, intent(out) :: mu(:) !! Response-scale means corresponding to `eta`.
      integer :: i

      allocate(mu(size(eta)))
      select case (family)
      case (gamm4_family_gaussian)
         mu = eta
      case (family_binomial)
         do i = 1, size(eta)
            mu(i) = 1.0_dp / (1.0_dp + exp(-max(-35.0_dp, min(35.0_dp, eta(i)))))
         end do
      case default
         do i = 1, size(eta)
            mu(i) = exp(max(-35.0_dp, min(35.0_dp, eta(i))))
         end do
      end select
   end subroutine inverse_link

   subroutine initialize_eta(eta, eta_s, terms, ns)
      real(dp), intent(out) :: eta(:) !! Initial variance/covariance parameter vector.
      integer, intent(in) :: eta_s(:) !! Smooth-to-variance-parameter index map.
      type(random_term_t), intent(in) :: terms(:) !! Ordinary random terms whose covariance parameters follow smooth parameters.
      integer, intent(in) :: ns !! Number of smooth terms used to delimit the leading parameter block.
      integer :: i, j, p0, q, np

      eta = 0.0_dp
      do i = 1, ns
         if (eta_s(i) > 0) eta(eta_s(i)) = log(0.5_dp)
      end do
      p0 = count(eta_s > 0)
      do i = 1, size(terms)
         q = terms(i)%n_coefficients()
         np = terms(i)%n_parameters()
         select case (terms(i)%covariance_structure)
         case (covariance_unstructured)
            j = 0
            do while (j < np)
               j = j + 1
               eta(p0 + j) = log(0.5_dp)
               if (j < np) then
                  if (mod(j, max(1, q)) /= 0) eta(p0 + j) = 0.0_dp
               end if
            end do
            call initialize_unstructured_eta(eta(p0 + 1:p0 + np), q)
         case (covariance_diagonal)
            eta(p0 + 1:p0 + np) = log(0.5_dp)
         case (covariance_compound_symmetry, covariance_ar1)
            eta(p0 + 1) = log(0.5_dp)
            eta(p0 + 2) = 0.0_dp
         end select
         p0 = p0 + np
      end do
   end subroutine initialize_eta

   subroutine initialize_unstructured_eta(values, q)
      real(dp), intent(out) :: values(:) !! Packed lower-triangular covariance-factor parameters to initialize.
      integer, intent(in) :: q !! Number of random coefficients per grouping level.
      integer :: i, j, k

      values = 0.0_dp
      k = 0
      do j = 1, q
         do i = j, q
            k = k + 1
            if (i == j) values(k) = log(0.5_dp)
         end do
      end do
   end subroutine initialize_unstructured_eta

   subroutine eta_bounds(index, qs, terms, eta, lower, upper, control)
      integer, intent(in) :: index !! One-based variance-parameter index whose optimization interval is requested.
      integer, intent(in) :: qs !! Number of smooth random-effect columns; retained for API clarity in parameter partitioning.
      type(random_term_t), intent(in) :: terms(:) !! Ordinary random-effect terms defining covariance parameter positions.
      real(dp), intent(in) :: eta(:) !! Current complete parameter vector, used only for dimensional consistency.
      real(dp), intent(out) :: lower !! Lower coordinate-search bound.
      real(dp), intent(out) :: upper !! Upper coordinate-search bound.
      type(gamm4_control_t), intent(in) :: control !! Log-SD and off-diagonal search bounds.
      integer :: i, k, np, nsmooth, p0, q, a, b

      if (qs < 0 .or. size(eta) < index) then
         lower = control%lower_log_sd
         upper = control%upper_log_sd
         return
      end if
      nsmooth = size(eta) - random_parameter_count(terms)
      if (index <= nsmooth) then
         lower = control%lower_log_sd
         upper = control%upper_log_sd
         return
      end if
      p0 = nsmooth
      do i = 1, size(terms)
         np = terms(i)%n_parameters()
         if (index <= p0 + np) then
            select case (terms(i)%covariance_structure)
            case (covariance_unstructured)
               q = terms(i)%n_coefficients()
               k = 0
               do a = 1, q
                  do b = a, q
                     k = k + 1
                     if (p0 + k == index) then
                        if (a == b) then
                           lower = control%lower_log_sd
                           upper = control%upper_log_sd
                        else
                           lower = control%lower_offdiag
                           upper = control%upper_offdiag
                        end if
                        return
                     end if
                  end do
               end do
            case (covariance_diagonal)
               lower = control%lower_log_sd
               upper = control%upper_log_sd
               return
            case default
               if (index == p0 + 1) then
                  lower = control%lower_log_sd
                  upper = control%upper_log_sd
               else
                  lower = control%lower_offdiag
                  upper = control%upper_offdiag
               end if
               return
            end select
         end if
         p0 = p0 + np
      end do
      lower = control%lower_log_sd
      upper = control%upper_log_sd
   end subroutine eta_bounds

   integer function random_parameter_count(terms) result(total)
      type(random_term_t), intent(in) :: terms(:) !! Ordinary random-effect terms whose covariance parameter counts are summed.
      integer :: i

      total = 0
      do i = 1, size(terms)
         total = total + terms(i)%n_parameters()
      end do
   end function random_parameter_count



   logical pure function supported_family(family) result(ok)
      integer, intent(in) :: family !! Candidate response-family code.

      ok = family == gamm4_family_gaussian .or. family == family_binomial .or. family == family_poisson .or. &
         family == family_gamma .or. family == family_inverse_gaussian .or. family == family_negative_binomial
   end function supported_family

   subroutine initialize_result(result)
      type(gamm4_result_t), intent(out) :: result !! Result object reset to a deterministic not-yet-fitted state.

      result%status = -1
      result%converged = .false.
      result%message = 'not fitted'
   end subroutine initialize_result

   subroutine fail_result(result, status, code, message)
      type(gamm4_result_t), intent(inout) :: result !! Result object marked as failed with the supplied code and message.
      integer, intent(out) :: status !! Procedure status set to `code`.
      integer, intent(in) :: code !! Nonzero failure code identifying the validation or numerical stage.
      character(len=*), intent(in) :: message !! Human-readable failure description stored in the result.

      status = code
      result%status = code
      result%converged = .false.
      result%message = message
   end subroutine fail_result

end module gamm4_fit_mod
