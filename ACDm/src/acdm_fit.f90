! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use acdm_kinds, only : dp, huge_penalty, ACDM_SUCCESS, ACDM_BAD_INPUT, &
                         ACDM_BAD_PARAMETER, ACDM_NOT_CONVERGED, &
                         ACDM_SINGULAR
  use acdm_math, only : rng_state, seed_rng, random_uniform, invert_matrix
  use acdm_distributions, only : distribution_parameter_count, &
       DIST_EXPONENTIAL, DIST_WEIBULL, DIST_BURR, DIST_GENGAMMA, &
       DIST_GENF, DIST_QWEIBULL, DIST_MIXQWE, DIST_MIXQWW, &
       DIST_MIXINVGAUSS, DIST_BIRNBAUM_SAUNDERS
  use acdm_models, only : acd_order, model_parameter_count, &
       default_model_parameters, acd_loglik, simulate_acd, &
       MODEL_ACD, MODEL_LACD1, MODEL_LACD2, MODEL_EXACD, MODEL_AMACD, &
       MODEL_ABACD, MODEL_AACD, MODEL_TACD, MODEL_BACD, MODEL_BCACD, &
       MODEL_SNIACD, MODEL_LSNIACD, MODEL_TAMACD
  implicit none
  private

  type, public :: acd_fit_options
    integer :: model = MODEL_ACD
    integer :: dist = DIST_EXPONENTIAL
    type(acd_order) :: order
    logical :: force_mean = .true.
    integer :: max_iterations = 4000
    integer :: restarts = 2
    real(dp) :: tolerance = 1.0e-7_dp
    integer :: seed = 12345
    logical :: compute_hessian = .true.
    logical :: compute_robust_se = .true.
  end type acd_fit_options

  type, public :: acd_fit_result
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: model_parameters(:)
    real(dp), allocatable :: distribution_parameters(:)
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: hessian(:, :)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: robust_covariance(:, :)
    real(dp), allocatable :: robust_standard_errors(:)
    real(dp) :: loglik = -huge_penalty
    real(dp) :: aic = huge_penalty
    real(dp) :: bic = huge_penalty
    real(dp) :: mse = huge_penalty
    integer :: status = ACDM_BAD_INPUT
    integer :: convergence = 1
    integer :: evaluations = 0
    integer :: iterations = 0
  end type acd_fit_result

  type :: fit_context
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: exogenous(:, :)
    real(dp), allocatable :: breakpoints(:)
    integer, allocatable :: new_day(:)
    integer :: model
    integer :: dist
    type(acd_order) :: order
    logical :: force_mean
    integer :: n_model
    integer :: n_dist
    integer :: n_total
    logical, allocatable :: fixed_mask(:)
    real(dp), allocatable :: template(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    integer :: evaluations = 0
  end type fit_context

  public :: acd_fit_model, acd_score_matrix, acd_bootstrap_se
  public :: default_distribution_parameters, default_parameter_bounds
  public :: forecast_acd

contains

  function default_distribution_parameters(dist) result(par)
    integer, intent(in) :: dist
    real(dp), allocatable :: par(:)

    select case (dist)
    case (DIST_EXPONENTIAL)
      allocate(par(0))
    case (DIST_WEIBULL)
      par = [0.8_dp]
    case (DIST_BURR)
      par = [1.2_dp, 0.3_dp]
    case (DIST_GENGAMMA)
      par = [2.0_dp, 0.8_dp]
    case (DIST_GENF)
      par = [0.8_dp, 2.0_dp, 1.2_dp]
    case (DIST_QWEIBULL)
      par = [1.2_dp, 1.15_dp]
    case (DIST_MIXQWE)
      par = [0.7_dp, 1.3_dp, 1.15_dp, 0.8_dp]
    case (DIST_MIXQWW)
      par = [0.7_dp, 1.3_dp, 1.15_dp, 0.8_dp, 1.2_dp]
    case (DIST_MIXINVGAUSS)
      par = [0.7_dp, 0.2_dp, 0.4_dp]
    case (DIST_BIRNBAUM_SAUNDERS)
      par = [1.0_dp]
    case default
      allocate(par(0))
    end select
  end function default_distribution_parameters

  subroutine acd_fit_model(x, options, result, start, lower, upper, &
                           breakpoints, exogenous, new_day, fixed_mask)
    real(dp), intent(in) :: x(:)
    type(acd_fit_options), intent(in) :: options
    type(acd_fit_result), intent(out) :: result
    real(dp), intent(in), optional :: start(:), lower(:), upper(:)
    real(dp), intent(in), optional :: breakpoints(:)
    real(dp), intent(in), optional :: exogenous(:, :)
    integer, intent(in), optional :: new_day(:)
    logical, intent(in), optional :: fixed_mask(:)

    type(fit_context) :: ctx
    real(dp), allocatable :: start_free(:), lo_free(:), hi_free(:)
    real(dp), allocatable :: best_free(:), best_full(:), hfree(:, :), hinv(:, :)
    real(dp), allocatable :: score(:, :), bread(:, :), meat(:, :)
    real(dp), allocatable :: mstart(:), dstart(:), mu(:), residual(:)
    type(rng_state) :: rng
    real(dp) :: best_obj, obj, mean_x
    integer :: nbreak, n_exo, nfree, i, j, st, restart, it, ev

    result%status = ACDM_BAD_INPUT
    if (size(x) < 5 .or. any(x <= 0.0_dp)) return
    mean_x = sum(x) / real(size(x), dp)
    nbreak = 0
    if (present(breakpoints)) nbreak = size(breakpoints)
    n_exo = 0
    if (present(exogenous)) then
      if (size(exogenous, 1) /= size(x)) return
      n_exo = size(exogenous, 2)
    end if

    ctx%model = options%model
    ctx%dist = options%dist
    ctx%order = options%order
    ctx%force_mean = options%force_mean
    ctx%n_model = model_parameter_count(options%model, options%order, nbreak) + n_exo
    ctx%n_dist = distribution_parameter_count(options%dist)
    if (ctx%n_model < 1 .or. ctx%n_dist < 0) return
    ctx%n_total = ctx%n_model + ctx%n_dist
    ctx%x = x
    if (present(exogenous)) ctx%exogenous = exogenous
    if (present(breakpoints)) ctx%breakpoints = breakpoints
    if (present(new_day)) ctx%new_day = new_day

    allocate(ctx%template(ctx%n_total), ctx%fixed_mask(ctx%n_total))
    allocate(ctx%lower(ctx%n_total), ctx%upper(ctx%n_total))
    mstart = default_model_parameters(options%model, options%order, mean_x, nbreak)
    dstart = default_distribution_parameters(options%dist)
    ctx%template(1:ctx%n_model - n_exo) = mstart
    if (n_exo > 0) ctx%template(ctx%n_model - n_exo + 1:ctx%n_model) = 0.0_dp
    if (ctx%n_dist > 0) ctx%template(ctx%n_model + 1:ctx%n_total) = dstart
    if (present(start)) then
      if (size(start) /= ctx%n_total) return
      ctx%template = start
    end if
    ctx%fixed_mask = .false.
    if (present(fixed_mask)) then
      if (size(fixed_mask) /= ctx%n_total) return
      ctx%fixed_mask = fixed_mask
    end if

    call default_parameter_bounds(options%model, options%order, options%dist, &
                                  mean_x, nbreak, n_exo, ctx%lower, ctx%upper)
    if (present(lower)) then
      if (size(lower) /= ctx%n_total) return
      ctx%lower = lower
    end if
    if (present(upper)) then
      if (size(upper) /= ctx%n_total) return
      ctx%upper = upper
    end if
    if (any(ctx%lower >= ctx%upper)) return
    ctx%template = min(ctx%upper, max(ctx%lower, ctx%template))

    nfree = count(.not. ctx%fixed_mask)
    if (nfree == 0) then
      best_full = ctx%template
      best_obj = full_objective(ctx, best_full)
      it = 0
      ev = ctx%evaluations
    else
      allocate(start_free(nfree), lo_free(nfree), hi_free(nfree), best_free(nfree))
      call pack_free(ctx%template, ctx%fixed_mask, start_free)
      call pack_free(ctx%lower, ctx%fixed_mask, lo_free)
      call pack_free(ctx%upper, ctx%fixed_mask, hi_free)
      call seed_rng(rng, options%seed)
      best_obj = huge_penalty
      best_free = start_free
      it = 0
      do restart = 0, max(0, options%restarts)
        if (restart > 0) then
          do i = 1, nfree
            start_free(i) = best_free(i) + 0.10_dp * &
              (hi_free(i) - lo_free(i)) * (random_uniform(rng) - 0.5_dp)
          end do
          start_free = min(hi_free, max(lo_free, start_free))
        end if
        call nelder_mead(ctx, start_free, lo_free, hi_free, &
                         options%max_iterations, options%tolerance, &
                         best_full_free=best_full, best_value=obj, &
                         iterations=j, evaluations=ev)
        it = it + j
        if (obj < best_obj) then
          best_obj = obj
          best_free = best_full
        end if
      end do
      call expand_free(ctx%template, ctx%fixed_mask, best_free, best_full)
    end if

    allocate(mu(size(x)), residual(size(x)))
    obj = evaluate_full(ctx, best_full, mu, residual, st)
    if (st /= ACDM_SUCCESS .or. obj >= huge_penalty / 10.0_dp) then
      result%status = ACDM_NOT_CONVERGED
      result%convergence = 1
      return
    end if

    result%parameters = best_full
    result%model_parameters = best_full(1:ctx%n_model)
    if (ctx%n_dist > 0) then
      result%distribution_parameters = best_full(ctx%n_model + 1:ctx%n_total)
    else
      allocate(result%distribution_parameters(0))
    end if
    result%mu = mu
    result%residuals = residual
    result%loglik = -obj
    result%aic = 2.0_dp * real(nfree, dp) + 2.0_dp * obj
    result%bic = real(nfree, dp) * log(real(size(x), dp)) + 2.0_dp * obj
    result%mse = sum((x - mu)**2) / real(size(x), dp)
    result%evaluations = ctx%evaluations
    result%iterations = it
    result%status = ACDM_SUCCESS
    result%convergence = 0

    allocate(result%standard_errors(ctx%n_total))
    result%standard_errors = 0.0_dp
    if (options%compute_hessian .and. nfree > 0) then
      allocate(hfree(nfree, nfree), hinv(nfree, nfree))
      call numerical_hessian(ctx, best_free, lo_free, hi_free, hfree)
      call invert_matrix(hfree, hinv, st)
      allocate(result%hessian(ctx%n_total, ctx%n_total))
      allocate(result%covariance(ctx%n_total, ctx%n_total))
      result%hessian = 0.0_dp
      result%covariance = 0.0_dp
      call scatter_matrix(hfree, ctx%fixed_mask, result%hessian)
      if (st == ACDM_SUCCESS) then
        call scatter_matrix(hinv, ctx%fixed_mask, result%covariance)
        do i = 1, ctx%n_total
          if (result%covariance(i, i) > 0.0_dp) then
            result%standard_errors(i) = sqrt(result%covariance(i, i))
          end if
        end do
      end if

      if (options%compute_robust_se .and. options%dist == DIST_EXPONENTIAL &
          .and. st == ACDM_SUCCESS) then
        call acd_score_matrix_context(ctx, best_free, lo_free, hi_free, score, st)
        if (st == ACDM_SUCCESS) then
          meat = matmul(transpose(score), score)
          bread = matmul(hinv, matmul(meat, hinv))
          allocate(result%robust_covariance(ctx%n_total, ctx%n_total))
          allocate(result%robust_standard_errors(ctx%n_total))
          result%robust_covariance = 0.0_dp
          result%robust_standard_errors = 0.0_dp
          call scatter_matrix(bread, ctx%fixed_mask, result%robust_covariance)
          do i = 1, ctx%n_total
            if (result%robust_covariance(i, i) > 0.0_dp) then
              result%robust_standard_errors(i) = &
                sqrt(result%robust_covariance(i, i))
            end if
          end do
        end if
      end if
    end if
  end subroutine acd_fit_model

  subroutine default_parameter_bounds(model, order, dist, mean_x, nbreak, &
                                      n_exo, lower, upper)
    integer, intent(in) :: model, dist, nbreak, n_exo
    type(acd_order), intent(in) :: order
    real(dp), intent(in) :: mean_x
    real(dp), intent(out) :: lower(:), upper(:)
    integer :: nm, nd, pos, jn

    nm = model_parameter_count(model, order, nbreak)
    nd = distribution_parameter_count(dist)
    lower = -5.0_dp
    upper = 5.0_dp
    if (size(lower) /= nm + n_exo + nd) return

    select case (model)
    case (MODEL_ACD, MODEL_AMACD, MODEL_ABACD, MODEL_AACD, MODEL_BACD, &
          MODEL_SNIACD)
      lower(1) = 1.0e-8_dp
      upper(1) = max(10.0_dp, 10.0_dp * mean_x)
    case (MODEL_TACD, MODEL_TAMACD)
      jn = nbreak + 1
      lower(1:jn) = 1.0e-8_dp
      upper(1:jn) = max(10.0_dp, 10.0_dp * mean_x)
    case default
      lower(1) = -10.0_dp
      upper(1) = 10.0_dp
    end select

    select case (model)
    case (MODEL_ABACD, MODEL_AACD)
      pos = 2 + order%p + order%q
      lower(pos) = -0.999_dp
      upper(pos) = 0.999_dp
      lower(pos + 1) = -5.0_dp
      upper(pos + 1) = 5.0_dp
      lower(pos + 2:pos + 3) = 0.05_dp
      upper(pos + 2:pos + 3) = 5.0_dp
    case (MODEL_BACD)
      pos = 2 + order%p + order%q
      lower(pos:pos + 1) = 0.05_dp
      upper(pos:pos + 1) = 5.0_dp
    case (MODEL_BCACD)
      pos = 2 + order%p + order%q
      lower(pos) = 0.05_dp
      upper(pos) = 5.0_dp
    end select

    pos = nm + n_exo + 1
    select case (dist)
    case (DIST_WEIBULL)
      lower(pos) = 0.1_dp
      upper(pos) = 15.0_dp
    case (DIST_BURR)
      lower(pos:pos + 1) = [0.2_dp, 0.01_dp]
      upper(pos:pos + 1) = [20.0_dp, 10.0_dp]
    case (DIST_GENGAMMA)
      lower(pos:pos + 1) = [0.1_dp, 0.1_dp]
      upper(pos:pos + 1) = [20.0_dp, 10.0_dp]
    case (DIST_GENF)
      lower(pos:pos + 2) = 0.1_dp
      upper(pos:pos + 2) = 20.0_dp
    case (DIST_QWEIBULL)
      lower(pos:pos + 1) = [0.1_dp, 0.50_dp]
      upper(pos:pos + 1) = [15.0_dp, 1.95_dp]
    case (DIST_MIXQWE)
      lower(pos:pos + 3) = [0.01_dp, 0.1_dp, 0.50_dp, 0.01_dp]
      upper(pos:pos + 3) = [0.99_dp, 15.0_dp, 1.95_dp, 10.0_dp]
    case (DIST_MIXQWW)
      lower(pos:pos + 4) = [0.01_dp, 0.1_dp, 0.50_dp, 0.01_dp, 0.1_dp]
      upper(pos:pos + 4) = [0.99_dp, 15.0_dp, 1.95_dp, 10.0_dp, 10.0_dp]
    case (DIST_MIXINVGAUSS)
      lower(pos:pos + 2) = [0.01_dp, 0.01_dp, 0.0_dp]
      upper(pos:pos + 2) = [10.0_dp, 10.0_dp, 10.0_dp]
    case (DIST_BIRNBAUM_SAUNDERS)
      lower(pos) = 0.05_dp
      upper(pos) = 10.0_dp
    end select
  end subroutine default_parameter_bounds

  function full_objective(ctx, theta) result(value)
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: theta(:)
    real(dp) :: value
    real(dp), allocatable :: mu(:), residual(:)
    integer :: st

    allocate(mu(size(ctx%x)), residual(size(ctx%x)))
    value = evaluate_full(ctx, theta, mu, residual, st)
  end function full_objective

  function free_objective(ctx, free) result(value)
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: free(:)
    real(dp) :: value
    real(dp), allocatable :: full(:)

    call expand_free(ctx%template, ctx%fixed_mask, free, full)
    value = full_objective(ctx, full)
  end function free_objective

  function evaluate_full(ctx, theta, mu, residual, status) result(obj)
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: mu(:), residual(:)
    integer, intent(out) :: status
    real(dp) :: obj, ll
    real(dp), allocatable :: mpar(:), dpar(:)

    ctx%evaluations = ctx%evaluations + 1
    if (size(theta) /= ctx%n_total .or. any(theta < ctx%lower) .or. &
        any(theta > ctx%upper)) then
      status = ACDM_BAD_PARAMETER
      obj = huge_penalty
      return
    end if
    mpar = theta(1:ctx%n_model)
    if (ctx%n_dist > 0) then
      dpar = theta(ctx%n_model + 1:ctx%n_total)
    else
      allocate(dpar(0))
    end if

    if (allocated(ctx%breakpoints)) then
      if (allocated(ctx%exogenous) .and. allocated(ctx%new_day)) then
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, &
             breakpoints=ctx%breakpoints, exogenous=ctx%exogenous, &
             new_day=ctx%new_day)
      else if (allocated(ctx%exogenous)) then
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, &
             breakpoints=ctx%breakpoints, exogenous=ctx%exogenous)
      else if (allocated(ctx%new_day)) then
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, &
             breakpoints=ctx%breakpoints, new_day=ctx%new_day)
      else
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, &
             breakpoints=ctx%breakpoints)
      end if
    else
      if (allocated(ctx%exogenous) .and. allocated(ctx%new_day)) then
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, &
             exogenous=ctx%exogenous, new_day=ctx%new_day)
      else if (allocated(ctx%exogenous)) then
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, exogenous=ctx%exogenous)
      else if (allocated(ctx%new_day)) then
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status, new_day=ctx%new_day)
      else
        ll = acd_loglik(ctx%x, ctx%model, ctx%order, mpar, ctx%dist, dpar, &
             ctx%force_mean, mu, residual, status)
      end if
    end if
    if (status /= ACDM_SUCCESS .or. .not. ieee_is_finite(ll)) then
      obj = huge_penalty
    else
      obj = -ll
    end if
  end function evaluate_full

  subroutine nelder_mead(ctx, start, lower, upper, maxit, tol, &
                         best_full_free, best_value, iterations, evaluations)
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: start(:), lower(:), upper(:)
    integer, intent(in) :: maxit
    real(dp), intent(in) :: tol
    real(dp), allocatable, intent(out) :: best_full_free(:)
    real(dp), intent(out) :: best_value
    integer, intent(out) :: iterations, evaluations

    real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: fr, fe, fc, spread, diam, alpha, gamma, rho, sigma, step
    integer :: n, i, j, ilo, ihi, inhi, ev0

    n = size(start)
    allocate(simplex(n, n + 1), f(n + 1), centroid(n), xr(n), xe(n), xc(n))
    alpha = 1.0_dp
    gamma = 2.0_dp
    rho = 0.5_dp
    sigma = 0.5_dp
    simplex(:, 1) = start
    do j = 1, n
      simplex(:, j + 1) = start
      step = max(0.05_dp * (upper(j) - lower(j)), &
                 0.05_dp * max(1.0_dp, abs(start(j))))
      simplex(j, j + 1) = min(upper(j), max(lower(j), start(j) + step))
      if (abs(simplex(j, j + 1) - start(j)) <= epsilon(1.0_dp)) then
        simplex(j, j + 1) = max(lower(j), start(j) - step)
      end if
    end do
    ev0 = ctx%evaluations
    do j = 1, n + 1
      f(j) = free_objective(ctx, simplex(:, j))
    end do

    iterations = 0
    do i = 1, maxit
      iterations = i
      call simplex_order(f, ilo, ihi, inhi)
      spread = maxval(abs(f - f(ilo))) / max(1.0_dp, abs(f(ilo)))
      diam = 0.0_dp
      do j = 1, n + 1
        diam = max(diam, maxval(abs(simplex(:, j) - simplex(:, ilo)) / &
                   max(1.0_dp, abs(simplex(:, ilo)))))
      end do
      if (spread < tol .and. diam < sqrt(tol)) exit

      centroid = (sum(simplex, dim=2) - simplex(:, ihi)) / real(n, dp)
      xr = centroid + alpha * (centroid - simplex(:, ihi))
      xr = min(upper, max(lower, xr))
      fr = free_objective(ctx, xr)

      if (fr < f(ilo)) then
        xe = centroid + gamma * (xr - centroid)
        xe = min(upper, max(lower, xe))
        fe = free_objective(ctx, xe)
        if (fe < fr) then
          simplex(:, ihi) = xe
          f(ihi) = fe
        else
          simplex(:, ihi) = xr
          f(ihi) = fr
        end if
      else if (fr < f(inhi)) then
        simplex(:, ihi) = xr
        f(ihi) = fr
      else
        if (fr < f(ihi)) then
          xc = centroid + rho * (xr - centroid)
        else
          xc = centroid - rho * (centroid - simplex(:, ihi))
        end if
        xc = min(upper, max(lower, xc))
        fc = free_objective(ctx, xc)
        if (fc < min(fr, f(ihi))) then
          simplex(:, ihi) = xc
          f(ihi) = fc
        else
          do j = 1, n + 1
            if (j /= ilo) then
              simplex(:, j) = simplex(:, ilo) + &
                               sigma * (simplex(:, j) - simplex(:, ilo))
              simplex(:, j) = min(upper, max(lower, simplex(:, j)))
              f(j) = free_objective(ctx, simplex(:, j))
            end if
          end do
        end if
      end if
    end do
    call simplex_order(f, ilo, ihi, inhi)
    best_full_free = simplex(:, ilo)
    best_value = f(ilo)
    evaluations = ctx%evaluations - ev0
  end subroutine nelder_mead

  pure subroutine simplex_order(f, ilo, ihi, inhi)
    real(dp), intent(in) :: f(:)
    integer, intent(out) :: ilo, ihi, inhi
    integer :: i

    ilo = minloc(f, dim=1)
    ihi = maxloc(f, dim=1)
    inhi = ilo
    do i = 1, size(f)
      if (i /= ihi) then
        if (inhi == ilo .or. f(i) > f(inhi)) inhi = i
      end if
    end do
  end subroutine simplex_order

  subroutine numerical_hessian(ctx, x, lower, upper, hess)
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: x(:), lower(:), upper(:)
    real(dp), intent(out) :: hess(:, :)
    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
    real(dp) :: f0, fp, fm, fpp, fpm, fmp, fmm, hi, hj
    integer :: n, i, j

    n = size(x)
    allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
    f0 = free_objective(ctx, x)
    hess = 0.0_dp
    do i = 1, n
      hi = min(1.0e-4_dp * max(1.0_dp, abs(x(i))), &
               0.25_dp * (upper(i) - lower(i)))
      hi = max(hi, 1.0e-6_dp)
      xp = x
      xm = x
      xp(i) = min(upper(i), x(i) + hi)
      xm(i) = max(lower(i), x(i) - hi)
      hi = 0.5_dp * (xp(i) - xm(i))
      fp = free_objective(ctx, xp)
      fm = free_objective(ctx, xm)
      hess(i, i) = (fp - 2.0_dp * f0 + fm) / (hi * hi)
      do j = i + 1, n
        hj = min(1.0e-4_dp * max(1.0_dp, abs(x(j))), &
                 0.25_dp * (upper(j) - lower(j)))
        hj = max(hj, 1.0e-6_dp)
        xpp = x
        xpm = x
        xmp = x
        xmm = x
        xpp(i) = min(upper(i), x(i) + hi)
        xpp(j) = min(upper(j), x(j) + hj)
        xpm(i) = min(upper(i), x(i) + hi)
        xpm(j) = max(lower(j), x(j) - hj)
        xmp(i) = max(lower(i), x(i) - hi)
        xmp(j) = min(upper(j), x(j) + hj)
        xmm(i) = max(lower(i), x(i) - hi)
        xmm(j) = max(lower(j), x(j) - hj)
        hi = 0.5_dp * (xpp(i) - xmp(i))
        hj = 0.5_dp * (xpp(j) - xpm(j))
        fpp = free_objective(ctx, xpp)
        fpm = free_objective(ctx, xpm)
        fmp = free_objective(ctx, xmp)
        fmm = free_objective(ctx, xmm)
        hess(i, j) = (fpp - fpm - fmp + fmm) / (4.0_dp * hi * hj)
        hess(j, i) = hess(i, j)
      end do
    end do
  end subroutine numerical_hessian

  subroutine acd_score_matrix(x, options, parameters, score, status, &
                              breakpoints, exogenous, new_day)
    real(dp), intent(in) :: x(:), parameters(:)
    type(acd_fit_options), intent(in) :: options
    real(dp), allocatable, intent(out) :: score(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: breakpoints(:)
    real(dp), intent(in), optional :: exogenous(:, :)
    integer, intent(in), optional :: new_day(:)
    type(fit_context) :: ctx
    real(dp), allocatable :: lower(:), upper(:)
    integer :: nb, ne

    nb = 0
    if (present(breakpoints)) nb = size(breakpoints)
    ne = 0
    if (present(exogenous)) ne = size(exogenous, 2)
    ctx%model = options%model
    ctx%dist = options%dist
    ctx%order = options%order
    ctx%force_mean = options%force_mean
    ctx%n_model = model_parameter_count(options%model, options%order, nb) + ne
    ctx%n_dist = distribution_parameter_count(options%dist)
    ctx%n_total = ctx%n_model + ctx%n_dist
    if (size(parameters) /= ctx%n_total) then
      status = ACDM_BAD_INPUT
      return
    end if
    ctx%x = x
    if (present(breakpoints)) ctx%breakpoints = breakpoints
    if (present(exogenous)) ctx%exogenous = exogenous
    if (present(new_day)) ctx%new_day = new_day
    allocate(ctx%template(ctx%n_total), ctx%fixed_mask(ctx%n_total))
    allocate(ctx%lower(ctx%n_total), ctx%upper(ctx%n_total))
    ctx%template = parameters
    ctx%fixed_mask = .false.
    call default_parameter_bounds(options%model, options%order, options%dist, &
         sum(x) / real(size(x), dp), nb, ne, ctx%lower, ctx%upper)
    lower = ctx%lower
    upper = ctx%upper
    call acd_score_matrix_context(ctx, parameters, lower, upper, score, status)
  end subroutine acd_score_matrix

  subroutine acd_score_matrix_context(ctx, theta_free, lower, upper, &
                                      score, status)
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: theta_free(:), lower(:), upper(:)
    real(dp), allocatable, intent(out) :: score(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: tp(:), tm(:), fp(:), fm(:)
    real(dp) :: h
    integer :: i, n, p

    n = size(ctx%x)
    p = size(theta_free)
    allocate(score(n, p), tp(p), tm(p), fp(n), fm(n))
    status = ACDM_SUCCESS
    do i = 1, p
      h = min(1.0e-5_dp * max(1.0_dp, abs(theta_free(i))), &
              0.2_dp * (upper(i) - lower(i)))
      h = max(h, 1.0e-7_dp)
      tp = theta_free
      tm = theta_free
      tp(i) = min(upper(i), tp(i) + h)
      tm(i) = max(lower(i), tm(i) - h)
      h = 0.5_dp * (tp(i) - tm(i))
      call loglik_contributions_free(ctx, tp, fp, status)
      if (status /= ACDM_SUCCESS) return
      call loglik_contributions_free(ctx, tm, fm, status)
      if (status /= ACDM_SUCCESS) return
      score(:, i) = (fp - fm) / (2.0_dp * h)
    end do
  end subroutine acd_score_matrix_context

  subroutine loglik_contributions_free(ctx, free, contrib, status)
    use acdm_distributions, only : distribution_logpdf
    type(fit_context), intent(inout) :: ctx
    real(dp), intent(in) :: free(:)
    real(dp), intent(out) :: contrib(:)
    integer, intent(out) :: status
    real(dp), allocatable :: full(:), mu(:), residual(:), dpar(:)
    real(dp) :: obj
    integer :: i

    call expand_free(ctx%template, ctx%fixed_mask, free, full)
    allocate(mu(size(ctx%x)), residual(size(ctx%x)))
    obj = evaluate_full(ctx, full, mu, residual, status)
    if (status /= ACDM_SUCCESS) return
    if (ctx%n_dist > 0) then
      dpar = full(ctx%n_model + 1:ctx%n_total)
    else
      allocate(dpar(0))
    end if
    do i = 1, size(ctx%x)
      contrib(i) = distribution_logpdf(residual(i), ctx%dist, dpar, &
                                       ctx%force_mean) - log(mu(i))
    end do
  end subroutine loglik_contributions_free

  subroutine pack_free(full, mask, free)
    real(dp), intent(in) :: full(:)
    logical, intent(in) :: mask(:)
    real(dp), intent(out) :: free(:)
    integer :: i, j

    j = 0
    do i = 1, size(full)
      if (.not. mask(i)) then
        j = j + 1
        free(j) = full(i)
      end if
    end do
  end subroutine pack_free

  subroutine expand_free(template, mask, free, full)
    real(dp), intent(in) :: template(:), free(:)
    logical, intent(in) :: mask(:)
    real(dp), allocatable, intent(out) :: full(:)
    integer :: i, j

    allocate(full(size(template)))
    full = template
    j = 0
    do i = 1, size(template)
      if (.not. mask(i)) then
        j = j + 1
        full(i) = free(j)
      end if
    end do
  end subroutine expand_free

  subroutine scatter_matrix(small, mask, full)
    real(dp), intent(in) :: small(:, :)
    logical, intent(in) :: mask(:)
    real(dp), intent(inout) :: full(:, :)
    integer :: i, j, ii, jj

    ii = 0
    do i = 1, size(mask)
      if (.not. mask(i)) then
        ii = ii + 1
        jj = 0
        do j = 1, size(mask)
          if (.not. mask(j)) then
            jj = jj + 1
            full(i, j) = small(ii, jj)
          end if
        end do
      end if
    end do
  end subroutine scatter_matrix

  subroutine forecast_acd(fit, options, n_ahead, forecast, status, &
                          breakpoints)
    type(acd_fit_result), intent(in) :: fit
    type(acd_fit_options), intent(in) :: options
    integer, intent(in) :: n_ahead
    real(dp), intent(out) :: forecast(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: breakpoints(:)
    type(rng_state) :: rng
    real(dp), allocatable :: errors(:), start_x(:), start_mu(:)
    integer :: k, nt

    if (size(forecast) /= n_ahead .or. .not. allocated(fit%parameters)) then
      status = ACDM_BAD_INPUT
      return
    end if
    k = max(options%order%p, options%order%q, options%order%r)
    if (k < 1) k = 1
    if (size(fit%mu) < k) then
      status = ACDM_BAD_INPUT
      return
    end if
    start_x = fit%residuals(size(fit%residuals) - k + 1:size(fit%residuals)) * &
              fit%mu(size(fit%mu) - k + 1:size(fit%mu))
    start_mu = fit%mu(size(fit%mu) - k + 1:size(fit%mu))
    nt = n_ahead + k
    allocate(errors(nt))
    errors = 1.0_dp
    call seed_rng(rng, options%seed)
    if (present(breakpoints)) then
      call simulate_acd(n_ahead, options%model, options%order, &
           fit%model_parameters, options%dist, fit%distribution_parameters, &
           forecast, status, rng, burn=k, errors=errors, &
           force_mean=options%force_mean, start_x=start_x, &
           start_mu=start_mu, breakpoints=breakpoints)
    else
      call simulate_acd(n_ahead, options%model, options%order, &
           fit%model_parameters, options%dist, fit%distribution_parameters, &
           forecast, status, rng, burn=k, errors=errors, &
           force_mean=options%force_mean, start_x=start_x, start_mu=start_mu)
    end if
  end subroutine forecast_acd

  subroutine acd_bootstrap_se(x, options, fit, nboot, bootstrap_se, status, &
                              breakpoints)
    real(dp), intent(in) :: x(:)
    type(acd_fit_options), intent(in) :: options
    type(acd_fit_result), intent(in) :: fit
    integer, intent(in) :: nboot
    real(dp), intent(out) :: bootstrap_se(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: breakpoints(:)
    type(acd_fit_result) :: bfit
    type(acd_fit_options) :: bopt
    type(rng_state) :: rng
    real(dp), allocatable :: sim(:), errors(:), pars(:, :), meanp(:)
    integer :: b, i, st, good

    if (nboot < 2 .or. size(bootstrap_se) /= size(fit%parameters)) then
      status = ACDM_BAD_INPUT
      return
    end if
    allocate(sim(size(x)), errors(size(x) + 50))
    allocate(pars(nboot, size(fit%parameters)))
    call seed_rng(rng, options%seed + 7919)
    bopt = options
    bopt%restarts = 0
    bopt%compute_hessian = .false.
    bopt%compute_robust_se = .false.
    good = 0
    do b = 1, nboot
      do i = 1, size(errors)
        errors(i) = fit%residuals(1 + int(random_uniform(rng) * &
                    real(size(fit%residuals), dp)))
      end do
      if (present(breakpoints)) then
        call simulate_acd(size(x), options%model, options%order, &
             fit%model_parameters, options%dist, fit%distribution_parameters, &
             sim, st, rng, burn=50, errors=errors, &
             force_mean=options%force_mean, breakpoints=breakpoints)
      else
        call simulate_acd(size(x), options%model, options%order, &
             fit%model_parameters, options%dist, fit%distribution_parameters, &
             sim, st, rng, burn=50, errors=errors, &
             force_mean=options%force_mean)
      end if
      if (st /= ACDM_SUCCESS) cycle
      if (present(breakpoints)) then
        call acd_fit_model(sim, bopt, bfit, start=fit%parameters, &
                           breakpoints=breakpoints)
      else
        call acd_fit_model(sim, bopt, bfit, start=fit%parameters)
      end if
      if (bfit%status == ACDM_SUCCESS) then
        good = good + 1
        pars(good, :) = bfit%parameters
      end if
    end do
    if (good < 2) then
      status = ACDM_NOT_CONVERGED
      return
    end if
    meanp = sum(pars(1:good, :), dim=1) / real(good, dp)
    do i = 1, size(bootstrap_se)
      bootstrap_se(i) = sqrt(sum((pars(1:good, i) - meanp(i))**2) / &
                             real(good - 1, dp))
    end do
    status = ACDM_SUCCESS
  end subroutine acd_bootstrap_se

end module acdm_fit
