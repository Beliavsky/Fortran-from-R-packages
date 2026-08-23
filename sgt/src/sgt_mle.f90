! SPDX-License-Identifier: GPL-3.0-or-later
module sgt_mle_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use sgt_kinds, only : dp
  use sgt_distribution, only : sgt_logpdf
  use sgt_special, only : invert_matrix, normal_cdf
  implicit none
  private
  integer, parameter :: name_len = 24

  type, public :: sgt_mle_result
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: varcov(:,:)
    real(dp), allocatable :: std_error(:)
    real(dp), allocatable :: z_score(:)
    real(dp), allocatable :: p_value(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: convcode = 10
    integer :: niter = 0
    character(len=name_len) :: best_method = ''
  end type sgt_mle_result

  type, public :: sgt_params
    real(dp) :: mu = 0.0_dp
    real(dp) :: sigma = 1.0_dp
    real(dp) :: lambda = 0.0_dp
    real(dp) :: p = 2.0_dp
    real(dp) :: q = 10.0_dp
  end type sgt_params

  abstract interface
    subroutine sgt_observation_model(theta, index, x, mu, sigma, lambda, p, q, status)
      import dp
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: index
      real(dp), intent(out) :: x, mu, sigma, lambda, p, q
      integer, intent(out) :: status
    end subroutine sgt_observation_model
  end interface

  public :: sgt_mle_model, sgt_mle_constant, sgt_observation_model

  real(dp), allocatable, save :: constant_data_context(:)
  type(sgt_params), save :: constant_start_context
  logical, save :: constant_free_context(5) = .false.
contains
  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  real(dp) function neg_loglik(theta, nobs, model, mean_cent, var_adj) result(value)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: nobs
    procedure(sgt_observation_model) :: model
    logical, intent(in) :: mean_cent, var_adj
    integer :: i, status
    real(dp) :: x, mu, sigma, lambda, p, q, lf
    value = 0.0_dp
    do i = 1, nobs
      call model(theta, i, x, mu, sigma, lambda, p, q, status)
      if (status /= 0) then
        value = huge(1.0_dp) / 16.0_dp
        return
      end if
      lf = sgt_logpdf(x, mu, sigma, lambda, p, q, mean_cent, var_adj)
      if (.not. ieee_is_finite(lf)) then
        value = huge(1.0_dp) / 16.0_dp
        return
      end if
      value = value - lf
      if (.not. ieee_is_finite(value)) then
        value = huge(1.0_dp) / 16.0_dp
        return
      end if
    end do
  end function neg_loglik

  subroutine numeric_gradient(theta, nobs, model, mean_cent, var_adj, grad)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: nobs
    procedure(sgt_observation_model) :: model
    logical, intent(in) :: mean_cent, var_adj
    real(dp), intent(out) :: grad(:)
    real(dp) :: xp(size(theta)), xm(size(theta)), hp, hm, h
    integer :: j
    do j = 1, size(theta)
      h = 1.0e-5_dp * max(1.0_dp, abs(theta(j)))
      xp = theta
      xm = theta
      xp(j) = xp(j) + h
      xm(j) = xm(j) - h
      hp = neg_loglik(xp, nobs, model, mean_cent, var_adj)
      hm = neg_loglik(xm, nobs, model, mean_cent, var_adj)
      grad(j) = (hp - hm) / (2.0_dp * h)
    end do
  end subroutine numeric_gradient

  subroutine numeric_hessian_loglik(theta, nobs, model, mean_cent, var_adj, hess)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: nobs
    procedure(sgt_observation_model) :: model
    logical, intent(in) :: mean_cent, var_adj
    real(dp), intent(out) :: hess(:,:)
    real(dp) :: xpp(size(theta)), xpm(size(theta)), xmp(size(theta)), xmm(size(theta))
    real(dp) :: xp(size(theta)), xm(size(theta)), f0, fp, fm, hi, hj
    integer :: i, j
    f0 = -neg_loglik(theta, nobs, model, mean_cent, var_adj)
    do i = 1, size(theta)
      hi = 1.0e-4_dp * max(1.0_dp, abs(theta(i)))
      xp = theta
      xm = theta
      xp(i) = xp(i) + hi
      xm(i) = xm(i) - hi
      fp = -neg_loglik(xp, nobs, model, mean_cent, var_adj)
      fm = -neg_loglik(xm, nobs, model, mean_cent, var_adj)
      hess(i,i) = (fp - 2.0_dp * f0 + fm) / (hi * hi)
      do j = 1, i - 1
        hj = 1.0e-4_dp * max(1.0_dp, abs(theta(j)))
        xpp = theta
        xpm = theta
        xmp = theta
        xmm = theta
        xpp(i) = xpp(i) + hi
        xpp(j) = xpp(j) + hj
        xpm(i) = xpm(i) + hi
        xpm(j) = xpm(j) - hj
        xmp(i) = xmp(i) - hi
        xmp(j) = xmp(j) + hj
        xmm(i) = xmm(i) - hi
        xmm(j) = xmm(j) - hj
        hess(i,j) = (-neg_loglik(xpp, nobs, model, mean_cent, var_adj) + &
          neg_loglik(xpm, nobs, model, mean_cent, var_adj) + &
          neg_loglik(xmp, nobs, model, mean_cent, var_adj) - &
          neg_loglik(xmm, nobs, model, mean_cent, var_adj)) / (4.0_dp * hi * hj)
        hess(j,i) = hess(i,j)
      end do
    end do
  end subroutine numeric_hessian_loglik

  subroutine nelder_mead(theta0, nobs, model, mean_cent, var_adj, max_iter, theta, value, iter, code)
    real(dp), intent(in) :: theta0(:)
    integer, intent(in) :: nobs, max_iter
    procedure(sgt_observation_model) :: model
    logical, intent(in) :: mean_cent, var_adj
    real(dp), intent(out) :: theta(:), value
    integer, intent(out) :: iter, code
    integer :: n, i, j, best, worst, second_worst
    real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: scale, fr, fe, fc, spread
    n = size(theta0)
    allocate(simplex(n,n+1), f(n+1), centroid(n), xr(n), xe(n), xc(n))
    simplex(:,1) = theta0
    do j = 1, n
      simplex(:,j+1) = theta0
      scale = 0.05_dp * max(1.0_dp, abs(theta0(j)))
      simplex(j,j+1) = simplex(j,j+1) + scale
    end do
    do i = 1, n + 1
      f(i) = neg_loglik(simplex(:,i), nobs, model, mean_cent, var_adj)
    end do
    code = 1
    do iter = 1, max_iter
      best = minloc(f, dim=1)
      worst = maxloc(f, dim=1)
      second_worst = best
      do i = 1, n + 1
        if (i == worst) cycle
        if (second_worst == worst .or. f(i) > f(second_worst)) second_worst = i
      end do
      spread = maxval(abs(f - f(best)))
      if (spread <= 1.0e-9_dp * max(1.0_dp, abs(f(best)))) then
        code = 0
        exit
      end if
      centroid = 0.0_dp
      do i = 1, n + 1
        if (i /= worst) centroid = centroid + simplex(:,i)
      end do
      centroid = centroid / real(n, dp)
      xr = centroid + (centroid - simplex(:,worst))
      fr = neg_loglik(xr, nobs, model, mean_cent, var_adj)
      if (fr < f(best)) then
        xe = centroid + 2.0_dp * (xr - centroid)
        fe = neg_loglik(xe, nobs, model, mean_cent, var_adj)
        if (fe < fr) then
          simplex(:,worst) = xe
          f(worst) = fe
        else
          simplex(:,worst) = xr
          f(worst) = fr
        end if
      else if (fr < f(second_worst)) then
        simplex(:,worst) = xr
        f(worst) = fr
      else
        if (fr < f(worst)) then
          xc = centroid + 0.5_dp * (xr - centroid)
        else
          xc = centroid + 0.5_dp * (simplex(:,worst) - centroid)
        end if
        fc = neg_loglik(xc, nobs, model, mean_cent, var_adj)
        if (fc < min(fr, f(worst))) then
          simplex(:,worst) = xc
          f(worst) = fc
        else
          do i = 1, n + 1
            if (i == best) cycle
            simplex(:,i) = simplex(:,best) + 0.5_dp * (simplex(:,i) - simplex(:,best))
            f(i) = neg_loglik(simplex(:,i), nobs, model, mean_cent, var_adj)
          end do
        end if
      end if
    end do
    best = minloc(f, dim=1)
    theta = simplex(:,best)
    value = f(best)
  end subroutine nelder_mead

  subroutine bfgs(theta0, nobs, model, mean_cent, var_adj, max_iter, theta, value, iter, code)
    real(dp), intent(in) :: theta0(:)
    integer, intent(in) :: nobs, max_iter
    procedure(sgt_observation_model) :: model
    logical, intent(in) :: mean_cent, var_adj
    real(dp), intent(out) :: theta(:), value
    integer, intent(out) :: iter, code
    integer :: n, i, ls
    real(dp), allocatable :: hmat(:,:), g(:), gnew(:), direction(:), trial(:), s(:), y(:), ident(:,:)
    real(dp) :: f, fnew, alpha, ys, gnorm
    n = size(theta0)
    allocate(hmat(n,n), g(n), gnew(n), direction(n), trial(n), s(n), y(n), ident(n,n))
    hmat = 0.0_dp
    ident = 0.0_dp
    do i = 1, n
      hmat(i,i) = 1.0_dp
      ident(i,i) = 1.0_dp
    end do
    theta = theta0
    f = neg_loglik(theta, nobs, model, mean_cent, var_adj)
    call numeric_gradient(theta, nobs, model, mean_cent, var_adj, g)
    code = 1
    do iter = 1, max_iter
      gnorm = maxval(abs(g))
      if (gnorm <= 1.0e-6_dp) then
        code = 0
        exit
      end if
      direction = -matmul(hmat, g)
      if (dot_product(direction, g) >= 0.0_dp) direction = -g
      alpha = 1.0_dp
      do ls = 1, 40
        trial = theta + alpha * direction
        fnew = neg_loglik(trial, nobs, model, mean_cent, var_adj)
        if (fnew <= f + 1.0e-4_dp * alpha * dot_product(g, direction)) exit
        alpha = 0.5_dp * alpha
      end do
      if (alpha < 1.0e-12_dp) exit
      s = trial - theta
      theta = trial
      fnew = neg_loglik(theta, nobs, model, mean_cent, var_adj)
      call numeric_gradient(theta, nobs, model, mean_cent, var_adj, gnew)
      y = gnew - g
      ys = dot_product(y, s)
      if (ys > 1.0e-12_dp * sqrt(max(dot_product(y,y) * dot_product(s,s), tiny(1.0_dp)))) then
        hmat = matmul(ident - outer_product(s,y) / ys, &
          matmul(hmat, ident - outer_product(y,s) / ys)) + outer_product(s,s) / ys
      else
        hmat = ident
      end if
      g = gnew
      f = fnew
      if (maxval(abs(s)) <= 1.0e-8_dp * max(1.0_dp, maxval(abs(theta)))) then
        code = 0
        exit
      end if
    end do
    value = neg_loglik(theta, nobs, model, mean_cent, var_adj)
  contains
    pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a), size(b))
      integer :: ii, jj
      do jj = 1, size(b)
        do ii = 1, size(a)
          c(ii,jj) = a(ii) * b(jj)
        end do
      end do
    end function outer_product
  end subroutine bfgs

  subroutine fill_inference(theta, nobs, model, mean_cent, var_adj, result)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: nobs
    procedure(sgt_observation_model) :: model
    logical, intent(in) :: mean_cent, var_adj
    type(sgt_mle_result), intent(inout) :: result
    real(dp), allocatable :: invh(:,:)
    integer :: i, status
    allocate(result%gradient(size(theta)), result%hessian(size(theta),size(theta)))
    allocate(result%varcov(size(theta),size(theta)), result%std_error(size(theta)))
    allocate(result%z_score(size(theta)), result%p_value(size(theta)))
    call numeric_gradient(theta, nobs, model, mean_cent, var_adj, result%gradient)
    result%gradient = -result%gradient
    call numeric_hessian_loglik(theta, nobs, model, mean_cent, var_adj, result%hessian)
    call invert_matrix(-result%hessian, invh, status)
    if (status == 0) then
      result%varcov = invh
      do i = 1, size(theta)
        if (result%varcov(i,i) > 0.0_dp) then
          result%std_error(i) = sqrt(result%varcov(i,i))
          result%z_score(i) = theta(i) / result%std_error(i)
          result%p_value(i) = 2.0_dp * normal_cdf(-abs(result%z_score(i)))
        else
          result%std_error(i) = nan_dp()
          result%z_score(i) = nan_dp()
          result%p_value(i) = nan_dp()
        end if
      end do
    else
      result%varcov = nan_dp()
      result%std_error = nan_dp()
      result%z_score = nan_dp()
      result%p_value = nan_dp()
    end if
  end subroutine fill_inference

  subroutine sgt_mle_model(nobs, theta_start, model, result, methods, mean_cent, var_adj, max_iter)
    integer, intent(in) :: nobs
    real(dp), intent(in) :: theta_start(:)
    procedure(sgt_observation_model) :: model
    type(sgt_mle_result), intent(out) :: result
    character(len=*), intent(in), optional :: methods(:)
    logical, intent(in), optional :: mean_cent, var_adj
    integer, intent(in), optional :: max_iter
    character(len=name_len), allocatable :: meth(:)
    real(dp), allocatable :: theta(:), best_theta(:)
    real(dp) :: value, best_value
    logical :: mc, va
    integer :: mi, iter, code, mit
    mc = .true.
    va = .true.
    mit = 1000
    if (present(mean_cent)) mc = mean_cent
    if (present(var_adj)) va = var_adj
    if (present(max_iter)) mit = max_iter
    if (present(methods)) then
      allocate(meth(size(methods)))
      meth = methods
    else
      allocate(meth(2))
      meth = [character(len=name_len) :: 'Nelder-Mead', 'BFGS']
    end if
    allocate(theta(size(theta_start)), best_theta(size(theta_start)))
    best_value = huge(1.0_dp)
    result%convcode = 10
    do mi = 1, size(meth)
      select case (trim(adjustl(meth(mi))))
      case ('Nelder-Mead', 'nelder-mead', 'NELDER-MEAD')
        call nelder_mead(theta_start, nobs, model, mc, va, mit, theta, value, iter, code)
      case ('BFGS', 'bfgs')
        call bfgs(theta_start, nobs, model, mc, va, mit, theta, value, iter, code)
      case default
        cycle
      end select
      if (ieee_is_finite(value) .and. value < best_value) then
        best_value = value
        best_theta = theta
        result%convcode = code
        result%niter = iter
        result%best_method = trim(meth(mi))
      end if
    end do
    allocate(result%estimate(size(theta_start)))
    result%estimate = best_theta
    result%loglik = -best_value
    call fill_inference(best_theta, nobs, model, mc, va, result)
  end subroutine sgt_mle_model

  subroutine constant_context_model(theta, index, x, mu, sigma, lambda, p, q, status)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: index
    real(dp), intent(out) :: x, mu, sigma, lambda, p, q
    integer, intent(out) :: status
    integer :: jj, kk
    mu = constant_start_context%mu
    sigma = constant_start_context%sigma
    lambda = constant_start_context%lambda
    p = constant_start_context%p
    q = constant_start_context%q
    kk = 0
    do jj = 1, 5
      if (.not. constant_free_context(jj)) cycle
      kk = kk + 1
      select case (jj)
      case (1)
        mu = theta(kk)
      case (2)
        sigma = exp(theta(kk))
      case (3)
        lambda = tanh(theta(kk))
      case (4)
        p = exp(theta(kk))
      case (5)
        q = exp(theta(kk))
      end select
    end do
    x = constant_data_context(index)
    status = 0
  end subroutine constant_context_model

  subroutine sgt_mle_constant(data, start, result, free, methods, mean_cent, var_adj, max_iter)
    real(dp), intent(in) :: data(:)
    type(sgt_params), intent(in) :: start
    type(sgt_mle_result), intent(out) :: result
    logical, intent(in), optional :: free(5)
    character(len=*), intent(in), optional :: methods(:)
    logical, intent(in), optional :: mean_cent, var_adj
    integer, intent(in), optional :: max_iter
    logical :: use_free(5), mc, va
    real(dp), allocatable :: theta0(:)
    integer :: nfree, j, k
    type(sgt_mle_result) :: tresult
    use_free = .true.
    if (present(free)) use_free = free
    if (.not. ieee_is_finite(start%p)) then
      use_free(3) = .false.
      use_free(4) = .false.
      use_free(5) = .false.
    else if (.not. ieee_is_finite(start%q)) then
      use_free(5) = .false.
    end if
    mc = .true.
    va = .true.
    if (present(mean_cent)) mc = mean_cent
    if (present(var_adj)) va = var_adj
    nfree = count(use_free)
    allocate(theta0(nfree))
    k = 0
    do j = 1, 5
      if (.not. use_free(j)) cycle
      k = k + 1
      select case (j)
      case (1)
        theta0(k) = start%mu
      case (2)
        theta0(k) = log(start%sigma)
      case (3)
        theta0(k) = atanh(max(-0.999999_dp, min(0.999999_dp, start%lambda)))
      case (4)
        theta0(k) = log(start%p)
      case (5)
        theta0(k) = log(start%q)
      end select
    end do
    if (allocated(constant_data_context)) deallocate(constant_data_context)
    allocate(constant_data_context(size(data)))
    constant_data_context = data
    constant_start_context = start
    constant_free_context = use_free
    if (present(methods)) then
      call sgt_mle_model(size(data), theta0, constant_context_model, tresult, methods, mc, va, max_iter)
    else
      call sgt_mle_model(size(data), theta0, constant_context_model, tresult, mean_cent=mc, &
        var_adj=va, max_iter=max_iter)
    end if
    call convert_constant_result(tresult, start, use_free, data, mc, va, result)
    if (allocated(constant_data_context)) deallocate(constant_data_context)
  end subroutine sgt_mle_constant

  subroutine convert_constant_result(tresult, start, free, data, mean_cent, var_adj, result)
    type(sgt_mle_result), intent(in) :: tresult
    type(sgt_params), intent(in) :: start
    logical, intent(in) :: free(5)
    real(dp), intent(in) :: data(:)
    logical, intent(in) :: mean_cent, var_adj
    type(sgt_mle_result), intent(out) :: result
    real(dp) :: raw(5)
    integer :: j, k, nf
    raw = [start%mu, start%sigma, start%lambda, start%p, start%q]
    k = 0
    do j = 1, 5
      if (.not. free(j)) cycle
      k = k + 1
      select case (j)
      case (1)
        raw(j) = tresult%estimate(k)
      case (2)
        raw(j) = exp(tresult%estimate(k))
      case (3)
        raw(j) = tanh(tresult%estimate(k))
      case (4,5)
        raw(j) = exp(tresult%estimate(k))
      end select
    end do
    nf = count(free)
    result%loglik = tresult%loglik
    result%convcode = tresult%convcode
    result%niter = tresult%niter
    result%best_method = tresult%best_method
    allocate(result%estimate(nf), result%gradient(nf), result%hessian(nf,nf))
    allocate(result%varcov(nf,nf), result%std_error(nf), result%z_score(nf), result%p_value(nf))
    k = 0
    do j = 1, 5
      if (.not. free(j)) cycle
      k = k + 1
      result%estimate(k) = raw(j)
    end do
    call constant_raw_inference(raw, free, data, mean_cent, var_adj, result)
  end subroutine convert_constant_result

  subroutine constant_raw_inference(raw, free, data, mean_cent, var_adj, result)
    real(dp), intent(in) :: raw(5), data(:)
    logical, intent(in) :: free(5), mean_cent, var_adj
    type(sgt_mle_result), intent(inout) :: result
    real(dp) :: xp(5), xm(5), xpp(5), xpm(5), xmp(5), xmm(5)
    real(dp), allocatable :: invh(:,:)
    real(dp) :: f0, fp, fm, hi, hj
    integer :: i, j, ii, jj, nf, status
    nf = count(free)
    f0 = raw_loglik(raw)
    ii = 0
    do i = 1, 5
      if (.not. free(i)) cycle
      ii = ii + 1
      hi = 1.0e-4_dp * max(1.0_dp, abs(raw(i)))
      if (i == 3) hi = min(hi, 0.1_dp * (1.0_dp - abs(raw(i))))
      if (i >= 2 .and. i /= 3) hi = min(hi, 0.2_dp * raw(i))
      xp = raw
      xm = raw
      xp(i) = xp(i) + hi
      xm(i) = xm(i) - hi
      fp = raw_loglik(xp)
      fm = raw_loglik(xm)
      result%gradient(ii) = (fp - fm) / (2.0_dp * hi)
      result%hessian(ii,ii) = (fp - 2.0_dp * f0 + fm) / (hi * hi)
      jj = 0
      do j = 1, i - 1
        if (.not. free(j)) cycle
        jj = jj + 1
        hj = 1.0e-4_dp * max(1.0_dp, abs(raw(j)))
        if (j == 3) hj = min(hj, 0.1_dp * (1.0_dp - abs(raw(j))))
        if (j >= 2 .and. j /= 3) hj = min(hj, 0.2_dp * raw(j))
        xpp = raw
        xpm = raw
        xmp = raw
        xmm = raw
        xpp(i) = xpp(i) + hi
        xpp(j) = xpp(j) + hj
        xpm(i) = xpm(i) + hi
        xpm(j) = xpm(j) - hj
        xmp(i) = xmp(i) - hi
        xmp(j) = xmp(j) + hj
        xmm(i) = xmm(i) - hi
        xmm(j) = xmm(j) - hj
        result%hessian(ii,jj) = (raw_loglik(xpp) - raw_loglik(xpm) - &
          raw_loglik(xmp) + raw_loglik(xmm)) / (4.0_dp * hi * hj)
        result%hessian(jj,ii) = result%hessian(ii,jj)
      end do
    end do
    call invert_matrix(-result%hessian, invh, status)
    if (status == 0) then
      result%varcov = invh
      do i = 1, nf
        if (result%varcov(i,i) > 0.0_dp) then
          result%std_error(i) = sqrt(result%varcov(i,i))
          result%z_score(i) = result%estimate(i) / result%std_error(i)
          result%p_value(i) = 2.0_dp * normal_cdf(-abs(result%z_score(i)))
        else
          result%std_error(i) = nan_dp()
          result%z_score(i) = nan_dp()
          result%p_value(i) = nan_dp()
        end if
      end do
    else
      result%varcov = nan_dp()
      result%std_error = nan_dp()
      result%z_score = nan_dp()
      result%p_value = nan_dp()
    end if
  contains
    real(dp) function raw_loglik(par) result(ll)
      real(dp), intent(in) :: par(5)
      integer :: obs
      real(dp) :: lf
      ll = 0.0_dp
      if (par(2) <= 0.0_dp .or. par(3) <= -1.0_dp .or. par(3) >= 1.0_dp .or. &
          par(4) <= 0.0_dp .or. par(5) <= 0.0_dp) then
        ll = -huge(1.0_dp) / 16.0_dp
        return
      end if
      do obs = 1, size(data)
        lf = sgt_logpdf(data(obs), par(1), par(2), par(3), par(4), par(5), mean_cent, var_adj)
        if (.not. ieee_is_finite(lf)) then
          ll = -huge(1.0_dp) / 16.0_dp
          return
        end if
        ll = ll + lf
      end do
    end function raw_loglik
  end subroutine constant_raw_inference
end module sgt_mle_mod
