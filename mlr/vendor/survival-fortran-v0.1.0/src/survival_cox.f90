! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_cox
  use survival_kinds, only : dp
  use survival_types, only : coxph_result
  use survival_linalg, only : solve_sym, invert_matrix, matrix_rank
  implicit none
  private
  public :: coxph_fit, coxph_fit_counting
  public :: cox_baseline, cox_martingale_residuals, cox_schoenfeld_residuals

contains

  subroutine coxph_fit(time, status, x, result, method, weights, offset, maxiter, eps)
    real(dp), intent(in) :: time(:), x(:,:)
    integer, intent(in) :: status(:)
    type(coxph_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: weights(:), offset(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: eps
    real(dp), allocatable :: start(:)

    allocate(start(size(time)))
    start = -huge(1.0_dp)
    call coxph_fit_counting(start, time, status, x, result, method, &
                            weights, offset, maxiter, eps)
  end subroutine coxph_fit

  subroutine coxph_fit_counting(start, stop, status, x, result, method, &
                                weights, offset, maxiter, eps)
    real(dp), intent(in) :: start(:), stop(:), x(:,:)
    integer, intent(in) :: status(:)
    type(coxph_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: weights(:), offset(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: eps

    integer :: n, p, iter, max_it
    real(dp) :: tolerance, ll, newll, alpha
    real(dp), allocatable :: beta(:), score(:), info(:,:), step(:)
    real(dp), allocatable :: w(:), off(:), means(:), vinv(:,:)
    logical :: ok
    character(len=12) :: tie_method

    n = size(stop)
    p = size(x,2)
    max_it = 20
    if (present(maxiter)) max_it = maxiter
    tolerance = 1.0e-9_dp
    if (present(eps)) tolerance = eps
    tie_method = 'efron'
    if (present(method)) tie_method = adjustl(method)

    allocate(beta(p), score(p), info(p,p), step(p), w(n), off(n))
    allocate(means(p), vinv(p,p))
    beta = 0.0_dp
    if (present(weights)) then
      w = weights
    else
      w = 1.0_dp
    end if
    if (present(offset)) then
      off = offset
    else
      off = 0.0_dp
    end if
    means = matmul(transpose(x), w) / sum(w)

    call cox_stats(start, stop, status, x, beta, w, off, tie_method, &
                   ll, score, info)
    result%loglik_initial = ll
    call solve_sym(info, score, step, ok)
    if (ok) result%score_test = dot_product(score, step)

    result%converged = .false.
    iter = 0
    do iter = 1, max_it
      call solve_sym(info, score, step, ok)
      if (.not. ok) exit
      alpha = 1.0_dp
      do
        call cox_stats(start, stop, status, x, beta + alpha*step, w, off, &
                       tie_method, newll, score, info)
        if (newll >= ll .or. alpha < 1.0e-8_dp) exit
        alpha = alpha / 2.0_dp
      end do
      beta = beta + alpha * step
      if (abs(newll - ll) <= tolerance * (1.0_dp + abs(newll))) then
        result%converged = .true.
        ll = newll
        exit
      end if
      ll = newll
      call cox_stats(start, stop, status, x, beta, w, off, tie_method, &
                     ll, score, info)
    end do

    call cox_stats(start, stop, status, x, beta, w, off, tie_method, &
                   ll, score, info)
    call invert_matrix(info, vinv, ok)
    if (.not. ok) vinv = 0.0_dp

    allocate(result%coef(p), result%var(p,p), result%score(p), result%means(p))
    result%coef = beta
    result%var = vinv
    result%score = score
    result%means = means
    result%loglik = ll
    result%iterations = min(iter, max_it)
    result%rank = matrix_rank(info)
  end subroutine coxph_fit_counting

  subroutine cox_stats(start, stop, status, x, beta, w, off, method, &
                       loglik, score, info)
    real(dp), intent(in) :: start(:), stop(:), x(:,:), beta(:), w(:), off(:)
    integer, intent(in) :: status(:)
    character(len=*), intent(in) :: method
    real(dp), intent(out) :: loglik, score(:), info(:,:)

    integer :: n, p, i, j, k, ndeath
    real(dp) :: t, denom, event_denom, event_weight, frac, eta, rr, part
    real(dp), allocatable :: event_times(:), a(:), a2(:,:), d1(:), d2(:,:)
    real(dp), allocatable :: meanx(:)

    n = size(stop)
    p = size(beta)
    allocate(event_times(n), a(p), a2(p,p), d1(p), d2(p,p), meanx(p))
    k = 0
    do i = 1, n
      if (status(i) /= 0 .and. .not. contains_time(event_times, k, stop(i))) then
        k = k + 1
        event_times(k) = stop(i)
      end if
    end do
    call sort_real(event_times(1:k))

    loglik = 0.0_dp
    score = 0.0_dp
    info = 0.0_dp
    do j = 1, k
      t = event_times(j)
      denom = 0.0_dp
      a = 0.0_dp
      a2 = 0.0_dp
      event_weight = 0.0_dp
      event_denom = 0.0_dp
      d1 = 0.0_dp
      d2 = 0.0_dp
      ndeath = 0

      do i = 1, n
        eta = dot_product(x(i,:), beta) + off(i)
        rr = w(i) * exp(min(eta, 700.0_dp))
        if (start(i) < t .and. stop(i) >= t) then
          denom = denom + rr
          a = a + rr * x(i,:)
          a2 = a2 + rr * outer(x(i,:), x(i,:))
        end if
        if (status(i) /= 0 .and. same_time(stop(i), t)) then
          ndeath = ndeath + 1
          event_weight = event_weight + w(i)
          loglik = loglik + w(i) * eta
          score = score + w(i) * x(i,:)
          event_denom = event_denom + rr
          d1 = d1 + rr * x(i,:)
          d2 = d2 + rr * outer(x(i,:), x(i,:))
        end if
      end do

      if (ndeath == 0 .or. denom <= 0.0_dp) cycle
      if (index(method, 'efron') == 1 .and. ndeath > 1) then
        do i = 0, ndeath - 1
          frac = real(i,dp) / real(ndeath,dp)
          part = denom - frac * event_denom
          meanx = (a - frac*d1) / part
          loglik = loglik - (event_weight/real(ndeath,dp)) * log(part)
          score = score - (event_weight/real(ndeath,dp)) * meanx
          info = info + (event_weight/real(ndeath,dp)) * &
                 ((a2-frac*d2)/part - outer(meanx,meanx))
        end do
      else
        meanx = a / denom
        loglik = loglik - event_weight * log(denom)
        score = score - event_weight * meanx
        info = info + event_weight * (a2/denom - outer(meanx,meanx))
      end if
    end do
  end subroutine cox_stats

  subroutine cox_baseline(time, status, x, beta, base_time, cumhaz, method, weights)
    real(dp), intent(in) :: time(:), x(:,:), beta(:)
    integer, intent(in) :: status(:)
    real(dp), allocatable, intent(out) :: base_time(:), cumhaz(:)
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: weights(:)

    real(dp), allocatable :: w(:), unique_time(:)
    real(dp) :: denom, event_weight, event_denom, inc, frac
    integer :: n, m, i, j, k, ndeath
    character(len=12) :: tie_method

    n = size(time)
    allocate(w(n), unique_time(n))
    if (present(weights)) then
      w = weights
    else
      w = 1.0_dp
    end if
    tie_method = 'breslow'
    if (present(method)) tie_method = method

    m = 0
    do i = 1, n
      if (status(i) /= 0 .and. .not. contains_time(unique_time,m,time(i))) then
        m = m + 1
        unique_time(m) = time(i)
      end if
    end do
    call sort_real(unique_time(1:m))
    allocate(base_time(m), cumhaz(m))
    base_time = unique_time(1:m)
    cumhaz = 0.0_dp

    do k = 1, m
      denom = 0.0_dp
      event_weight = 0.0_dp
      event_denom = 0.0_dp
      ndeath = 0
      do i = 1, n
        if (time(i) >= base_time(k)) then
          denom = denom + w(i) * exp(min(dot_product(x(i,:),beta),700.0_dp))
        end if
        if (status(i) /= 0 .and. same_time(time(i),base_time(k))) then
          event_weight = event_weight + w(i)
          event_denom = event_denom + w(i) * &
                        exp(min(dot_product(x(i,:),beta),700.0_dp))
          ndeath = ndeath + 1
        end if
      end do
      inc = 0.0_dp
      if (index(tie_method,'efron') == 1 .and. ndeath > 1) then
        do j = 0, ndeath - 1
          frac = real(j,dp) / real(ndeath,dp)
          inc = inc + (event_weight/real(ndeath,dp)) / &
                (denom - frac*event_denom)
        end do
      else if (denom > 0.0_dp) then
        inc = event_weight / denom
      end if
      if (k == 1) then
        cumhaz(k) = inc
      else
        cumhaz(k) = cumhaz(k-1) + inc
      end if
    end do
  end subroutine cox_baseline

  subroutine cox_martingale_residuals(time, status, x, beta, resid, weights)
    real(dp), intent(in) :: time(:), x(:,:), beta(:)
    integer, intent(in) :: status(:)
    real(dp), intent(out) :: resid(:)
    real(dp), intent(in), optional :: weights(:)

    real(dp), allocatable :: base_time(:), cumhaz(:), w(:)
    real(dp) :: hazard
    integer :: i, k, n

    n = size(time)
    allocate(w(n))
    if (present(weights)) then
      w = weights
    else
      w = 1.0_dp
    end if
    call cox_baseline(time, status, x, beta, base_time, cumhaz, 'breslow', w)
    do i = 1, n
      hazard = 0.0_dp
      do k = 1, size(base_time)
        if (base_time(k) <= time(i)) hazard = cumhaz(k)
      end do
      resid(i) = w(i)*real(status(i),dp) - &
                 w(i)*exp(dot_product(x(i,:),beta))*hazard
    end do
  end subroutine cox_martingale_residuals

  subroutine cox_schoenfeld_residuals(time, status, x, beta, resid)
    real(dp), intent(in) :: time(:), x(:,:), beta(:)
    integer, intent(in) :: status(:)
    real(dp), allocatable, intent(out) :: resid(:,:)

    integer :: n, p, i, j, k, nevent
    real(dp) :: denom, rr
    real(dp), allocatable :: average(:)

    n = size(time)
    p = size(beta)
    nevent = count(status /= 0)
    allocate(resid(nevent,p), average(p))
    k = 0
    do i = 1, n
      if (status(i) == 0) cycle
      k = k + 1
      denom = 0.0_dp
      average = 0.0_dp
      do j = 1, n
        if (time(j) >= time(i)) then
          rr = exp(dot_product(x(j,:), beta))
          denom = denom + rr
          average = average + rr * x(j,:)
        end if
      end do
      if (denom > 0.0_dp) average = average / denom
      resid(k,:) = x(i,:) - average
    end do
  end subroutine cox_schoenfeld_residuals

  pure function outer(a,b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a),size(b))
    c = spread(a,2,size(b)) * spread(b,1,size(a))
  end function outer

  logical function contains_time(x,nused,value) result(found)
    real(dp), intent(in) :: x(:), value
    integer, intent(in) :: nused
    integer :: i
    found = .false.
    do i = 1, nused
      if (same_time(x(i),value)) then
        found = .true.
        return
      end if
    end do
  end function contains_time

  pure logical function same_time(a,b) result(equal)
    real(dp), intent(in) :: a,b
    real(dp) :: scale
    scale = max(1.0_dp,abs(a),abs(b))
    equal = abs(a-b) <= 8.0_dp*epsilon(1.0_dp)*scale
  end function same_time

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: value
    do i=2,size(x)
      value=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=value) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=value
    end do
  end subroutine sort_real

end module survival_cox
