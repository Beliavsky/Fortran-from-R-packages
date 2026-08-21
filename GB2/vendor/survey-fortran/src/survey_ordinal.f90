! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_ordinal
  use survey_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use survey_types, only : survey_design_t, olr_result_t, &
    OLR_LOGISTIC, OLR_PROBIT, OLR_CLOGLOG, OLR_CAUCHIT
  use survey_taylor, only : svyrecvar
  use survey_linalg, only : sym_pinv
  use minqa_module, only : minqa_result_t, minqa_control_t, uobyqa
  use numderiv, only : hessian, nd_success
  implicit none
  private
  public :: svy_olr, olr_predict_proba

  real(dp), allocatable :: active_x(:,:), active_w(:)
  integer, allocatable :: active_y(:)
  integer :: active_method = OLR_LOGISTIC
  integer :: active_k = 0

contains

  subroutine svy_olr(x, y, design, result, method, start, maxfun)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: y(:)
    type(survey_design_t), intent(in) :: design
    type(olr_result_t), intent(out) :: result
    integer, intent(in), optional :: method, maxfun
    real(dp), intent(in), optional :: start(:)
    type(minqa_result_t) :: opt
    type(minqa_control_t) :: ctrl
    real(dp), allocatable :: par(:), hh(:,:), vtheta(:,:), amat(:,:), &
      score(:,:), infl(:,:), probs(:,:), theta(:), zeta(:)
    integer :: n, p, q, kcat, i, j, rank, info, stat
    character(:), allocatable :: msg

    n = size(x,1)
    p = size(x,2)
    if (size(y) /= n .or. design%n /= n) error stop 'svy_olr: shape mismatch'
    if (n < 1) error stop 'svy_olr: no observations'
    kcat = maxval(y)
    if (minval(y) < 1 .or. kcat < 3) error stop 'svy_olr: y must use categories 1..K with K >= 3'
    if (any([(count(y == i) == 0, i=1,kcat)])) error stop 'svy_olr: every category must be observed'
    q = kcat - 1

    if (allocated(active_x)) error stop 'svy_olr: nested/concurrent calls are not supported'
    allocate(active_x, source=x)
    allocate(active_y, source=y)
    allocate(active_w, source=design%weight)
    active_method = OLR_LOGISTIC
    if (present(method)) active_method = method
    if (active_method < OLR_LOGISTIC .or. active_method > OLR_CAUCHIT) &
      error stop 'svy_olr: unsupported method'
    active_k = kcat

    allocate(par(p+q), theta(q), zeta(q))
    if (present(start)) then
      if (size(start) /= p+q) error stop 'svy_olr: start size mismatch'
      par = start
    else
      call ordinal_start(x, y, design%weight, active_method, par)
    end if

    ctrl%maxfun = 20000
    if (present(maxfun)) ctrl%maxfun = maxfun
    call uobyqa(ordinal_objective, par, opt, ctrl)
    par = opt%x

    call hessian(ordinal_objective, par, hh, status=stat, message=msg)
    if (stat /= nd_success) then
      call clear_active()
      error stop 'svy_olr: numerical Hessian failed'
    end if
    allocate(vtheta(p+q,p+q))
    call sym_pinv(hh, vtheta, rank, info=info)

    theta = par(p+1:p+q)
    call theta_to_zeta(theta, zeta)
    allocate(amat(p+q,p+q))
    amat = 0.0_dp
    do i = 1, p+q
      amat(i,i) = 1.0_dp
    end do
    if (q > 0) then
      amat(p+1:p+q,p+1:p+q) = 0.0_dp
      amat(p+1:p+q,p+1) = 1.0_dp
      do j = 2, q
        do i = j, q
          amat(p+i,p+j) = exp(theta(j))
        end do
      end do
    end if

    allocate(result%coef(p), result%zeta(q), result%vcov(p+q,p+q), &
      result%naive_vcov(p+q,p+q), result%fitted(n,kcat))
    result%coef = par(1:p)
    result%zeta = zeta
    result%naive_vcov = matmul(amat, matmul(vtheta, transpose(amat)))

    allocate(score(n,p+q), infl(n,p+q))
    call ordinal_scores_untransformed(par(1:p), zeta, score)
    infl = matmul(score, result%naive_vcov)
    result%vcov = svyrecvar(infl, design)

    allocate(probs(n,kcat))
    call olr_predict_proba(x, result%coef, result%zeta, active_method, probs)
    result%fitted = probs
    result%deviance = 2.0_dp*ordinal_objective(par)/(sum(design%weight)/real(n,dp))
    result%iterations = opt%evaluations
    result%converged = (opt%status == 0)
    result%method = active_method
    call clear_active()
  end subroutine svy_olr

  subroutine ordinal_start(x, y, w, method, par)
    real(dp), intent(in) :: x(:,:), w(:)
    integer, intent(in) :: y(:), method
    real(dp), intent(out) :: par(:)
    integer :: p, q, j
    real(dp) :: sw, cp, prev, z
    real(dp), allocatable :: zeta(:)
    p = size(x,2)
    q = maxval(y)-1
    par = 0.0_dp
    allocate(zeta(q))
    sw = sum(w)
    do j = 1, q
      cp = sum(w, mask=y <= j)/sw
      cp = max(1.0e-5_dp, min(1.0_dp-1.0e-5_dp, cp))
      zeta(j) = inverse_link_cdf(cp, method)
    end do
    do j = 2, q
      if (zeta(j) <= zeta(j-1)+1.0e-6_dp) zeta(j) = zeta(j-1)+0.25_dp
    end do
    par(p+1) = zeta(1)
    prev = zeta(1)
    do j = 2, q
      z = max(zeta(j)-prev, 1.0e-6_dp)
      par(p+j) = log(z)
      prev = zeta(j)
    end do
  end subroutine ordinal_start

  real(dp) function ordinal_objective(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp), allocatable :: zeta(:)
    real(dp) :: eta, lower, upper, pr
    integer :: i, p, q
    if (.not.allocated(active_x)) error stop 'survey_ordinal: no active fit'
    p = size(active_x,2)
    q = active_k-1
    allocate(zeta(q))
    call theta_to_zeta(par(p+1:p+q), zeta)
    value = 0.0_dp
    do i = 1, size(active_y)
      eta = 0.0_dp
      if (p > 0) eta = dot_product(active_x(i,:), par(1:p))
      if (active_y(i) == 1) then
        lower = 0.0_dp
      else
        lower = olr_cdf(zeta(active_y(i)-1)-eta, active_method)
      end if
      if (active_y(i) == active_k) then
        upper = 1.0_dp
      else
        upper = olr_cdf(zeta(active_y(i))-eta, active_method)
      end if
      pr = upper-lower
      if (pr <= tiny(1.0_dp) .or. ieee_is_nan(pr)) then
        value = huge(1.0_dp)/100.0_dp
        return
      end if
      value = value-active_w(i)*log(pr)
    end do
  end function ordinal_objective

  subroutine ordinal_scores_untransformed(beta, zeta, score)
    real(dp), intent(in) :: beta(:), zeta(:)
    real(dp), intent(out) :: score(:,:)
    real(dp) :: eta, fup, flo, pup, plo, pr
    integer :: n, p, q, i, j, yi
    n = size(active_y)
    p = size(beta)
    q = size(zeta)
    if (any(shape(score) /= [n,p+q])) error stop 'ordinal_scores_untransformed: shape mismatch'
    score = 0.0_dp
    do i = 1, n
      yi = active_y(i)
      eta = 0.0_dp
      if (p > 0) eta = dot_product(active_x(i,:), beta)
      if (yi == 1) then
        plo = 0.0_dp
        flo = 0.0_dp
      else
        plo = olr_cdf(zeta(yi-1)-eta, active_method)
        flo = olr_pdf(zeta(yi-1)-eta, active_method)
      end if
      if (yi == active_k) then
        pup = 1.0_dp
        fup = 0.0_dp
      else
        pup = olr_cdf(zeta(yi)-eta, active_method)
        fup = olr_pdf(zeta(yi)-eta, active_method)
      end if
      pr = max(pup-plo, tiny(1.0_dp))
      do j = 1, p
        score(i,j) = active_w(i)*active_x(i,j)*(fup-flo)/pr
      end do
      if (yi > 1) score(i,p+yi-1) = score(i,p+yi-1)+active_w(i)*flo/pr
      if (yi <= q) score(i,p+yi) = score(i,p+yi)-active_w(i)*fup/pr
    end do
  end subroutine ordinal_scores_untransformed

  subroutine theta_to_zeta(theta, zeta)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: zeta(:)
    integer :: j
    if (size(theta) /= size(zeta)) error stop 'theta_to_zeta: shape mismatch'
    if (size(theta) == 0) return
    zeta(1) = theta(1)
    do j = 2, size(theta)
      zeta(j) = zeta(j-1)+exp(theta(j))
    end do
  end subroutine theta_to_zeta

  subroutine olr_predict_proba(x, beta, zeta, method, prob)
    real(dp), intent(in) :: x(:,:), beta(:), zeta(:)
    integer, intent(in) :: method
    real(dp), intent(out) :: prob(:,:)
    real(dp) :: eta, prev, cur
    integer :: i, j, n, k
    n = size(x,1)
    k = size(zeta)+1
    if (size(x,2) /= size(beta) .or. any(shape(prob) /= [n,k])) &
      error stop 'olr_predict_proba: shape mismatch'
    do i = 1, n
      eta = 0.0_dp
      if (size(beta) > 0) eta = dot_product(x(i,:), beta)
      prev = 0.0_dp
      do j = 1, k-1
        cur = olr_cdf(zeta(j)-eta, method)
        prob(i,j) = max(0.0_dp, cur-prev)
        prev = cur
      end do
      prob(i,k) = max(0.0_dp, 1.0_dp-prev)
      if (sum(prob(i,:)) > 0.0_dp) prob(i,:) = prob(i,:)/sum(prob(i,:))
    end do
  end subroutine olr_predict_proba

  real(dp) function olr_cdf(z, method) result(p)
    real(dp), intent(in) :: z
    integer, intent(in) :: method
    real(dp), parameter :: pi = acos(-1.0_dp)
    select case(method)
    case(OLR_LOGISTIC)
      if (z >= 0.0_dp) then
        p = 1.0_dp/(1.0_dp+exp(-min(z,700.0_dp)))
      else
        p = exp(max(z,-700.0_dp))/(1.0_dp+exp(max(z,-700.0_dp)))
      end if
    case(OLR_PROBIT)
      p = 0.5_dp*erfc(-z/sqrt(2.0_dp))
    case(OLR_CLOGLOG)
      if (z < -40.0_dp) then
        p = 0.0_dp
      else if (z > 40.0_dp) then
        p = 1.0_dp
      else
        p = exp(-exp(-z))
      end if
    case(OLR_CAUCHIT)
      p = 0.5_dp+atan(z)/pi
    case default
      error stop 'olr_cdf: unsupported method'
    end select
    p = max(0.0_dp,min(1.0_dp,p))
  end function olr_cdf

  real(dp) function olr_pdf(z, method) result(d)
    real(dp), intent(in) :: z
    integer, intent(in) :: method
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: p
    select case(method)
    case(OLR_LOGISTIC)
      p = olr_cdf(z,method)
      d = p*(1.0_dp-p)
    case(OLR_PROBIT)
      d = exp(-0.5_dp*z*z)/sqrt(2.0_dp*pi)
    case(OLR_CLOGLOG)
      if (abs(z) > 40.0_dp) then
        d = 0.0_dp
      else
        d = exp(-z-exp(-z))
      end if
    case(OLR_CAUCHIT)
      d = 1.0_dp/(pi*(1.0_dp+z*z))
    case default
      error stop 'olr_pdf: unsupported method'
    end select
  end function olr_pdf

  real(dp) function inverse_link_cdf(p, method) result(z)
    real(dp), intent(in) :: p
    integer, intent(in) :: method
    real(dp) :: lo, hi, mid
    integer :: it
    lo = -50.0_dp
    hi = 50.0_dp
    do it = 1, 160
      mid = 0.5_dp*(lo+hi)
      if (olr_cdf(mid,method) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    z = 0.5_dp*(lo+hi)
  end function inverse_link_cdf

  subroutine clear_active()
    if (allocated(active_x)) deallocate(active_x)
    if (allocated(active_y)) deallocate(active_y)
    if (allocated(active_w)) deallocate(active_w)
    active_k = 0
  end subroutine clear_active

end module survey_ordinal
