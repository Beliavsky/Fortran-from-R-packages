! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_pi
  use peerperformance_kinds, only: dp
  use peerperformance_math, only: normal_pdf, normal_cdf, normal_quantile, &
                                  finite_value, missing_value, set_random_seed, &
                                  random_integer, clamp_probability
  implicit none
  private
  public :: adjust_pi, compute_pizero, optimal_lambda, compute_peer_ratios

contains

  pure real(dp) function forward_adjustment(pi0, n, lambda) result(value)
    real(dp), intent(in) :: pi0, lambda
    integer, intent(in) :: n
    real(dp) :: nlambda, s, zcrit
    nlambda = pi0*real(n,dp)*(1.0_dp-lambda)
    if (nlambda >= real(n,dp) .or. nlambda <= 0.0_dp) then
      value = pi0
      return
    end if
    s = sqrt(nlambda*(real(n,dp)-nlambda) / &
        (real(n,dp)**3*(1.0_dp-lambda)**2))
    if (s <= tiny(1.0_dp)) then
      value = pi0
      return
    end if
    zcrit = (1.0_dp-pi0)/s
    value = pi0+s*(-normal_pdf(zcrit)+(1.0_dp-normal_cdf(zcrit))*zcrit)
  end function forward_adjustment

  pure real(dp) function adjust_pi(pi_hat, n, lambda, fast) result(value)
    real(dp), intent(in) :: pi_hat, lambda
    integer, intent(in) :: n
    logical, intent(in), optional :: fast
    real(dp) :: lo, hi, mid, flo, fhi
    integer :: iter, niter
    logical :: accurate
    if (.not. finite_value(pi_hat) .or. n < 1 .or. lambda < 0.0_dp .or. lambda >= 1.0_dp) then
      value = pi_hat
      return
    end if
    accurate = .true.
    if (present(fast)) accurate = fast
    lo = 1.0e-5_dp
    hi = 1.5_dp
    flo = forward_adjustment(lo,n,lambda)
    fhi = forward_adjustment(hi,n,lambda)
    if (pi_hat < flo .or. pi_hat > fhi) then
      value = clamp_probability(pi_hat)
      return
    end if
    if (accurate) then
      niter = 50
    else
      niter = 15
    end if
    do iter = 1, niter
      mid = 0.5_dp*(lo+hi)
      if (forward_adjustment(mid,n,lambda) < pi_hat) then
        lo = mid
      else
        hi = mid
      end if
    end do
    value = clamp_probability(0.5_dp*(lo+hi))
  end function adjust_pi

  pure real(dp) function compute_pizero(pvalue, lambda, adjust, fast, n_trials) result(value)
    real(dp), intent(in) :: pvalue(:), lambda
    logical, intent(in), optional :: adjust, fast
    integer, intent(in), optional :: n_trials
    logical :: use_adjust, accurate
    integer :: nvalid, n
    use_adjust = .true.; if (present(adjust)) use_adjust = adjust
    accurate = .true.; if (present(fast)) accurate = fast
    nvalid = count(finite_value(pvalue))
    if (nvalid <= 0 .or. lambda < 0.0_dp .or. lambda >= 1.0_dp) then
      value = missing_value()
      return
    end if
    value = real(count(finite_value(pvalue) .and. pvalue >= lambda),dp) / &
            real(nvalid,dp)/(1.0_dp-lambda)
    value = min(1.0_dp,value)
    n = size(pvalue); if (present(n_trials)) n = n_trials
    if (use_adjust) value = adjust_pi(value,n,lambda,accurate)
  end function compute_pizero

  real(dp) function optimal_lambda(pvalue, n_boot, seed, adjust, fast) result(value)
    real(dp), intent(in) :: pvalue(:)
    integer, intent(in), optional :: n_boot, seed
    logical, intent(in), optional :: adjust, fast
    real(dp), parameter :: grid(5) = [0.3_dp,0.4_dp,0.5_dp,0.6_dp,0.7_dp]
    real(dp), allocatable :: p(:), sample(:)
    real(dp) :: estimates(5), mse(5), target
    integer :: nb, sd, n, i, j, b, best
    logical :: use_adjust, accurate
    nb = 499; if (present(n_boot)) nb = n_boot
    sd = 12345; if (present(seed)) sd = seed
    use_adjust = .true.; if (present(adjust)) use_adjust = adjust
    accurate = .true.; if (present(fast)) accurate = fast
    p = pack(pvalue,finite_value(pvalue))
    n = size(p)
    if (n <= 0) then
      value = 0.5_dp
      return
    end if
    do j = 1, 5
      estimates(j) = compute_pizero(p,grid(j),use_adjust,accurate,n)
    end do
    target = minval(estimates)
    mse = 0.0_dp
    allocate(sample(n))
    call set_random_seed(sd)
    do b = 1, nb
      do i = 1, n
        sample(i) = p(random_integer(n))
      end do
      do j = 1, 5
        mse(j) = mse(j)+(compute_pizero(sample,grid(j),use_adjust,accurate,n)-target)**2
      end do
    end do
    best = 1
    do j = 2, 5
      if (mse(j) < mse(best)) best = j
    end do
    value = grid(best)
  end function optimal_lambda

  pure integer function half_to_even_count(n) result(value)
    integer, intent(in) :: n
    integer :: lower
    lower = n/2
    if (mod(n,2) == 0) then
      value = lower
    else if (mod(lower,2) == 0) then
      value = lower
    else
      value = lower+1
    end if
  end function half_to_even_count

  subroutine compute_peer_ratios(pvalue, difference, tstat, pizero, pipos, pineg, &
                                 lambda_used, lambda, n_boot, gamma_pos, gamma_neg, &
                                 seed, adjust, fast)
    real(dp), intent(in) :: pvalue(:,:), difference(:,:), tstat(:,:)
    real(dp), intent(out) :: pizero(size(pvalue,1)), pipos(size(pvalue,1))
    real(dp), intent(out) :: pineg(size(pvalue,1)), lambda_used(size(pvalue,1))
    real(dp), intent(in), optional :: lambda(:), gamma_pos, gamma_neg
    integer, intent(in), optional :: n_boot, seed
    logical, intent(in), optional :: adjust, fast
    logical :: mask(size(pvalue,2)), use_adjust, accurate
    real(dp) :: lp, ln, qpos, qneg, ni0
    integer :: i, n, hn, nb, sd
    use_adjust=.true.; if (present(adjust)) use_adjust=adjust
    accurate=.true.; if (present(fast)) accurate=fast
    nb=499; if (present(n_boot)) nb=n_boot
    sd=12345; if (present(seed)) sd=seed
    lp=0.4_dp; if (present(gamma_pos)) lp=gamma_pos
    ln=0.6_dp; if (present(gamma_neg)) ln=gamma_neg
    qpos=normal_quantile(lp); qneg=normal_quantile(ln)
    pizero=missing_value(); pipos=missing_value(); pineg=missing_value(); lambda_used=missing_value()
    do i=1,size(pvalue,1)
      mask=finite_value(pvalue(i,:)) .and. finite_value(difference(i,:)) .and. finite_value(tstat(i,:))
      n=count(mask)
      if (n <= 1) cycle
      if (present(lambda)) then
        if (size(lambda)==1) then
          lambda_used(i)=lambda(1)
        else if (size(lambda)==size(pvalue,1)) then
          lambda_used(i)=lambda(i)
        else
          cycle
        end if
      else
        lambda_used(i)=optimal_lambda(pack(pvalue(i,:),finite_value(pvalue(i,:))),nb,sd+i,use_adjust,accurate)
      end if
      pizero(i)=compute_pizero(pack(pvalue(i,:),finite_value(pvalue(i,:))),lambda_used(i), &
                               use_adjust,accurate,size(pvalue,2))
      ni0=pizero(i)*real(n,dp)
      hn=half_to_even_count(n)
      pipos(i)=0.0_dp; pineg(i)=0.0_dp
      if (count(mask .and. difference(i,:) >= 0.0_dp) >= hn) then
        pipos(i)=min(real(n,dp)-ni0,max(real(count(mask .and. tstat(i,:) >= qpos),dp)- &
                   ni0*(1.0_dp-lp),0.0_dp))/real(n,dp)
        pineg(i)=clamp_probability(1.0_dp-pizero(i)-pipos(i))
      else
        pineg(i)=min(real(n,dp)-ni0,max(real(count(mask .and. tstat(i,:) <= qneg),dp)- &
                   ni0*ln,0.0_dp))/real(n,dp)
        pipos(i)=clamp_probability(1.0_dp-pizero(i)-pineg(i))
      end if
    end do
  end subroutine compute_peer_ratios

end module peerperformance_pi
