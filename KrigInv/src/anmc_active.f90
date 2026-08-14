! SPDX-License-Identifier: GPL-3.0-only
module anmc_active
  use anmc_kinds, only : dp
  use anmc_types, only : active_dims_result
  use anmc_math, only : pmvnorm, probability_control, probability_result, genz_bretz
  use anmc_utils, only : normal_cdf_local, weighted_sample_without_replacement, weighted_sample_one, &
                         pairwise_distances, positive_infinity, negative_infinity, gather_vector, gather_matrix
  implicit none
  private
  public :: select_active_dims, select_q_dims

contains

  function select_active_dims(q, e, threshold, mu, sigma, pn, method) result(ind_q)
    integer, intent(in) :: q
    real(dp), intent(in) :: e(:,:), threshold, mu(:), sigma(:,:)
    real(dp), intent(in), optional :: pn(:)
    integer, intent(in), optional :: method
    integer, allocatable :: ind_q(:)

    real(dp), allocatable :: p(:), weights(:), distances(:,:), dd(:)
    integer :: n, m, i, k
    real(dp) :: den

    n = size(mu)
    m = 1
    if (present(method)) m = method
    if (q < 1 .or. q > n .or. size(e,1) /= n .or. size(sigma,1) /= n .or. size(sigma,2) /= n) then
      allocate(ind_q(0))
      return
    end if

    allocate(p(n))
    if (present(pn)) then
      if (size(pn) /= n) then
        allocate(ind_q(0))
        return
      end if
      p = max(0.0_dp, pn)
    else
      do i = 1, n
        if (sigma(i,i) > 0.0_dp) then
          p(i) = normal_cdf_local((mu(i)-threshold)/sqrt(sigma(i,i)))
        else
          p(i) = merge(1.0_dp, 0.0_dp, mu(i) > threshold)
        end if
      end do
    end if

    select case (m)
    case (0)
      allocate(ind_q(q))
      if (q == 1) then
        ind_q(1) = 1
      else
        do i = 1, q
          ! R's as.integer(seq(1,n,length.out=q)) truncates toward zero.
          ind_q(i) = int(1.0_dp + real(i-1,dp)*real(n-1,dp)/real(q-1,dp))
        end do
      end if
    case (1)
      ind_q = weighted_sample_without_replacement(p, q)
    case (2)
      ind_q = weighted_sample_without_replacement(p*(1.0_dp-p), q)
    case (3,4)
      distances = pairwise_distances(e)
      allocate(ind_q(q), dd(n), weights(n))
      if (m == 3) then
        weights = p
      else
        weights = p*(1.0_dp-p)
      end if
      k = weighted_sample_one(weights)
      if (k == 0) k = 1
      ind_q(1) = k
      dd = 1.0_dp
      do i = 2, q
        dd = dd**0.8_dp * distances(ind_q(i-1),:)
        den = maxval(dd) - minval(dd)
        if (den > tiny(1.0_dp)) then
          dd = dd / den
        else
          dd = 1.0_dp
        end if
        do k = 1, i-1
          dd(ind_q(k)) = 0.0_dp
        end do
        if (m == 3) then
          weights = dd*p
        else
          weights = dd*p*(1.0_dp-p)
        end if
        k = weighted_sample_one(weights)
        if (k == 0) then
          weights = 1.0_dp
          do k = 1, i-1
            weights(ind_q(k)) = 0.0_dp
          end do
          k = weighted_sample_one(weights)
        end if
        ind_q(i) = k
      end do
    case (5)
      allocate(weights(n)); weights = 1.0_dp
      ind_q = weighted_sample_without_replacement(weights, q)
    case default
      allocate(ind_q(q))
      if (q == 1) then
        ind_q(1) = 1
      else
        do i = 1, q
          ind_q(i) = int(1.0_dp + real(i-1,dp)*real(n-1,dp)/real(q-1,dp))
        end do
      end if
    end select

    call sort_integer(ind_q)
  end function select_active_dims

  function select_q_dims(e, threshold, mu, sigma, pn, method, limits, prob_control, reduced_return) result(res)
    real(dp), intent(in) :: e(:,:), threshold, mu(:), sigma(:,:)
    real(dp), intent(in), optional :: pn(:)
    integer, intent(in), optional :: method, limits(2)
    type(probability_control), intent(in), optional :: prob_control
    logical, intent(in), optional :: reduced_return
    type(active_dims_result) :: res

    type(probability_control) :: ctl, step_ctl
    type(probability_result) :: pr
    integer :: n, q0, qmax, qinc, q, m, flag
    integer, allocatable :: idx(:)
    real(dp), allocatable :: eq(:,:), meq(:), keq(:,:), lower(:), upper(:)
    real(dp) :: pprime, temp, err, delta_p
    logical :: reduced

    n = size(mu)
    m = 1
    if (present(method)) m = method
    reduced = .true.
    if (present(reduced_return)) reduced = reduced_return
    allocate(idx(0))
    ctl = genz_bretz()
    if (present(prob_control)) ctl = prob_control

    if (present(limits)) then
      q0 = max(2, min(limits(1), n))
      qmax = min(min(limits(2), n), 300)
    else
      q0 = min(10, n)
      qmax = min(n, 300)
    end if
    if (qmax < 1) then
      res%ok = .false.
      res%message = 'empty problem'
      allocate(res%ind_q(0))
      return
    end if
    q0 = min(q0, qmax)

    if (present(pn)) then
      idx = select_active_dims(q0,e,threshold,mu,sigma,pn,m)
    else
      idx = select_active_dims(q0,e,threshold,mu,sigma,method=m)
    end if
    call active_subproblem(idx,e,mu,sigma,eq,meq,keq)
    allocate(lower(q0), upper(q0))
    lower = negative_infinity(); upper = threshold
    pr = pmvnorm(lower,upper,meq,keq,ctl)
    pprime = 1.0_dp - pr%value
    err = pr%error
    delta_p = 1.0_dp
    qinc = min(10, ceiling(0.01_dp*real(n,dp)))
    qinc = max(1, qinc)
    flag = 0
    q = q0

    do while (flag < 2)
      q = min(q + qinc, qmax)
      if (present(pn)) then
        idx = select_active_dims(q,e,threshold,mu,sigma,pn,m)
      else
        idx = select_active_dims(q,e,threshold,mu,sigma,method=m)
      end if
      call active_subproblem(idx,e,mu,sigma,eq,meq,keq)
      if (allocated(lower)) deallocate(lower,upper)
      allocate(lower(q),upper(q)); lower=negative_infinity(); upper=threshold
      step_ctl = ctl
      step_ctl%abseps = 0.01_dp
      pr = pmvnorm(lower,upper,meq,keq,step_ctl)
      temp = 1.0_dp - pr%value
      err = pr%error
      delta_p = abs(temp-pprime)/(temp+1.0_dp)
      if (delta_p <= err) flag = flag + 1
      if (q == qmax) flag = flag + 2
      pprime = temp
    end do

    res%ind_q = idx
    res%pq = pprime
    res%error = err
    if (.not. reduced) then
      res%eq = eq
      res%mu_eq = meq
      res%k_eq = keq
      pr = pmvnorm(lower,upper,meq,keq,ctl)
      res%pq = 1.0_dp - pr%value
      res%error = pr%error
    end if
    res%ok = (pr%inform == 0 .or. pr%inform == 1)
    res%message = trim(pr%message)
  end function select_q_dims

  subroutine active_subproblem(idx,e,mu,sigma,eq,meq,keq)
    integer, intent(in) :: idx(:)
    real(dp), intent(in) :: e(:,:), mu(:), sigma(:,:)
    real(dp), allocatable, intent(out) :: eq(:,:), meq(:), keq(:,:)
    integer :: i
    eq = gather_matrix(e,idx,[(i,i=1,size(e,2))])
    meq = gather_vector(mu,idx)
    keq = gather_matrix(sigma,idx,idx)
  end subroutine active_subproblem

  subroutine sort_integer(x)
    integer, intent(inout) :: x(:)
    integer :: i,j,key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_integer

end module anmc_active
