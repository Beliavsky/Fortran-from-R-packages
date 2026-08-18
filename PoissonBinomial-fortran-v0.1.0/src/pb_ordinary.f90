! SPDX-License-Identifier: GPL-3.0-only
module pb_ordinary
  use pb_kinds, only : dp, pi
  use pb_math, only : normal_cdf, normal_pdf, binom_pmf, binom_cdf, &
                      poisson_pmf, poisson_cdf
  use pb_numerics, only : normalize_pmf, convolve_fft, expand_probs, sample_from_pmf
  implicit none
  private
  public :: dpbinom, dpbinom_at, dpbinom_values
  public :: ppbinom, ppbinom_at, ppbinom_values
  public :: qpbinom, qpbinom_values, rpbinom
  public :: dpb_convolve, dpb_dividefft, dpb_characteristic, dpb_recursive
  public :: dpb_mean, dpb_geomean, dpb_geomean_counter, dpb_poisson
  public :: dpb_normal, ppb_normal

contains

  recursive function divide_pmf(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:), left(:), right(:)
    integer :: n, m
    n = size(probs)
    if (n == 0) then
      allocate(pmf(0:0))
      pmf(0) = 1.0_dp
    else if (n == 1) then
      allocate(pmf(0:1))
      pmf(0) = 1.0_dp - probs(1)
      pmf(1) = probs(1)
    else if (n <= 32) then
      pmf = dpb_convolve(probs)
    else
      m = n/2
      left = divide_pmf(probs(:m))
      right = divide_pmf(probs(m+1:))
      pmf = convolve_fft(left,right)
      call normalize_pmf(pmf)
    end if
  end function divide_pmf

  function dpb_convolve(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    integer :: i, k, n
    if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
      error stop "probabilities must lie in [0,1]"
    end if
    n = size(probs)
    allocate(pmf(0:n))
    pmf = 0.0_dp
    pmf(0) = 1.0_dp
    do i = 1, n
      do k = i, 1, -1
        pmf(k) = pmf(k)*(1.0_dp-probs(i)) + pmf(k-1)*probs(i)
      end do
      pmf(0) = pmf(0)*(1.0_dp-probs(i))
    end do
    call normalize_pmf(pmf)
  end function dpb_convolve

  function dpb_dividefft(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
      error stop "probabilities must lie in [0,1]"
    end if
    pmf = divide_pmf(probs)
  end function dpb_dividefft

  function dpb_characteristic(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    complex(dp), allocatable :: phi(:)
    complex(dp) :: z, prod, accum
    integer :: n, nn, j, l, k
    real(dp) :: theta
    if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
      error stop "probabilities must lie in [0,1]"
    end if
    n = size(probs)
    nn = n + 1
    allocate(phi(0:n), pmf(0:n))
    do l = 0, n
      theta = 2.0_dp*pi*real(l,dp)/real(nn,dp)
      z = cmplx(cos(theta),sin(theta),kind=dp)
      prod = cmplx(1.0_dp,0.0_dp,kind=dp)
      do j = 1, n
        prod = prod*((1.0_dp-probs(j)) + probs(j)*z)
      end do
      phi(l) = prod
    end do
    do k = 0, n
      accum = cmplx(0.0_dp,0.0_dp,kind=dp)
      do l = 0, n
        theta = -2.0_dp*pi*real(k*l,dp)/real(nn,dp)
        accum = accum + phi(l)*cmplx(cos(theta),sin(theta),kind=dp)
      end do
      pmf(k) = real(accum,dp)/real(nn,dp)
    end do
    where (pmf < 2.22e-16_dp) pmf = 0.0_dp
    where (pmf > 1.0_dp) pmf = 1.0_dp
    call normalize_pmf(pmf)
  end function dpb_characteristic

  function dpb_recursive(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:), dist(:,:)
    integer :: n, i, j, col_new, col_old, tmp
    if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
      error stop "probabilities must lie in [0,1]"
    end if
    n = size(probs)
    allocate(pmf(0:n))
    if (n == 0) then
      pmf(0) = 1.0_dp
      return
    end if
    allocate(dist(0:n,2))
    dist = 0.0_dp
    col_new = 1
    col_old = 2
    dist(0,col_new) = 1.0_dp
    do j = 1, n
      dist(j,col_new) = (1.0_dp-probs(j))*dist(j-1,col_new)
    end do
    pmf(0) = dist(n,col_new)
    do i = 1, n
      tmp = col_old
      col_old = col_new
      col_new = tmp
      dist(:,col_new) = 0.0_dp
      do j = i, n
        dist(j,col_new) = (1.0_dp-probs(j))*dist(j-1,col_new) + &
                          probs(j)*dist(j-1,col_old)
      end do
      pmf(i) = dist(n,col_new)
    end do
    call normalize_pmf(pmf)
  end function dpb_recursive

  function binomial_approx(n, p) result(pmf)
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    real(dp), allocatable :: pmf(:)
    integer :: k
    allocate(pmf(0:n))
    do k = 0, n
      pmf(k) = binom_pmf(k,n,p)
    end do
    call normalize_pmf(pmf)
  end function binomial_approx

  function dpb_mean(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    real(dp) :: p
    if (size(probs) == 0) then
      allocate(pmf(0:0)); pmf(0)=1.0_dp; return
    end if
    p = sum(probs)/real(size(probs),dp)
    pmf = binomial_approx(size(probs),p)
  end function dpb_mean

  function dpb_geomean(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    real(dp) :: p
    if (size(probs) == 0) then
      allocate(pmf(0:0)); pmf(0)=1.0_dp; return
    end if
    if (any(probs <= 0.0_dp)) then
      p = 0.0_dp
    else
      p = exp(sum(log(probs))/real(size(probs),dp))
    end if
    pmf = binomial_approx(size(probs),p)
  end function dpb_geomean

  function dpb_geomean_counter(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    real(dp) :: p
    if (size(probs) == 0) then
      allocate(pmf(0:0)); pmf(0)=1.0_dp; return
    end if
    if (any(probs >= 1.0_dp)) then
      p = 1.0_dp
    else
      p = 1.0_dp-exp(sum(log(1.0_dp-probs))/real(size(probs),dp))
    end if
    pmf = binomial_approx(size(probs),p)
  end function dpb_geomean_counter

  function dpb_poisson(probs) result(pmf)
    real(dp), intent(in) :: probs(:)
    real(dp), allocatable :: pmf(:)
    integer :: n, k
    real(dp) :: lambda
    n = size(probs)
    lambda = sum(probs)
    allocate(pmf(0:n))
    if (n == 0) then
      pmf(0) = 1.0_dp
      return
    end if
    do k = 0, n-1
      pmf(k) = poisson_pmf(k,lambda)
    end do
    pmf(n) = poisson_pmf(n,lambda) + (1.0_dp-poisson_cdf(n,lambda))
    call normalize_pmf(pmf)
  end function dpb_poisson

  function ppb_normal(probs, refined, lower_tail) result(cdf)
    real(dp), intent(in) :: probs(:)
    logical, intent(in), optional :: refined, lower_tail
    real(dp), allocatable :: cdf(:)
    logical :: r, lower
    integer :: n, k
    real(dp) :: mu, sigma, gamma3, z, val
    real(dp), allocatable :: pq(:)
    r = .true.; if (present(refined)) r = refined
    lower = .true.; if (present(lower_tail)) lower = lower_tail
    n = size(probs)
    allocate(cdf(0:n))
    if (n == 0) then
      cdf(0) = merge(1.0_dp,0.0_dp,lower)
      return
    end if
    mu = sum(probs)
    allocate(pq(n)); pq = probs*(1.0_dp-probs)
    sigma = sqrt(sum(pq))
    if (sigma <= tiny(1.0_dp)) then
      do k = 0, n
        if (lower) then
          cdf(k) = merge(1.0_dp,0.0_dp,real(k,dp) >= mu)
        else
          cdf(k) = merge(0.0_dp,1.0_dp,real(k,dp) >= mu)
        end if
      end do
      return
    end if
    gamma3 = sum(pq*(1.0_dp-2.0_dp*probs))
    do k = 0, n
      z = (real(k,dp)+0.5_dp-mu)/sigma
      val = normal_cdf(z)
      if (.not. lower) val = 1.0_dp-val
      if (r) then
        if (lower) then
          val = val + gamma3/(6.0_dp*sigma**3)*(1.0_dp-z*z)*normal_pdf(z)
        else
          val = val - gamma3/(6.0_dp*sigma**3)*(1.0_dp-z*z)*normal_pdf(z)
        end if
      end if
      cdf(k) = max(0.0_dp,min(1.0_dp,val))
    end do
    cdf(n) = merge(1.0_dp,0.0_dp,lower)
  end function ppb_normal

  function dpb_normal(probs, refined) result(pmf)
    real(dp), intent(in) :: probs(:)
    logical, intent(in), optional :: refined
    real(dp), allocatable :: pmf(:), cl(:), cu(:)
    integer :: n, k, mid
    logical :: r
    r = .true.; if (present(refined)) r = refined
    n = size(probs)
    allocate(pmf(0:n))
    if (n == 0) then
      pmf(0)=1.0_dp; return
    end if
    cl = ppb_normal(probs,r,.true.)
    cu = ppb_normal(probs,r,.false.)
    mid = floor(sum(probs)+0.5_dp)
    pmf(0) = cl(0)
    do k = 1, n
      if (k <= mid) then
        pmf(k) = cl(k)-cl(k-1)
      else
        pmf(k) = cu(k-1)-cu(k)
      end if
    end do
    where (pmf < 0.0_dp .and. pmf > -1.0e-12_dp) pmf = 0.0_dp
    call normalize_pmf(pmf)
  end function dpb_normal

  function reduced_method_pmf(probs, method) result(pmf)
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in) :: method
    real(dp), allocatable :: pmf(:)
    character(len=:), allocatable :: m
    integer :: n
    n = size(probs)
    if (n == 0) then
      allocate(pmf(0:0)); pmf(0)=1.0_dp; return
    end if
    if (maxval(abs(probs-probs(1))) <= 16.0_dp*epsilon(1.0_dp)) then
      pmf = binomial_approx(n,probs(1))
      return
    end if
    m = trim(adjustl(method))
    select case (m)
    case ("DivideFFT","dividefft","DIVIDEFFT")
      pmf = dpb_dividefft(probs)
    case ("Convolve","convolve","CONVOLVE")
      pmf = dpb_convolve(probs)
    case ("Characteristic","characteristic","CHARACTERISTIC")
      pmf = dpb_characteristic(probs)
    case ("Recursive","recursive","RECURSIVE")
      pmf = dpb_recursive(probs)
    case ("Mean","mean","MEAN")
      pmf = dpb_mean(probs)
    case ("GeoMean","geomean","GEOMEAN")
      pmf = dpb_geomean(probs)
    case ("GeoMeanCounter","geomeancounter","GEOMEANCOUNTER")
      pmf = dpb_geomean_counter(probs)
    case ("Poisson","poisson","POISSON")
      pmf = dpb_poisson(probs)
    case ("Normal","normal","NORMAL")
      pmf = dpb_normal(probs,.false.)
    case ("RefinedNormal","refinednormal","REFINEDNORMAL")
      pmf = dpb_normal(probs,.true.)
    case default
      error stop "unknown Poisson-binomial method"
    end select
  end function reduced_method_pmf

  function dpbinom(probs, method, wts) result(pmf)
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    real(dp), allocatable :: pmf(:), pexp(:), pred(:), pr(:)
    character(len=32) :: m
    integer :: n, n1, i, j
    m = "DivideFFT"; if (present(method)) m = method
    call expand_probs(probs,wts,pexp)
    if (any(pexp < 0.0_dp) .or. any(pexp > 1.0_dp)) then
      error stop "probabilities must lie in [0,1]"
    end if
    n = size(pexp)
    n1 = count(pexp >= 1.0_dp)
    allocate(pred(count(pexp > 0.0_dp .and. pexp < 1.0_dp)))
    j = 0
    do i = 1, n
      if (pexp(i) > 0.0_dp .and. pexp(i) < 1.0_dp) then
        j = j + 1
        pred(j) = pexp(i)
      end if
    end do
    pr = reduced_method_pmf(pred,m)
    allocate(pmf(0:n)); pmf = 0.0_dp
    do i = 1, size(pr)
      pmf(n1+i-1) = pr(i)
    end do
  end function dpbinom

  real(dp) function dpbinom_at(x, probs, method, wts, log_p) result(p)
    integer, intent(in) :: x
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: log_p
    real(dp), allocatable :: pmf(:)
    logical :: lp
    lp=.false.; if(present(log_p)) lp=log_p
    pmf = dpbinom(probs,method,wts)
    if (x < 0 .or. x > size(pmf)-1) then
      p = 0.0_dp
    else
      p = pmf(x+1)
    end if
    if (lp) p = log(p)
  end function dpbinom_at

  function ppbinom(probs, method, wts, lower_tail) result(cdf)
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: lower_tail
    real(dp), allocatable :: cdf(:), pmf(:)
    logical :: lower
    integer :: k, n
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    pmf=dpbinom(probs,method,wts)
    n=size(pmf)-1
    allocate(cdf(0:n))
    if(lower) then
      cdf(0)=pmf(1)
      do k=1,n; cdf(k)=cdf(k-1)+pmf(k+1); end do
      cdf(n)=1.0_dp
    else
      cdf(n)=0.0_dp
      do k=n-1,0,-1; cdf(k)=cdf(k+1)+pmf(k+2); end do
    end if
  end function ppbinom

  real(dp) function ppbinom_at(x, probs, method, wts, lower_tail, log_p) result(p)
    integer, intent(in) :: x
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: lower_tail, log_p
    real(dp), allocatable :: cdf(:)
    logical :: lower, lp
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    cdf=ppbinom(probs,method,wts,lower)
    if(x<0) then
      p=merge(0.0_dp,1.0_dp,lower)
    else if(x>size(cdf)-1) then
      p=merge(1.0_dp,0.0_dp,lower)
    else
      p=cdf(x+1)
    end if
    if(lp) p=log(p)
  end function ppbinom_at

  integer function qpbinom(prob, probs, method, wts, lower_tail, log_p) result(q)
    real(dp), intent(in) :: prob
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: lower_tail, log_p
    real(dp), allocatable :: cdf(:), pmf(:)
    logical :: lower, lp
    real(dp) :: target
    integer :: k, lo, hi
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    target=prob; if(lp) target=exp(prob)
    if(target<0.0_dp .or. target>1.0_dp) error stop "quantile probability outside [0,1]"
    pmf=dpbinom(probs,method,wts)
    lo=0; hi=size(pmf)-1
    do while (lo < hi)
      if (abs(pmf(lo+1)) > tiny(1.0_dp)) exit
      lo = lo + 1
    end do
    do while (hi > lo)
      if (abs(pmf(hi+1)) > tiny(1.0_dp)) exit
      hi = hi - 1
    end do
    if(lower) then
      if(target<=0.0_dp) then; q=lo; return; end if
      if(target>=1.0_dp) then; q=hi; return; end if
      cdf=ppbinom(probs,method,wts,.true.)
      q=hi
      do k=lo,hi
        if(cdf(k+1)>=target) then; q=k; exit; end if
      end do
    else
      if(target>=1.0_dp) then; q=lo; return; end if
      if(target<=0.0_dp) then; q=hi; return; end if
      cdf=ppbinom(probs,method,wts,.false.)
      q=hi
      do k=lo,hi
        if(cdf(k+1)<=target) then; q=k; exit; end if
      end do
    end if
  end function qpbinom

  function dpbinom_values(x, probs, method, wts, log_p) result(v)
    integer, intent(in) :: x(:)
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: log_p
    real(dp), allocatable :: v(:), pmf(:)
    logical :: lp
    integer :: i, n
    lp=.false.; if(present(log_p)) lp=log_p
    pmf=dpbinom(probs,method,wts)
    n=size(pmf)-1
    allocate(v(size(x)))
    do i=1,size(x)
      if(x(i)<0 .or. x(i)>n) then
        v(i)=0.0_dp
      else
        v(i)=pmf(x(i)+1)
      end if
    end do
    if(lp) v=log(v)
  end function dpbinom_values

  function ppbinom_values(x, probs, method, wts, lower_tail, log_p) result(v)
    integer, intent(in) :: x(:)
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: lower_tail, log_p
    real(dp), allocatable :: v(:), cdf(:)
    logical :: lower, lp
    integer :: i, n
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    cdf=ppbinom(probs,method,wts,lower)
    n=size(cdf)-1
    allocate(v(size(x)))
    do i=1,size(x)
      if(x(i)<0) then
        v(i)=merge(0.0_dp,1.0_dp,lower)
      else if(x(i)>n) then
        v(i)=merge(1.0_dp,0.0_dp,lower)
      else
        v(i)=cdf(x(i)+1)
      end if
    end do
    if(lp) v=log(v)
  end function ppbinom_values

  function qpbinom_values(prob, probs, method, wts, lower_tail, log_p) result(q)
    real(dp), intent(in) :: prob(:)
    real(dp), intent(in) :: probs(:)
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: lower_tail, log_p
    integer, allocatable :: q(:)
    integer :: i
    allocate(q(size(prob)))
    do i=1,size(prob)
      q(i)=qpbinom(prob(i),probs,method,wts,lower_tail,log_p)
    end do
  end function qpbinom_values

  subroutine rpbinom(n, probs, draws, method, wts, generator)
    integer, intent(in) :: n
    real(dp), intent(in) :: probs(:)
    integer, allocatable, intent(out) :: draws(:)
    character(len=*), intent(in), optional :: method, generator
    integer, intent(in), optional :: wts(:)
    real(dp), allocatable :: pexp(:), pmf(:)
    character(len=32) :: m, g
    integer :: i,j,s
    real(dp) :: u
    if(n<0) error stop "n must be nonnegative"
    m="DivideFFT"; if(present(method)) m=method
    g="Sample"; if(present(generator)) g=generator
    call expand_probs(probs,wts,pexp)
    if(any(pexp<0.0_dp).or.any(pexp>1.0_dp)) error stop "probabilities must lie in [0,1]"
    allocate(draws(n))
    select case(trim(adjustl(g)))
    case("Bernoulli","bernoulli","BERNOULLI")
      do j=1,n
        s=0
        do i=1,size(pexp)
          call random_number(u)
          if(u<pexp(i)) s=s+1
        end do
        draws(j)=s
      end do
    case("Sample","sample","SAMPLE")
      pmf=dpbinom(pexp,m)
      do j=1,n; draws(j)=sample_from_pmf(pmf,0); end do
    case default
      error stop "unknown generator"
    end select
  end subroutine rpbinom

end module pb_ordinary
