! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_distributions
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mixtools_kinds, only : dp, pi
  use mixtools_status, only : MIXTOOLS_SUCCESS, MIXTOOLS_DIMENSION_ERROR, MIXTOOLS_INVALID_ARGUMENT
  use mixtools_linalg, only : cholesky_lower, solve_lower, inverse_spd
  use mixtools_rng, only : rng_state, random_normal, random_gamma, random_poisson
  use mixtools_rng, only : random_weibull, random_uniform
  use r_stability, only : r_log_sum_exp
  implicit none
  private
  public :: normal_pdf, normal_logpdf, gamma_logpdf, poisson_logpmf
  public :: binomial_logpmf, multinomial_logpmf, ddirichlet
  public :: dmvnorm, logdmvnorm, logistic, logsumexp, normalize_logweights
  public :: digamma_approx, trigamma_approx, rexpmix, dexpmixt
  public :: rnormmix, rmvnorm, rmvnormmix, rweibullmix, rlnormscalemix
contains
  elemental function normal_logpdf(x, mu, sigma) result(v)
    real(dp), intent(in) :: x, mu, sigma
    real(dp) :: v
    if (sigma <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      v = -0.5_dp*log(2.0_dp*pi) - log(sigma) - 0.5_dp*((x-mu)/sigma)**2
    end if
  end function normal_logpdf

  elemental function normal_pdf(x, mu, sigma) result(v)
    real(dp), intent(in) :: x, mu, sigma
    real(dp) :: v
    v = exp(normal_logpdf(x,mu,sigma))
  end function normal_pdf

  elemental function gamma_logpdf(x, shape, scale) result(v)
    real(dp), intent(in) :: x, shape, scale
    real(dp) :: v
    if (x <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      v = (shape-1.0_dp)*log(x) - x/scale - log_gamma(shape) - shape*log(scale)
    end if
  end function gamma_logpdf

  elemental function poisson_logpmf(y, mean) result(v)
    real(dp), intent(in) :: y, mean
    real(dp) :: v
    if (y < 0.0_dp .or. mean <= 0.0_dp) then
      v = -huge(1.0_dp)
    else
      v = y*log(mean) - mean - log_gamma(y+1.0_dp)
    end if
  end function poisson_logpmf

  elemental function binomial_logpmf(y, n, prob) result(v)
    real(dp), intent(in) :: y, n, prob
    real(dp) :: v, p
    p = min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),prob))
    if (y < 0.0_dp .or. y > n .or. n < 0.0_dp) then
      v = -huge(1.0_dp)
    else
      v = log_gamma(n+1.0_dp)-log_gamma(y+1.0_dp)-log_gamma(n-y+1.0_dp) &
        + y*log(p)+(n-y)*log(1.0_dp-p)
    end if
  end function binomial_logpmf

  function multinomial_logpmf(y, theta) result(v)
    real(dp), intent(in) :: y(:), theta(:)
    real(dp) :: v
    integer :: j
    if (size(y) /= size(theta) .or. any(y < 0.0_dp) .or. any(theta <= 0.0_dp)) then
      v = -huge(1.0_dp); return
    end if
    v = log_gamma(sum(y)+1.0_dp)
    do j = 1, size(y)
      v = v - log_gamma(y(j)+1.0_dp) + y(j)*log(theta(j))
    end do
  end function multinomial_logpmf

  function ddirichlet(x, alpha) result(v)
    real(dp), intent(in) :: x(:), alpha(:)
    real(dp) :: v
    if (size(x) /= size(alpha) .or. any(x <= 0.0_dp) .or. any(alpha <= 0.0_dp)) then
      v = 0.0_dp; return
    end if
    if (abs(sum(x)-1.0_dp) > 1.0e-8_dp) then
      v = 0.0_dp; return
    end if
    v = exp(log_gamma(sum(alpha)) - sum(log_gamma(alpha)) + sum((alpha-1.0_dp)*log(x)))
  end function ddirichlet

  function logdmvnorm(y, mu, sigma, status) result(v)
    real(dp), intent(in) :: y(:), mu(:), sigma(:,:)
    integer, intent(out), optional :: status
    real(dp) :: v, logdet
    real(dp), allocatable :: inv(:,:), d(:)
    integer :: stat, p
    p = size(y)
    if (size(mu) /= p .or. size(sigma,1) /= p .or. size(sigma,2) /= p) then
      v = -huge(1.0_dp); if (present(status)) status=MIXTOOLS_DIMENSION_ERROR; return
    end if
    call inverse_spd(sigma, inv, stat, logdet)
    if (stat /= MIXTOOLS_SUCCESS) then
      v = -huge(1.0_dp); if (present(status)) status=stat; return
    end if
    allocate(d(p)); d = y-mu
    v = -0.5_dp*(real(p,dp)*log(2.0_dp*pi)+logdet+dot_product(d,matmul(inv,d)))
    if (present(status)) status=MIXTOOLS_SUCCESS
  end function logdmvnorm

  function dmvnorm(y, mu, sigma, status) result(v)
    real(dp), intent(in) :: y(:), mu(:), sigma(:,:)
    integer, intent(out), optional :: status
    real(dp) :: v
    integer :: stat
    v = exp(logdmvnorm(y,mu,sigma,stat))
    if (present(status)) status=stat
  end function dmvnorm

  elemental function logistic(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    if (x >= 0.0_dp) then
      p = 1.0_dp/(1.0_dp+exp(-x))
    else
      p = exp(x)/(1.0_dp+exp(x))
    end if
  end function logistic

  function logsumexp(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v = r_log_sum_exp(x)
  end function logsumexp

  subroutine normalize_logweights(logw, prob, lognorm)
    real(dp), intent(in) :: logw(:)
    real(dp), intent(out) :: prob(size(logw))
    real(dp), intent(out) :: lognorm
    lognorm = logsumexp(logw)
    if (ieee_is_finite(lognorm)) then
      prob = exp(logw-lognorm)
    else
      prob = 1.0_dp/real(size(logw),dp)
    end if
  end subroutine normalize_logweights

  elemental function digamma_approx(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, y, inv, inv2
    y = x; v = 0.0_dp
    if (y <= 0.0_dp) then
      v = -huge(1.0_dp); return
    end if
    do while (y < 8.0_dp)
      v = v - 1.0_dp/y; y = y + 1.0_dp
    end do
    inv = 1.0_dp/y; inv2 = inv*inv
    v = v + log(y) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp &
      - inv2*(1.0_dp/252.0_dp)))
  end function digamma_approx

  elemental function trigamma_approx(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, y, inv, inv2
    y=x; v=0.0_dp
    if (y <= 0.0_dp) then
      v=huge(1.0_dp); return
    end if
    do while (y < 8.0_dp)
      v = v + 1.0_dp/(y*y); y=y+1.0_dp
    end do
    inv=1.0_dp/y; inv2=inv*inv
    v = v + inv + 0.5_dp*inv2 + inv*inv2/6.0_dp - inv*inv2*inv2/30.0_dp &
      + inv*inv2*inv2*inv2/42.0_dp
  end function trigamma_approx

  subroutine rexpmix(rng, n, lambda, rate, x, component)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: lambda(:), rate(:)
    real(dp), intent(out) :: x(n)
    integer, intent(out), optional :: component(n)
    integer :: i,j,c
    real(dp) :: u,s
    do i=1,n
      u=random_uniform(rng); s=0.0_dp; c=size(lambda)
      do j=1,size(lambda)
        s=s+lambda(j)
        if (u<=s) then; c=j; exit; end if
      end do
      x(i)=-log(random_uniform(rng))/rate(c)
      if (present(component)) component(i)=c
    end do
  end subroutine rexpmix

  function dexpmixt(t, lambda, rate) result(v)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: lambda(:), rate(:)
    real(dp) :: v
    integer :: j
    v=0.0_dp
    if (t<0.0_dp) return
    do j=1,size(lambda)
      v=v+lambda(j)*rate(j)*exp(-rate(j)*t)
    end do
  end function dexpmixt

  subroutine rnormmix(rng,n,lambda,mu,sigma,x,component)
    type(rng_state), intent(inout) :: rng
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),mu(:),sigma(:)
    real(dp),intent(out)::x(n)
    integer,intent(out),optional::component(n)
    integer::i,j,c
    real(dp)::u,s
    do i=1,n
      u=random_uniform(rng);s=0.0_dp;c=size(lambda)
      do j=1,size(lambda);s=s+lambda(j);if(u<=s)then;c=j;exit;end if;end do
      x(i)=mu(c)+sigma(c)*random_normal(rng)
      if(present(component))component(i)=c
    end do
  end subroutine rnormmix

  subroutine rmvnorm(rng,n,mu,sigma,x,status)
    type(rng_state),intent(inout)::rng
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),sigma(:,:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer,intent(out)::status
    real(dp),allocatable::l(:,:),z(:)
    integer::i,j,p
    p=size(mu);call cholesky_lower(sigma,l,status)
    if(status/=MIXTOOLS_SUCCESS)then;allocate(x(0,0));return;end if
    allocate(x(n,p),z(p))
    do i=1,n
      do j=1,p;z(j)=random_normal(rng);end do
      x(i,:)=mu+matmul(l,z)
    end do
  end subroutine rmvnorm

  subroutine rmvnormmix(rng,n,lambda,mu,sigma,x,component,status)
    type(rng_state),intent(inout)::rng
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),mu(:,:),sigma(:,:,:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer,allocatable,intent(out),optional::component(:)
    integer,intent(out)::status
    integer::i,j,c,k,p
    real(dp)::u,s
    real(dp),allocatable::one(:,:)
    k=size(lambda);p=size(mu,1)
    if(size(mu,2)/=k.or.size(sigma,1)/=p.or.size(sigma,2)/=p.or.size(sigma,3)/=k)then
      allocate(x(0,0));status=MIXTOOLS_DIMENSION_ERROR;return
    end if
    allocate(x(n,p));if(present(component))allocate(component(n))
    do i=1,n
      u=random_uniform(rng);s=0.0_dp;c=k
      do j=1,k;s=s+lambda(j);if(u<=s)then;c=j;exit;end if;end do
      call rmvnorm(rng,1,mu(:,c),sigma(:,:,c),one,status)
      if(status/=MIXTOOLS_SUCCESS)return
      x(i,:)=one(1,:);if(present(component))component(i)=c
    end do
  end subroutine rmvnormmix

  subroutine rweibullmix(rng,n,lambda,shape,scale,x,component)
    type(rng_state),intent(inout)::rng
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),shape(:),scale(:)
    real(dp),intent(out)::x(n)
    integer,intent(out),optional::component(n)
    integer::i,j,c
    real(dp)::u,s
    do i=1,n
      u=random_uniform(rng);s=0.0_dp;c=size(lambda)
      do j=1,size(lambda);s=s+lambda(j);if(u<=s)then;c=j;exit;end if;end do
      x(i)=random_weibull(rng,shape(c),scale(c));if(present(component))component(i)=c
    end do
  end subroutine rweibullmix

  subroutine rlnormscalemix(rng,n,lambda,meanlog,sdlog,scale,x,component)
    type(rng_state),intent(inout)::rng
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),meanlog(:),sdlog(:),scale(:)
    real(dp),intent(out)::x(n)
    integer,intent(out),optional::component(n)
    integer::i,j,c
    real(dp)::u,s
    do i=1,n
      u=random_uniform(rng);s=0.0_dp;c=size(lambda)
      do j=1,size(lambda);s=s+lambda(j);if(u<=s)then;c=j;exit;end if;end do
      x(i)=scale(c)*exp(meanlog(c)+sdlog(c)*random_normal(rng))
      if(present(component))component(i)=c
    end do
  end subroutine rlnormscalemix
end module mixtools_distributions
