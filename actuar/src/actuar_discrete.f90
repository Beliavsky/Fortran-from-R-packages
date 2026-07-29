! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_discrete
  use actuar_kinds, only : dp, huge_dp
  use actuar_special, only : nan_dp
  use actuar_rng, only : runif, rpois, rbinom, rnbinom, rinvgauss_rng
  implicit none
  private
  public :: dlogarithmic, plogarithmic, qlogarithmic, rlogarithmic
  public :: dztpois, pztpois, qztpois, rztpois
  public :: dztgeom, pztgeom, qztgeom, rztgeom
  public :: dztbinom, pztbinom, qztbinom, rztbinom
  public :: dztnbinom, pztnbinom, qztnbinom, rztnbinom
  public :: dzmpois, pzmpois, qzmpois, rzmpois
  public :: dzmgeom, pzmgeom, qzmgeom, rzmgeom
  public :: dzmbinom, pzmbinom, qzmbinom, rzmbinom
  public :: dzmnbinom, pzmnbinom, qzmnbinom, rzmnbinom
  public :: dzmlogarithmic, pzmlogarithmic, qzmlogarithmic, rzmlogarithmic
  public :: dpoisinvgauss, ppoisinvgauss, qpoisinvgauss, rpoisinvgauss
  public :: poisson_pmf, binomial_pmf, nbinomial_pmf
contains

  pure function is_integer_value(x) result(ok)
    real(dp), intent(in) :: x
    logical :: ok
    ok = abs(x-anint(x)) <= 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x))
  end function is_integer_value

  pure function poisson_pmf(k, lambda) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    real(dp) :: p
    if (k < 0 .or. lambda < 0.0_dp) then
      p = 0.0_dp
    else if (lambda == 0.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==0)
    else
      p = exp(-lambda+real(k,dp)*log(lambda)-log_gamma(real(k+1,dp)))
    end if
  end function poisson_pmf

  pure function poisson_cdf(k, lambda) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    real(dp) :: p, term
    integer :: i
    if (k < 0) then
      p = 0.0_dp; return
    end if
    term = exp(-lambda); p = term
    do i = 1, k
      term = term*lambda/real(i,dp)
      p = p+term
    end do
    p = min(1.0_dp,p)
  end function poisson_cdf

  function poisson_quantile(prob, lambda) result(k)
    real(dp), intent(in) :: prob, lambda
    integer :: k
    real(dp) :: cdf, term
    if (prob <= 0.0_dp) then
      k = 0; return
    else if (prob >= 1.0_dp) then
      k = huge(0); return
    end if
    k = 0; term = exp(-lambda); cdf = term
    do while (cdf < prob .and. k < 10000000)
      k = k+1
      term = term*lambda/real(k,dp)
      cdf = cdf+term
    end do
  end function poisson_quantile

  pure function binomial_pmf(k, n, prob) result(p)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: prob
    real(dp) :: p
    if (k < 0 .or. k > n .or. n < 0 .or. prob < 0.0_dp .or. prob > 1.0_dp) then
      p = 0.0_dp
    else if (prob == 0.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==0)
    else if (prob == 1.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==n)
    else
      p = exp(log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))- &
          log_gamma(real(n-k+1,dp))+real(k,dp)*log(prob)+ &
          real(n-k,dp)*log(1.0_dp-prob))
    end if
  end function binomial_pmf

  pure function binomial_cdf(k, n, prob) result(p)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: prob
    real(dp) :: p
    integer :: i
    p = 0.0_dp
    do i = 0, min(k,n)
      p = p+binomial_pmf(i,n,prob)
    end do
    p = min(1.0_dp,p)
  end function binomial_cdf

  function binomial_quantile(q, n, prob) result(k)
    real(dp), intent(in) :: q, prob
    integer, intent(in) :: n
    integer :: k
    if (q <= 0.0_dp) then
      k = 0; return
    end if
    do k = 0, n
      if (binomial_cdf(k,n,prob) >= q) return
    end do
    k = n
  end function binomial_quantile

  pure function nbinomial_pmf(k, size, prob) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: size, prob
    real(dp) :: p
    if (k < 0 .or. size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) then
      p = 0.0_dp
    else if (prob == 1.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==0)
    else
      p = exp(log_gamma(real(k,dp)+size)-log_gamma(size)- &
          log_gamma(real(k+1,dp))+size*log(prob)+real(k,dp)*log(1.0_dp-prob))
    end if
  end function nbinomial_pmf

  pure function nbinomial_cdf(k, size, prob) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: size, prob
    real(dp) :: p
    integer :: i
    p = 0.0_dp
    do i = 0, k
      p = p+nbinomial_pmf(i,size,prob)
    end do
    p = min(1.0_dp,p)
  end function nbinomial_cdf

  function nbinomial_quantile(q, size, prob) result(k)
    real(dp), intent(in) :: q, size, prob
    integer :: k
    if (q <= 0.0_dp) then
      k = 0; return
    end if
    do k = 0, 10000000
      if (nbinomial_cdf(k,size,prob) >= q) return
    end do
    k = huge(0)
  end function nbinomial_quantile

  pure function dlogarithmic(x, prob) result(p)
    real(dp), intent(in) :: x, prob
    real(dp) :: p
    integer :: k
    if (prob < 0.0_dp .or. prob >= 1.0_dp) then
      p = nan_dp(); return
    end if
    if (x < 1.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    if (prob == 0.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==1)
    else
      p = -prob**k/(real(k,dp)*log(1.0_dp-prob))
    end if
  end function dlogarithmic

  pure function plogarithmic(x, prob) result(p)
    real(dp), intent(in) :: x, prob
    real(dp) :: p, term
    integer :: k, upper
    if (prob < 0.0_dp .or. prob >= 1.0_dp) then
      p = nan_dp(); return
    end if
    if (x < 1.0_dp) then
      p = 0.0_dp; return
    end if
    upper = floor(x)
    if (prob == 0.0_dp) then
      p = 1.0_dp; return
    end if
    term = -prob/log(1.0_dp-prob); p = term
    do k = 1, upper-1
      term = term*prob*real(k,dp)/real(k+1,dp)
      p = p+term
    end do
    p = min(1.0_dp,p)
  end function plogarithmic

  function qlogarithmic(q, prob) result(k)
    real(dp), intent(in) :: q, prob
    integer :: k
    real(dp) :: p, term
    if (q <= 0.0_dp) then
      k = 1; return
    else if (q >= 1.0_dp) then
      k = huge(0); return
    end if
    if (prob == 0.0_dp) then
      k = 1; return
    end if
    k = 1; term = -prob/log(1.0_dp-prob); p = term
    do while (p < q .and. k < 10000000)
      term = term*prob*real(k,dp)/real(k+1,dp)
      k = k+1
      p = p+term
    end do
  end function qlogarithmic

  function rlogarithmic(prob) result(k)
    real(dp), intent(in) :: prob
    integer :: k
    k = qlogarithmic(runif(),prob)
  end function rlogarithmic

  pure function dztpois(x, lambda) result(p)
    real(dp), intent(in) :: x, lambda
    real(dp) :: p, den
    integer :: k
    if (lambda < 0.0_dp) then
      p = nan_dp(); return
    end if
    if (x < 1.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    if (lambda == 0.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==1)
    else
      den = 1.0_dp-exp(-lambda)
      p = poisson_pmf(k,lambda)/den
    end if
  end function dztpois

  pure function pztpois(x, lambda) result(p)
    real(dp), intent(in) :: x, lambda
    real(dp) :: p, p0
    if (x < 1.0_dp) then
      p = 0.0_dp
    else if (lambda == 0.0_dp) then
      p = 1.0_dp
    else
      p0 = exp(-lambda)
      p = (poisson_cdf(floor(x),lambda)-p0)/(1.0_dp-p0)
    end if
  end function pztpois

  function qztpois(q, lambda) result(k)
    real(dp), intent(in) :: q, lambda
    integer :: k
    if (lambda == 0.0_dp) then
      k = 1
    else
      k = poisson_quantile(exp(-lambda)+(1.0_dp-exp(-lambda))*q,lambda)
      k = max(1,k)
    end if
  end function qztpois

  function rztpois(lambda) result(k)
    real(dp), intent(in) :: lambda
    integer :: k
    k = qztpois(runif(),lambda)
  end function rztpois

  pure function dztgeom(x, prob) result(p)
    real(dp), intent(in) :: x, prob
    real(dp) :: p
    integer :: k
    if (prob <= 0.0_dp .or. prob > 1.0_dp) then
      p = nan_dp(); return
    end if
    if (x < 1.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    p = prob*(1.0_dp-prob)**(k-1)
  end function dztgeom

  pure function pztgeom(x, prob) result(p)
    real(dp), intent(in) :: x, prob
    real(dp) :: p
    if (x < 1.0_dp) then
      p = 0.0_dp
    else
      p = 1.0_dp-(1.0_dp-prob)**floor(x)
    end if
  end function pztgeom

  pure function qztgeom(q, prob) result(k)
    real(dp), intent(in) :: q, prob
    integer :: k
    if (q <= 0.0_dp .or. prob == 1.0_dp) then
      k = 1
    else if (q >= 1.0_dp) then
      k = huge(0)
    else
      k = ceiling(log(1.0_dp-q)/log(1.0_dp-prob))
      k = max(1,k)
    end if
  end function qztgeom

  function rztgeom(prob) result(k)
    real(dp), intent(in) :: prob
    integer :: k
    k = qztgeom(runif(),prob)
  end function rztgeom

  pure function dztbinom(x, n, prob) result(p)
    real(dp), intent(in) :: x, prob
    integer, intent(in) :: n
    real(dp) :: p, p0
    integer :: k
    if (x < 1.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x); p0 = (1.0_dp-prob)**n
    if (p0 >= 1.0_dp) then
      p = merge(1.0_dp,0.0_dp,k==1)
    else
      p = binomial_pmf(k,n,prob)/(1.0_dp-p0)
    end if
  end function dztbinom

  pure function pztbinom(x, n, prob) result(p)
    real(dp), intent(in) :: x, prob
    integer, intent(in) :: n
    real(dp) :: p, p0
    if (x < 1.0_dp) then
      p = 0.0_dp; return
    end if
    p0 = (1.0_dp-prob)**n
    if (p0 >= 1.0_dp) then
      p = 1.0_dp
    else
      p = (binomial_cdf(floor(x),n,prob)-p0)/(1.0_dp-p0)
    end if
  end function pztbinom

  function qztbinom(q, n, prob) result(k)
    real(dp), intent(in) :: q, prob
    integer, intent(in) :: n
    integer :: k
    real(dp) :: p0
    p0 = (1.0_dp-prob)**n
    if (p0 >= 1.0_dp) then
      k = 1
    else
      k = max(1,binomial_quantile(p0+(1.0_dp-p0)*q,n,prob))
    end if
  end function qztbinom

  function rztbinom(n, prob) result(k)
    integer, intent(in) :: n
    real(dp), intent(in) :: prob
    integer :: k
    k = qztbinom(runif(),n,prob)
  end function rztbinom

  pure function dztnbinom(x, size, prob) result(p)
    real(dp), intent(in) :: x, size, prob
    real(dp) :: p, p0
    integer :: k
    if (x < 1.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x); p0 = prob**size
    p = nbinomial_pmf(k,size,prob)/(1.0_dp-p0)
  end function dztnbinom

  pure function pztnbinom(x, size, prob) result(p)
    real(dp), intent(in) :: x, size, prob
    real(dp) :: p, p0
    if (x < 1.0_dp) then
      p = 0.0_dp; return
    end if
    p0 = prob**size
    p = (nbinomial_cdf(floor(x),size,prob)-p0)/(1.0_dp-p0)
  end function pztnbinom

  function qztnbinom(q, size, prob) result(k)
    real(dp), intent(in) :: q, size, prob
    integer :: k
    real(dp) :: p0
    p0 = prob**size
    k = max(1,nbinomial_quantile(p0+(1.0_dp-p0)*q,size,prob))
  end function qztnbinom

  function rztnbinom(size, prob) result(k)
    real(dp), intent(in) :: size, prob
    integer :: k
    k = qztnbinom(runif(),size,prob)
  end function rztnbinom

  pure function zero_modified_pmf(base, p0, p0m, k) result(p)
    real(dp), intent(in) :: base, p0, p0m
    integer, intent(in) :: k
    real(dp) :: p
    if (k == 0) then
      p = p0m
    else
      p = (1.0_dp-p0m)*base/(1.0_dp-p0)
    end if
  end function zero_modified_pmf

  pure function dzmpois(x, lambda, p0m) result(p)
    real(dp), intent(in) :: x, lambda, p0m
    real(dp) :: p
    integer :: k
    if (x < 0.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    if (lambda == 0.0_dp) then
      if (k == 0) then; p = p0m
      else if (k == 1) then; p = 1.0_dp-p0m
      else; p = 0.0_dp
      end if
    else
      p = zero_modified_pmf(poisson_pmf(k,lambda),exp(-lambda),p0m,k)
    end if
  end function dzmpois

  pure function pzmpois(x, lambda, p0m) result(p)
    real(dp), intent(in) :: x, lambda, p0m
    real(dp) :: p, p0
    if (x < 0.0_dp) then
      p = 0.0_dp; return
    else if (x < 1.0_dp) then
      p = p0m; return
    end if
    if (lambda == 0.0_dp) then
      p = 1.0_dp
    else
      p0 = exp(-lambda)
      p = p0m+(1.0_dp-p0m)*(poisson_cdf(floor(x),lambda)-p0)/(1.0_dp-p0)
    end if
  end function pzmpois

  function qzmpois(q, lambda, p0m) result(k)
    real(dp), intent(in) :: q, lambda, p0m
    integer :: k
    if (q <= p0m) then
      k = 0
    else
      k = qztpois((q-p0m)/(1.0_dp-p0m),lambda)
    end if
  end function qzmpois

  function rzmpois(lambda, p0m) result(k)
    real(dp), intent(in) :: lambda, p0m
    integer :: k
    k = qzmpois(runif(),lambda,p0m)
  end function rzmpois

  pure function dzmgeom(x, prob, p0m) result(p)
    real(dp), intent(in) :: x, prob, p0m
    real(dp) :: p
    integer :: k
    if (x < 0.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    if (k == 0) then
      p = p0m
    else
      p = (1.0_dp-p0m)*prob*(1.0_dp-prob)**(k-1)
    end if
  end function dzmgeom

  pure function pzmgeom(x, prob, p0m) result(p)
    real(dp), intent(in) :: x, prob, p0m
    real(dp) :: p
    if (x < 0.0_dp) then
      p = 0.0_dp
    else if (x < 1.0_dp) then
      p = p0m
    else
      p = p0m+(1.0_dp-p0m)*pztgeom(x,prob)
    end if
  end function pzmgeom

  pure function qzmgeom(q, prob, p0m) result(k)
    real(dp), intent(in) :: q, prob, p0m
    integer :: k
    if (q <= p0m) then
      k = 0
    else
      k = qztgeom((q-p0m)/(1.0_dp-p0m),prob)
    end if
  end function qzmgeom

  function rzmgeom(prob, p0m) result(k)
    real(dp), intent(in) :: prob, p0m
    integer :: k
    k = qzmgeom(runif(),prob,p0m)
  end function rzmgeom

  pure function dzmbinom(x, n, prob, p0m) result(p)
    real(dp), intent(in) :: x, prob, p0m
    integer, intent(in) :: n
    real(dp) :: p
    integer :: k
    if (x < 0.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    p = zero_modified_pmf(binomial_pmf(k,n,prob),(1.0_dp-prob)**n,p0m,k)
  end function dzmbinom

  pure function pzmbinom(x, n, prob, p0m) result(p)
    real(dp), intent(in) :: x, prob, p0m
    integer, intent(in) :: n
    real(dp) :: p, p0
    if (x < 0.0_dp) then
      p = 0.0_dp; return
    else if (x < 1.0_dp) then
      p = p0m; return
    end if
    p0 = (1.0_dp-prob)**n
    p = p0m+(1.0_dp-p0m)*(binomial_cdf(floor(x),n,prob)-p0)/(1.0_dp-p0)
  end function pzmbinom

  function qzmbinom(q, n, prob, p0m) result(k)
    real(dp), intent(in) :: q, prob, p0m
    integer, intent(in) :: n
    integer :: k
    if (q <= p0m) then
      k = 0
    else
      k = qztbinom((q-p0m)/(1.0_dp-p0m),n,prob)
    end if
  end function qzmbinom

  function rzmbinom(n, prob, p0m) result(k)
    integer, intent(in) :: n
    real(dp), intent(in) :: prob, p0m
    integer :: k
    k = qzmbinom(runif(),n,prob,p0m)
  end function rzmbinom

  pure function dzmnbinom(x, size, prob, p0m) result(p)
    real(dp), intent(in) :: x, size, prob, p0m
    real(dp) :: p
    integer :: k
    if (x < 0.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    p = zero_modified_pmf(nbinomial_pmf(k,size,prob),prob**size,p0m,k)
  end function dzmnbinom

  pure function pzmnbinom(x, size, prob, p0m) result(p)
    real(dp), intent(in) :: x, size, prob, p0m
    real(dp) :: p, p0
    if (x < 0.0_dp) then
      p = 0.0_dp; return
    else if (x < 1.0_dp) then
      p = p0m; return
    end if
    p0 = prob**size
    p = p0m+(1.0_dp-p0m)*(nbinomial_cdf(floor(x),size,prob)-p0)/(1.0_dp-p0)
  end function pzmnbinom

  function qzmnbinom(q, size, prob, p0m) result(k)
    real(dp), intent(in) :: q, size, prob, p0m
    integer :: k
    if (q <= p0m) then
      k = 0
    else
      k = qztnbinom((q-p0m)/(1.0_dp-p0m),size,prob)
    end if
  end function qzmnbinom

  function rzmnbinom(size, prob, p0m) result(k)
    real(dp), intent(in) :: size, prob, p0m
    integer :: k
    k = qzmnbinom(runif(),size,prob,p0m)
  end function rzmnbinom

  pure function dzmlogarithmic(x, prob, p0m) result(p)
    real(dp), intent(in) :: x, prob, p0m
    real(dp) :: p
    integer :: k
    if (x < 0.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    if (k == 0) then
      p = p0m
    else
      p = (1.0_dp-p0m)*dlogarithmic(x,prob)
    end if
  end function dzmlogarithmic

  pure function pzmlogarithmic(x, prob, p0m) result(p)
    real(dp), intent(in) :: x, prob, p0m
    real(dp) :: p
    if (x < 0.0_dp) then
      p = 0.0_dp
    else if (x < 1.0_dp) then
      p = p0m
    else
      p = p0m+(1.0_dp-p0m)*plogarithmic(x,prob)
    end if
  end function pzmlogarithmic

  function qzmlogarithmic(q, prob, p0m) result(k)
    real(dp), intent(in) :: q, prob, p0m
    integer :: k
    if (q <= p0m) then
      k = 0
    else
      k = qlogarithmic((q-p0m)/(1.0_dp-p0m),prob)
    end if
  end function qzmlogarithmic

  function rzmlogarithmic(prob, p0m) result(k)
    real(dp), intent(in) :: prob, p0m
    integer :: k
    k = qzmlogarithmic(runif(),prob,p0m)
  end function rzmlogarithmic

  pure function dpoisinvgauss(x, mu, phi) result(p)
    real(dp), intent(in) :: x, mu, phi
    real(dp) :: p, pim2, pim1, a, b
    integer :: k, i
    if (mu <= 0.0_dp .or. phi <= 0.0_dp) then
      p = nan_dp(); return
    end if
    if (x < 0.0_dp .or. .not. is_integer_value(x)) then
      p = 0.0_dp; return
    end if
    k = nint(x)
    p = exp((1.0_dp-sqrt(1.0_dp+2.0_dp*phi*mu*mu))/(phi*mu))
    if (k == 0) return
    pim2 = p
    p = mu*p/sqrt(1.0_dp+2.0_dp*phi*mu*mu)
    if (k == 1) return
    pim1 = p
    a = 1.0_dp/(1.0_dp+1.0_dp/(2.0_dp*phi*mu*mu))
    b = mu*mu/(1.0_dp+2.0_dp*phi*mu*mu)
    do i = 2, k
      p = a*(1.0_dp-1.5_dp/real(i,dp))*pim1 + &
          b*pim2/(real(i,dp)*real(i-1,dp))
      pim2 = pim1; pim1 = p
    end do
  end function dpoisinvgauss

  pure function ppoisinvgauss(x, mu, phi) result(p)
    real(dp), intent(in) :: x, mu, phi
    real(dp) :: p
    integer :: i
    if (x < 0.0_dp) then
      p = 0.0_dp; return
    end if
    p = 0.0_dp
    do i = 0, floor(x)
      p = p+dpoisinvgauss(real(i,dp),mu,phi)
    end do
    p = min(1.0_dp,p)
  end function ppoisinvgauss

  function qpoisinvgauss(q, mu, phi) result(k)
    real(dp), intent(in) :: q, mu, phi
    integer :: k
    if (q <= 0.0_dp) then
      k = 0; return
    else if (q >= 1.0_dp) then
      k = huge(0); return
    end if
    do k = 0, 10000000
      if (ppoisinvgauss(real(k,dp),mu,phi) >= q) return
    end do
    k = huge(0)
  end function qpoisinvgauss

  function rpoisinvgauss(mu, phi) result(k)
    real(dp), intent(in) :: mu, phi
    integer :: k
    k = rpois(rinvgauss_rng(mu,phi))
  end function rpoisinvgauss

end module actuar_discrete
