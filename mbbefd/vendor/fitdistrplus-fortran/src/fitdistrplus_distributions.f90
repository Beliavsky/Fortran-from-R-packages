! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_distributions
  use fitdistrplus_kinds, only : dp
  use fitdistrplus_types, only : distribution_model
  use fitdistrplus_math, only : pi_dp, normal_cdf, normal_quantile, regularized_gamma_p, &
    regularized_beta, random_normal, random_gamma, random_poisson, clamp_probability
  implicit none
  private

  public :: make_normal, make_lognormal, make_exponential, make_gamma
  public :: make_weibull, make_uniform, make_logistic, make_cauchy
  public :: make_beta, make_poisson, make_geometric, make_negative_binomial
  public :: make_distribution

contains

  subroutine set_common(dist, name, npar, discrete, names, lower, upper)
    type(distribution_model), intent(out) :: dist
    character(len=*), intent(in) :: name
    integer, intent(in) :: npar
    logical, intent(in) :: discrete
    character(len=*), intent(in) :: names(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer :: i
    dist%name = name
    dist%npar = npar
    dist%discrete = discrete
    allocate(dist%parameter_names(npar), dist%default_lower(npar), dist%default_upper(npar))
    do i = 1, npar
      dist%parameter_names(i) = names(i)
    end do
    dist%default_lower = lower
    dist%default_upper = upper
  end subroutine set_common

  subroutine make_distribution(name, dist, status)
    character(len=*), intent(in) :: name
    type(distribution_model), intent(out) :: dist
    integer, intent(out), optional :: status
    character(len=32) :: key
    integer :: st
    key = lower_ascii(trim(name))
    st = 0
    select case (trim(key))
    case ("norm", "normal")
      call make_normal(dist)
    case ("lnorm", "lognormal")
      call make_lognormal(dist)
    case ("exp", "exponential")
      call make_exponential(dist)
    case ("gamma")
      call make_gamma(dist)
    case ("weibull")
      call make_weibull(dist)
    case ("unif", "uniform")
      call make_uniform(dist)
    case ("logis", "logistic")
      call make_logistic(dist)
    case ("cauchy")
      call make_cauchy(dist)
    case ("beta")
      call make_beta(dist)
    case ("pois", "poisson")
      call make_poisson(dist)
    case ("geom", "geometric")
      call make_geometric(dist)
    case ("nbinom", "negative-binomial", "negative_binomial")
      call make_negative_binomial(dist)
    case default
      st = 1
      dist%name = "unknown"
      dist%npar = 0
    end select
    if (present(status)) status = st
  end subroutine make_distribution

  subroutine make_normal(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"normal",2,.false.,[character(len=8)::"mean","sd"], &
      [-huge(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => normal_logpdf
    dist%cdf => normal_cdf_model
    dist%quantile => normal_quantile_model
    dist%raw_moment => normal_moment
    dist%random_value => normal_random
  end subroutine make_normal

  subroutine make_lognormal(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"lognormal",2,.false.,[character(len=8)::"meanlog","sdlog"], &
      [-huge(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => lognormal_logpdf
    dist%cdf => lognormal_cdf
    dist%quantile => lognormal_quantile
    dist%raw_moment => lognormal_moment
    dist%random_value => lognormal_random
  end subroutine make_lognormal

  subroutine make_exponential(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"exponential",1,.false.,[character(len=8)::"rate"], &
      [tiny(1.0_dp)],[huge(1.0_dp)])
    dist%logpdf => exponential_logpdf
    dist%cdf => exponential_cdf
    dist%quantile => exponential_quantile
    dist%raw_moment => exponential_moment
    dist%random_value => exponential_random
  end subroutine make_exponential

  subroutine make_gamma(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"gamma",2,.false.,[character(len=8)::"shape","rate"], &
      [tiny(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => gamma_logpdf
    dist%cdf => gamma_cdf
    dist%quantile => gamma_quantile
    dist%raw_moment => gamma_moment
    dist%random_value => gamma_random
  end subroutine make_gamma

  subroutine make_weibull(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"weibull",2,.false.,[character(len=8)::"shape","scale"], &
      [tiny(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => weibull_logpdf
    dist%cdf => weibull_cdf
    dist%quantile => weibull_quantile
    dist%raw_moment => weibull_moment
    dist%random_value => weibull_random
  end subroutine make_weibull

  subroutine make_uniform(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"uniform",2,.false.,[character(len=8)::"minimum","maximum"], &
      [-huge(1.0_dp),-huge(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => uniform_logpdf
    dist%cdf => uniform_cdf
    dist%quantile => uniform_quantile
    dist%raw_moment => uniform_moment
    dist%random_value => uniform_random
  end subroutine make_uniform

  subroutine make_logistic(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"logistic",2,.false.,[character(len=8)::"location","scale"], &
      [-huge(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => logistic_logpdf
    dist%cdf => logistic_cdf
    dist%quantile => logistic_quantile
    dist%raw_moment => logistic_moment
    dist%random_value => logistic_random
  end subroutine make_logistic

  subroutine make_cauchy(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"cauchy",2,.false.,[character(len=8)::"location","scale"], &
      [-huge(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => cauchy_logpdf
    dist%cdf => cauchy_cdf
    dist%quantile => cauchy_quantile
    dist%raw_moment => cauchy_moment
    dist%random_value => cauchy_random
  end subroutine make_cauchy

  subroutine make_beta(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"beta",2,.false.,[character(len=8)::"shape1","shape2"], &
      [tiny(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => beta_logpdf
    dist%cdf => beta_cdf
    dist%quantile => beta_quantile
    dist%raw_moment => beta_moment
    dist%random_value => beta_random
  end subroutine make_beta

  subroutine make_poisson(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"poisson",1,.true.,[character(len=8)::"lambda"], &
      [tiny(1.0_dp)],[huge(1.0_dp)])
    dist%logpdf => poisson_logpdf
    dist%cdf => poisson_cdf
    dist%quantile => poisson_quantile
    dist%raw_moment => poisson_moment
    dist%random_value => poisson_random
  end subroutine make_poisson

  subroutine make_geometric(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"geometric",1,.true.,[character(len=8)::"prob"], &
      [tiny(1.0_dp)],[1.0_dp])
    dist%logpdf => geometric_logpdf
    dist%cdf => geometric_cdf
    dist%quantile => geometric_quantile
    dist%raw_moment => geometric_moment
    dist%random_value => geometric_random
  end subroutine make_geometric

  subroutine make_negative_binomial(dist)
    type(distribution_model), intent(out) :: dist
    call set_common(dist,"negative-binomial",2,.true., &
      [character(len=8)::"size","mu"], &
      [tiny(1.0_dp),tiny(1.0_dp)],[huge(1.0_dp),huge(1.0_dp)])
    dist%logpdf => nbinom_logpdf
    dist%cdf => nbinom_cdf
    dist%quantile => nbinom_quantile
    dist%raw_moment => nbinom_moment
    dist%random_value => nbinom_random
  end subroutine make_negative_binomial

  pure function lower_ascii(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: i, code
    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) value(i:i)=achar(code+32)
    end do
  end function lower_ascii

  pure function normal_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, z
    if (size(par) /= 2 .or. par(2) <= 0.0_dp) then
      value = -huge(1.0_dp); return
    end if
    z = (x-par(1))/par(2)
    value = -0.5_dp*log(2.0_dp*pi_dp)-log(par(2))-0.5_dp*z*z
  end function normal_logpdf

  pure function normal_cdf_model(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (size(par) /= 2 .or. par(2) <= 0.0_dp) then
      value = 0.0_dp
    else
      value = normal_cdf((x-par(1))/par(2))
    end if
  end function normal_cdf_model

  function normal_quantile_model(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value = par(1)+par(2)*normal_quantile(prob)
  end function normal_quantile_model

  pure function normal_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    select case(order)
    case(0); value=1.0_dp
    case(1); value=par(1)
    case(2); value=par(1)**2+par(2)**2
    case(3); value=par(1)**3+3.0_dp*par(1)*par(2)**2
    case(4); value=par(1)**4+6.0_dp*par(1)**2*par(2)**2+3.0_dp*par(2)**4
    case default; value=huge(1.0_dp)
    end select
  end function normal_moment

  function normal_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value = par(1)+par(2)*random_normal()
  end function normal_random

  pure function lognormal_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, z
    if (x <= 0.0_dp .or. size(par)/=2 .or. par(2)<=0.0_dp) then
      value=-huge(1.0_dp); return
    end if
    z=(log(x)-par(1))/par(2)
    value=-log(x)-log(par(2))-0.5_dp*log(2.0_dp*pi_dp)-0.5_dp*z*z
  end function lognormal_logpdf

  pure function lognormal_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=0.0_dp) then
      value=0.0_dp
    else
      value=normal_cdf((log(x)-par(1))/par(2))
    end if
  end function lognormal_cdf

  function lognormal_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=exp(par(1)+par(2)*normal_quantile(prob))
  end function lognormal_quantile

  pure function lognormal_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=exp(real(order,dp)*par(1)+0.5_dp*real(order*order,dp)*par(2)**2)
  end function lognormal_moment

  function lognormal_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=exp(par(1)+par(2)*random_normal())
  end function lognormal_random

  pure function exponential_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<0.0_dp .or. par(1)<=0.0_dp) then
      value=-huge(1.0_dp)
    else
      value=log(par(1))-par(1)*x
    end if
  end function exponential_logpdf

  pure function exponential_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=0.0_dp) then
      value=0.0_dp
    else
      value=1.0_dp-exp(-par(1)*x)
    end if
  end function exponential_cdf

  pure function exponential_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    if (prob>=1.0_dp) then
      value=huge(1.0_dp)
    else
      value=-log(1.0_dp-prob)/par(1)
    end if
  end function exponential_quantile

  pure function exponential_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=gamma(real(order+1,dp))/par(1)**order
  end function exponential_moment

  function exponential_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, u
    call random_number(u)
    value=-log(max(1.0_dp-u,tiny(1.0_dp)))/par(1)
  end function exponential_random

  pure function gamma_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<0.0_dp .or. any(par<=0.0_dp)) then
      value=-huge(1.0_dp)
    else if (x<=tiny(1.0_dp) .and. par(1)<1.0_dp) then
      value=huge(1.0_dp)
    else
      value=par(1)*log(par(2))-log_gamma(par(1)) + &
        (par(1)-1.0_dp)*log(max(x,tiny(1.0_dp)))-par(2)*x
    end if
  end function gamma_logpdf

  pure function gamma_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=0.0_dp) then
      value=0.0_dp
    else
      value=regularized_gamma_p(par(1),par(2)*x)
    end if
  end function gamma_cdf

  function gamma_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=continuous_quantile(prob,par,gamma_cdf,0.0_dp,max(1.0_dp,par(1)/par(2)))
  end function gamma_quantile

  pure function gamma_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=exp(log_gamma(par(1)+real(order,dp))-log_gamma(par(1)))/par(2)**order
  end function gamma_moment

  function gamma_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=random_gamma(par(1),1.0_dp/par(2))
  end function gamma_random

  pure function weibull_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=0.0_dp .or. any(par<=0.0_dp)) then
      value=-huge(1.0_dp)
    else
      value=log(par(1))-log(par(2))+(par(1)-1.0_dp)*log(x/par(2))-(x/par(2))**par(1)
    end if
  end function weibull_logpdf

  pure function weibull_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=0.0_dp) then
      value=0.0_dp
    else
      value=1.0_dp-exp(-(x/par(2))**par(1))
    end if
  end function weibull_cdf

  pure function weibull_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=par(2)*(-log(1.0_dp-prob))**(1.0_dp/par(1))
  end function weibull_quantile

  pure function weibull_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=par(2)**order*gamma(1.0_dp+real(order,dp)/par(1))
  end function weibull_moment

  function weibull_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, u
    call random_number(u)
    value=weibull_quantile(u,par)
  end function weibull_random

  pure function uniform_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (par(2)<=par(1) .or. x<par(1) .or. x>par(2)) then
      value=-huge(1.0_dp)
    else
      value=-log(par(2)-par(1))
    end if
  end function uniform_logpdf

  pure function uniform_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=par(1)) then
      value=0.0_dp
    else if (x>=par(2)) then
      value=1.0_dp
    else
      value=(x-par(1))/(par(2)-par(1))
    end if
  end function uniform_cdf

  pure function uniform_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=par(1)+prob*(par(2)-par(1))
  end function uniform_quantile

  pure function uniform_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=(par(2)**(order+1)-par(1)**(order+1))/ &
      (real(order+1,dp)*(par(2)-par(1)))
  end function uniform_moment

  function uniform_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, u
    call random_number(u)
    value=uniform_quantile(u,par)
  end function uniform_random

  pure function logistic_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, z
    if (par(2)<=0.0_dp) then
      value=-huge(1.0_dp); return
    end if
    z=(x-par(1))/par(2)
    if (z>=0.0_dp) then
      value=-log(par(2))-z-2.0_dp*log(1.0_dp+exp(-z))
    else
      value=-log(par(2))+z-2.0_dp*log(1.0_dp+exp(z))
    end if
  end function logistic_logpdf

  pure function logistic_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, z
    z=(x-par(1))/par(2)
    if (z>=0.0_dp) then
      value=1.0_dp/(1.0_dp+exp(-z))
    else
      value=exp(z)/(1.0_dp+exp(z))
    end if
  end function logistic_cdf

  pure function logistic_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=par(1)+par(2)*log(prob/(1.0_dp-prob))
  end function logistic_quantile

  pure function logistic_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    select case(order)
    case(0); value=1.0_dp
    case(1); value=par(1)
    case(2); value=par(1)**2+pi_dp*pi_dp*par(2)**2/3.0_dp
    case default; value=huge(1.0_dp)
    end select
  end function logistic_moment

  function logistic_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, u
    call random_number(u)
    value=logistic_quantile(clamp_probability(u),par)
  end function logistic_random

  pure function cauchy_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, z
    if (par(2)<=0.0_dp) then
      value=-huge(1.0_dp); return
    end if
    z=(x-par(1))/par(2)
    value=-log(pi_dp*par(2))-log(1.0_dp+z*z)
  end function cauchy_logpdf

  pure function cauchy_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    value=0.5_dp+atan((x-par(1))/par(2))/pi_dp
  end function cauchy_cdf

  pure function cauchy_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=par(1)+par(2)*tan(pi_dp*(prob-0.5_dp))
  end function cauchy_quantile

  pure function cauchy_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    if (order==0) then
      value=1.0_dp + 0.0_dp*sum(par)
    else
      value=huge(1.0_dp)
    end if
  end function cauchy_moment

  function cauchy_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, u
    call random_number(u)
    value=cauchy_quantile(clamp_probability(u),par)
  end function cauchy_random

  pure function beta_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<=0.0_dp .or. x>=1.0_dp .or. any(par<=0.0_dp)) then
      value=-huge(1.0_dp)
    else
      value=(par(1)-1.0_dp)*log(x)+(par(2)-1.0_dp)*log(1.0_dp-x) + &
        log_gamma(par(1)+par(2))-log_gamma(par(1))-log_gamma(par(2))
    end if
  end function beta_logpdf

  pure function beta_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    value=regularized_beta(x,par(1),par(2))
  end function beta_cdf

  function beta_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    value=continuous_quantile(prob,par,beta_cdf,0.0_dp,1.0_dp)
  end function beta_quantile

  pure function beta_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=exp(log_gamma(par(1)+real(order,dp))-log_gamma(par(1)) + &
      log_gamma(par(1)+par(2))-log_gamma(par(1)+par(2)+real(order,dp)))
  end function beta_moment

  function beta_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, x, y
    x=random_gamma(par(1),1.0_dp)
    y=random_gamma(par(2),1.0_dp)
    value=x/(x+y)
  end function beta_random

  pure function poisson_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<0.0_dp .or. abs(x-nint(x))>1.0e-10_dp .or. par(1)<=0.0_dp) then
      value=-huge(1.0_dp)
    else
      value=x*log(par(1))-par(1)-log_gamma(x+1.0_dp)
    end if
  end function poisson_logpdf

  pure function poisson_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    integer :: k
    if (x<0.0_dp) then
      value=0.0_dp
    else
      k=floor(x)
      value=1.0_dp-regularized_gamma_p(real(k+1,dp),par(1))
      value=min(max(value,0.0_dp),1.0_dp)
    end if
  end function poisson_cdf

  function poisson_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    integer :: k
    k=0
    do while(poisson_cdf(real(k,dp),par)<prob .and. k<1000000)
      k=k+1
    end do
    value=real(k,dp)
  end function poisson_quantile

  pure function poisson_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value, l
    l=par(1)
    select case(order)
    case(0); value=1.0_dp
    case(1); value=l
    case(2); value=l*l+l
    case(3); value=l**3+3.0_dp*l*l+l
    case(4); value=l**4+6.0_dp*l**3+7.0_dp*l*l+l
    case default; value=huge(1.0_dp)
    end select
  end function poisson_moment

  function poisson_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value
    value=real(random_poisson(par(1)),dp)
  end function poisson_random

  pure function geometric_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<0.0_dp .or. abs(x-nint(x))>1.0e-10_dp .or. par(1)<=0.0_dp .or. par(1)>1.0_dp) then
      value=-huge(1.0_dp)
    else
      value=log(par(1))+x*log(1.0_dp-par(1))
    end if
  end function geometric_logpdf

  pure function geometric_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value
    if (x<0.0_dp) then
      value=0.0_dp
    else
      value=1.0_dp-(1.0_dp-par(1))**(floor(x)+1.0_dp)
    end if
  end function geometric_cdf

  pure function geometric_quantile(prob, par) result(value)
    real(dp), intent(in) :: prob, par(:)
    real(dp) :: value
    if (prob<=par(1)) then
      value=0.0_dp
    else
      value=ceiling(log(1.0_dp-prob)/log(1.0_dp-par(1))-1.0_dp)
    end if
  end function geometric_quantile

  pure function geometric_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value, p
    p=par(1)
    select case(order)
    case(0); value=1.0_dp
    case(1); value=(1.0_dp-p)/p
    case(2); value=(1.0_dp-p)*(2.0_dp-p)/p**2
    case default; value=huge(1.0_dp)
    end select
  end function geometric_moment

  function geometric_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, u
    call random_number(u)
    value=geometric_quantile(u,par)
  end function geometric_random

  pure function nbinom_logpdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, sizep, mu
    sizep=par(1); mu=par(2)
    if (x<0.0_dp .or. abs(x-nint(x))>1.0e-10_dp .or. sizep<=0.0_dp .or. mu<=0.0_dp) then
      value=-huge(1.0_dp)
    else
      value=log_gamma(sizep+x)-log_gamma(sizep)-log_gamma(x+1.0_dp) + &
        sizep*log(sizep/(sizep+mu))+x*log(mu/(sizep+mu))
    end if
  end function nbinom_logpdf

  pure function nbinom_cdf(x, par) result(value)
    real(dp), intent(in) :: x, par(:)
    real(dp) :: value, prob
    integer :: k
    if (x<0.0_dp) then
      value=0.0_dp
    else
      k=floor(x)
      prob=par(1)/(par(1)+par(2))
      value=regularized_beta(prob,par(1),real(k+1,dp))
    end if
  end function nbinom_cdf

  function nbinom_quantile(probability, par) result(value)
    real(dp), intent(in) :: probability, par(:)
    real(dp) :: value
    integer :: k
    k=0
    do while(nbinom_cdf(real(k,dp),par)<probability .and. k<1000000)
      k=k+1
    end do
    value=real(k,dp)
  end function nbinom_quantile

  pure function nbinom_moment(order, par) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: par(:)
    real(dp) :: value, mu, sizep
    sizep=par(1); mu=par(2)
    select case(order)
    case(0); value=1.0_dp
    case(1); value=mu
    case(2); value=mu*mu+mu+mu*mu/sizep
    case default; value=huge(1.0_dp)
    end select
  end function nbinom_moment

  function nbinom_random(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: value, lambda
    lambda=random_gamma(par(1),par(2)/par(1))
    value=real(random_poisson(lambda),dp)
  end function nbinom_random

  function continuous_quantile(prob, par, cdf_fun, lower, initial_upper) result(value)
    real(dp), intent(in) :: prob, par(:), lower, initial_upper
    interface
      pure function cdf_fun(x, par) result(cdf_value)
        import dp
        real(dp), intent(in) :: x, par(:)
        real(dp) :: cdf_value
      end function cdf_fun
    end interface
    real(dp) :: value, lo, hi, mid
    integer :: iter
    if (prob<=0.0_dp) then
      value=lower; return
    else if (prob>=1.0_dp) then
      value=huge(1.0_dp); return
    end if
    lo=lower; hi=max(initial_upper,lower+1.0_dp)
    do while(cdf_fun(hi,par)<prob .and. hi<huge(1.0_dp)/4.0_dp)
      hi=max(2.0_dp*hi,hi+1.0_dp)
    end do
    do iter=1,200
      mid=0.5_dp*(lo+hi)
      if (cdf_fun(mid,par)<prob) then
        lo=mid
      else
        hi=mid
      end if
      if (abs(hi-lo)<=1.0e-12_dp*max(1.0_dp,abs(mid))) exit
    end do
    value=0.5_dp*(lo+hi)
  end function continuous_quantile

end module fitdistrplus_distributions
