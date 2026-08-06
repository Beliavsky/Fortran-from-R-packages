! SPDX-License-Identifier: GPL-3.0-only
module mass_distribution
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rrcov_kinds, only : dp
  use rrcov_random, only : seed_random, random_normal
  use mass_types, only : density_fit_result, mass_success, mass_invalid_argument, mass_no_convergence
  use mass_math, only : pi_dp, normal_pdf, student_t_pdf, bfgs_minimize, numerical_hessian, &
    covariance_from_hessian, sample_variance, sample_sd, median_mass, type7_quantile, &
    digamma_mass, trigamma_mass
  implicit none
  private
  public :: fit_distribution, rnegbin, theta_ml, theta_mm, theta_md
  public :: negative_binomial_logpmf, gamma_shape_estimate
contains

  subroutine fit_distribution(x,distribution,result,start,maxit,tolerance)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in)::distribution
    type(density_fit_result),intent(out)::result
    real(dp),intent(in),optional::start(:),tolerance
    integer,intent(in),optional::maxit
    character(len=32)::dist
    real(dp),allocatable::par(:),hess(:,:),covt(:,:),jac(:,:),natural(:)
    real(dp)::fval,m,v,s,tol
    integer::status,it,mit,n,k
    dist=lower_string(trim(distribution));n=size(x);mit=400;if(present(maxit))mit=maxit
    tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
    if(n<1 .or. any(.not.ieee_is_finite(x)))then;result%status=mass_invalid_argument;return;end if
    m=sum(x)/real(n,dp);v=sum((x-m)**2)/real(n,dp);s=sqrt(max(v,tiny(1.0_dp)))
    select case(trim(dist))
    case("normal")
      allocate(result%estimates(2),result%covariance(2,2));result%estimates=[m,s]
      result%covariance=0.0_dp;result%covariance(1,1)=s*s/real(n,dp);result%covariance(2,2)=s*s/(2.0_dp*real(n,dp))
      result%log_likelihood=sum(-log(s)-0.5_dp*log(2.0_dp*pi_dp)-0.5_dp*((x-m)/s)**2)
      result%status=mass_success
    case("lognormal","log-normal")
      if(any(x<=0.0_dp))then;result%status=mass_invalid_argument;return;end if
      m=sum(log(x))/real(n,dp);s=sqrt(sum((log(x)-m)**2)/real(n,dp))
      allocate(result%estimates(2),result%covariance(2,2));result%estimates=[m,s];result%covariance=0.0_dp
      result%covariance(1,1)=s*s/real(n,dp);result%covariance(2,2)=s*s/(2.0_dp*real(n,dp))
      result%log_likelihood=sum(-log(x)-log(s)-0.5_dp*log(2.0_dp*pi_dp)-0.5_dp*((log(x)-m)/s)**2)
      result%status=mass_success
    case("poisson")
      if(any(x<0.0_dp))then;result%status=mass_invalid_argument;return;end if
      allocate(result%estimates(1),result%covariance(1,1));result%estimates=[m];result%covariance(1,1)=m/real(n,dp)
      result%log_likelihood=sum(x*log(max(m,tiny(1.0_dp)))-m-log_gamma(x+1.0_dp));result%status=mass_success
    case("exponential")
      if(any(x<0.0_dp) .or. m<=0.0_dp)then;result%status=mass_invalid_argument;return;end if
      allocate(result%estimates(1),result%covariance(1,1));result%estimates=[1.0_dp/m]
      result%covariance(1,1)=result%estimates(1)**2/real(n,dp)
      result%log_likelihood=real(n,dp)*log(result%estimates(1))-result%estimates(1)*sum(x);result%status=mass_success
    case("geometric")
      if(any(x<0.0_dp))then;result%status=mass_invalid_argument;return;end if
      allocate(result%estimates(1),result%covariance(1,1));result%estimates=[1.0_dp/(1.0_dp+m)]
      result%covariance(1,1)=result%estimates(1)**2*(1.0_dp-result%estimates(1))/real(n,dp)
      result%log_likelihood=sum(log(result%estimates(1))+x*log(1.0_dp-result%estimates(1)));result%status=mass_success
    case default
      call initial_parameters(dist,x,m,v,s,par,status)
      if(status/=mass_success)then;result%status=status;return;end if
      if(present(start))then
        if(size(start)/=size(par))then;result%status=mass_invalid_argument;return;end if
        call natural_to_internal(dist,start,par)
      end if
      call bfgs_minimize(objective,par,fval,status,it,mit,tol)
      call internal_to_natural(dist,par,natural,jac)
      call numerical_hessian(objective,par,hess);covt=covariance_from_hessian(hess,k)
      result%estimates=natural;result%covariance=matmul(jac,matmul(covt,transpose(jac)))
      result%log_likelihood=-fval;result%iterations=it;result%status=status
    end select
    result%distribution=dist;result%aic=-2.0_dp*result%log_likelihood+2.0_dp*real(size(result%estimates),dp)
  contains
    function objective(p) result(value)
      real(dp),intent(in)::p(:);real(dp)::value
      value=-distribution_log_likelihood(dist,x,p)
      if(.not.ieee_is_finite(value))value=huge(1.0_dp)/100.0_dp
    end function objective
  end subroutine fit_distribution

  subroutine initial_parameters(dist,x,m,v,s,par,status)
    character(len=*),intent(in)::dist
    real(dp),intent(in)::x(:),m,v,s
    real(dp),allocatable,intent(out)::par(:)
    integer,intent(out)::status
    real(dp)::shape,scale,q25,q75
    status=mass_success;q25=type7_quantile(x,0.25_dp);q75=type7_quantile(x,0.75_dp)
    select case(trim(dist))
    case("gamma")
      if(any(x<0.0_dp) .or. m<=0.0_dp .or. v<=0.0_dp)then;status=mass_invalid_argument;return;end if
      par=[log(m*m/v),log(m/v)]
    case("weibull")
      if(any(x<=0.0_dp))then;status=mass_invalid_argument;return;end if
      scale=max(sample_sd(log(x)),0.1_dp);shape=1.2_dp/scale;par=[log(shape),log(exp(sum(log(x))/real(size(x),dp)+0.572_dp/shape))]
    case("beta")
      if(any(x<=0.0_dp) .or. any(x>=1.0_dp) .or. v<=0.0_dp)then;status=mass_invalid_argument;return;end if
      shape=max(m*(1.0_dp-m)/v-1.0_dp,0.1_dp);par=[log(max(m*shape,0.1_dp)),log(max((1.0_dp-m)*shape,0.1_dp))]
    case("cauchy","logistic")
      par=[median_mass(x),log(max(0.5_dp*(q75-q25),0.1_dp*s))]
    case("negative binomial","negative-binomial")
      shape=merge(m*m/(v-m),100.0_dp,v>m);par=[log(max(shape,1.0e-3_dp)),log(max(m,1.0e-3_dp))]
    case("t","student-t")
      par=[median_mass(x),log(max(0.5_dp*(q75-q25),0.1_dp*s)),log(8.0_dp)]
    case("chi-squared","chisquared","chi-square")
      if(any(x<0.0_dp))then;status=mass_invalid_argument;return;end if
      par=[log(max(m,0.1_dp))]
    case("f")
      if(any(x<=0.0_dp))then;status=mass_invalid_argument;return;end if
      par=[log(5.0_dp),log(10.0_dp)]
    case default
      status=mass_invalid_argument
    end select
  end subroutine initial_parameters

  function distribution_log_likelihood(dist, x, p) result(value)
    character(len=*), intent(in) :: dist
    real(dp), intent(in) :: x(:), p(:)
    real(dp) :: value, a, b, loc, scale, nu, sizep, mu, z
    integer :: i

    value = 0.0_dp
    select case (trim(dist))
    case ("gamma")
      a = exp(p(1))
      b = exp(p(2))
      if (any(x < 0.0_dp)) then
        value = -huge(1.0_dp)
        return
      end if
      value = sum(a * log(b) - log_gamma(a) + &
        (a - 1.0_dp) * log(max(x, tiny(1.0_dp))) - b * x)
    case ("weibull")
      a = exp(p(1))
      b = exp(p(2))
      if (any(x <= 0.0_dp)) then
        value = -huge(1.0_dp)
        return
      end if
      value = sum(log(a) - log(b) + (a - 1.0_dp) * log(x / b) - &
        (x / b)**a)
    case ("beta")
      a = exp(p(1))
      b = exp(p(2))
      if (any(x <= 0.0_dp) .or. any(x >= 1.0_dp)) then
        value = -huge(1.0_dp)
        return
      end if
      value = sum((a - 1.0_dp) * log(x) + &
        (b - 1.0_dp) * log(1.0_dp - x)) + real(size(x), dp) * &
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b))
    case ("cauchy")
      loc = p(1)
      scale = exp(p(2))
      do i = 1, size(x)
        z = (x(i) - loc) / scale
        value = value - log(pi_dp * scale) - log(1.0_dp + z * z)
      end do
    case ("logistic")
      loc = p(1)
      scale = exp(p(2))
      do i = 1, size(x)
        z = (x(i) - loc) / scale
        if (z >= 0.0_dp) then
          value = value - log(scale) - z - &
            2.0_dp * log(1.0_dp + exp(-z))
        else
          value = value - log(scale) + z - &
            2.0_dp * log(1.0_dp + exp(z))
        end if
      end do
    case ("negative binomial", "negative-binomial")
      sizep = exp(p(1))
      mu = exp(p(2))
      do i = 1, size(x)
        value = value + negative_binomial_logpmf(x(i), mu, sizep)
      end do
    case ("t", "student-t")
      loc = p(1)
      scale = exp(p(2))
      nu = 2.0_dp + exp(p(3))
      do i = 1, size(x)
        z = (x(i) - loc) / scale
        value = value + log(max(student_t_pdf(z, nu), tiny(1.0_dp))) - &
          log(scale)
      end do
    case ("chi-squared", "chisquared", "chi-square")
      nu = exp(p(1))
      if (any(x < 0.0_dp)) then
        value = -huge(1.0_dp)
        return
      end if
      value = sum((0.5_dp * nu - 1.0_dp) * &
        log(max(x, tiny(1.0_dp))) - 0.5_dp * x - &
        0.5_dp * nu * log(2.0_dp) - log_gamma(0.5_dp * nu))
    case ("f")
      a = exp(p(1))
      b = exp(p(2))
      if (any(x <= 0.0_dp)) then
        value = -huge(1.0_dp)
        return
      end if
      value = sum(0.5_dp * a * log(a / b) + &
        (0.5_dp * a - 1.0_dp) * log(x) - &
        0.5_dp * (a + b) * log(1.0_dp + a * x / b) + &
        log_gamma(0.5_dp * (a + b)) - log_gamma(0.5_dp * a) - &
        log_gamma(0.5_dp * b))
    case default
      value = -huge(1.0_dp)
    end select
  end function distribution_log_likelihood

  subroutine internal_to_natural(dist,p,natural,jac)
    character(len=*),intent(in)::dist
    real(dp),intent(in)::p(:)
    real(dp),allocatable,intent(out)::natural(:),jac(:,:)
    integer::i,n
    n=size(p);allocate(natural(n),jac(n,n));jac=0.0_dp
    select case(trim(dist))
    case("cauchy","logistic")
      natural=[p(1),exp(p(2))];jac(1,1)=1.0_dp;jac(2,2)=natural(2)
    case("t","student-t")
      natural=[p(1),exp(p(2)),2.0_dp+exp(p(3))];jac(1,1)=1.0_dp;jac(2,2)=natural(2);jac(3,3)=natural(3)-2.0_dp
    case default
      natural=exp(p);do i=1,n;jac(i,i)=natural(i);end do
    end select
  end subroutine internal_to_natural

  subroutine natural_to_internal(dist,natural,p)
    character(len=*),intent(in)::dist
    real(dp),intent(in)::natural(:)
    real(dp),intent(out)::p(:)
    select case(trim(dist))
    case("cauchy","logistic")
      p=[natural(1),log(max(natural(2),tiny(1.0_dp)))]
    case("t","student-t")
      p=[natural(1),log(max(natural(2),tiny(1.0_dp))),log(max(natural(3)-2.0_dp,tiny(1.0_dp)))]
    case default
      p=log(max(natural,tiny(1.0_dp)))
    end select
  end subroutine natural_to_internal

  pure function negative_binomial_logpmf(y,mu,theta) result(value)
    real(dp),intent(in)::y,mu,theta
    real(dp)::value
    if(y<0.0_dp .or. mu<=0.0_dp .or. theta<=0.0_dp)then;value=-huge(1.0_dp);return;end if
    value=log_gamma(theta+y)-log_gamma(theta)-log_gamma(y+1.0_dp)+theta*log(theta/(theta+mu))+y*log(mu/(theta+mu))
  end function negative_binomial_logpmf

  subroutine rnegbin(n,mu,theta,values,seed,status)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),theta
    integer,allocatable,intent(out)::values(:)
    integer,intent(in),optional::seed
    integer,intent(out),optional::status
    real(dp)::m,lambda
    integer::i
    if(n<1 .or. size(mu)<1 .or. theta<=0.0_dp .or. any(mu<0.0_dp))then
      allocate(values(0));if(present(status))status=mass_invalid_argument;return
    end if
    call seed_random(seed);allocate(values(n))
    do i=1,n
      m=mu(1+mod(i-1,size(mu)));lambda=random_gamma(theta,m/theta);values(i)=random_poisson(lambda)
    end do
    if(present(status))status=mass_success
  end subroutine rnegbin

  recursive function random_gamma(shape, scale) result(value)
    real(dp),intent(in)::shape,scale
    real(dp)::value,d,c,x,v,u
    if(shape<=0.0_dp .or. scale<0.0_dp)then;value=0.0_dp;return;end if
    if(shape<1.0_dp)then
      call random_number(u);value=random_gamma(shape+1.0_dp,scale)*u**(1.0_dp/shape);return
    end if
    d=shape-1.0_dp/3.0_dp;c=1.0_dp/sqrt(9.0_dp*d)
    do
      x=random_normal();v=(1.0_dp+c*x)**3;if(v<=0.0_dp)cycle
      call random_number(u)
      if(u<1.0_dp-0.0331_dp*x**4 .or. log(u)<0.5_dp*x*x+d*(1.0_dp-v+log(v)))exit
    end do
    value=scale*d*v
  end function random_gamma

  function random_poisson(lambda) result(value)
    real(dp),intent(in)::lambda
    integer::value,k
    real(dp)::l,p,u,z
    if(lambda<=0.0_dp)then;value=0;return;end if
    if(lambda<30.0_dp)then
      l=exp(-lambda);k=0;p=1.0_dp
      do;k=k+1;call random_number(u);p=p*u;if(p<=l)exit;end do
      value=k-1
    else
      do
        z=lambda+sqrt(lambda)*random_normal();value=nint(z)
        if(value>=0)exit
      end do
    end if
  end function random_poisson

  subroutine theta_ml(y,mu,theta,se,status,weights,limit,tolerance)
    real(dp),intent(in)::y(:),mu(:)
    real(dp),intent(out)::theta,se
    integer,intent(out)::status
    real(dp),intent(in),optional::weights(:),tolerance
    integer,intent(in),optional::limit
    real(dp),allocatable::w(:)
    real(dp)::eps,del,score,info,nw
    integer::it,mit
    if (size(y) /= size(mu) .or. any(mu <= 0.0_dp) .or. &
        any(y < 0.0_dp)) then
      theta = 0.0_dp
      se = 0.0_dp
      status = mass_invalid_argument
      return
    end if
    allocate(w(size(y)));w=1.0_dp;if(present(weights))w=weights
    mit=10;if(present(limit))mit=limit;eps=epsilon(1.0_dp)**0.25_dp;if(present(tolerance))eps=tolerance;nw=sum(w)
    theta=nw/max(sum(w*(y/mu-1.0_dp)**2),tiny(1.0_dp));info=1.0_dp
    do it=1,mit
      theta=abs(theta)
      score=sum(w*(digamma_mass(theta+y)-digamma_mass(theta)+log(theta)+1.0_dp-log(theta+mu)-(y+theta)/(mu+theta)))
      info=sum(w*(-trigamma_mass(theta+y)+trigamma_mass(theta)-1.0_dp/theta+2.0_dp/(mu+theta)-(y+theta)/(mu+theta)**2))
      del=score/max(info,tiny(1.0_dp));theta=theta+del;if(abs(del)<=eps)exit
    end do
    theta=max(theta,0.0_dp);se=sqrt(1.0_dp/max(info,tiny(1.0_dp)));status=merge(mass_success,mass_no_convergence,it<=mit)
  end subroutine theta_ml

  subroutine theta_mm(y,mu,dfr,theta,status,weights,limit,tolerance)
    real(dp),intent(in)::y(:),mu(:),dfr
    real(dp),intent(out)::theta
    integer,intent(out)::status
    real(dp),intent(in),optional::weights(:),tolerance
    integer,intent(in),optional::limit
    real(dp),allocatable::w(:)
    real(dp)::eps,del,nw
    integer::it,mit
    if(size(y)/=size(mu).or.any(mu<=0.0_dp))then;theta=0.0_dp;status=mass_invalid_argument;return;end if
    allocate(w(size(y)));w=1.0_dp;if(present(weights))w=weights
    mit=10;if(present(limit))mit=limit;eps=epsilon(1.0_dp)**0.25_dp;if(present(tolerance))eps=tolerance;nw=sum(w)
    theta=nw/max(sum(w*(y/mu-1.0_dp)**2),tiny(1.0_dp))
    do it=1,mit
      theta=abs(theta);del=(sum(w*(y-mu)**2/(mu+mu*mu/theta))-dfr)/max(sum(w*(y-mu)**2/(mu+theta)**2),tiny(1.0_dp))
      theta=theta-del;if(abs(del)<=eps)exit
    end do
    theta=max(theta,0.0_dp);status=merge(mass_success,mass_no_convergence,it<=mit)
  end subroutine theta_mm

  subroutine theta_md(y,mu,dfr,theta,status,weights,limit,tolerance)
    real(dp),intent(in)::y(:),mu(:),dfr
    real(dp),intent(out)::theta
    integer,intent(out)::status
    real(dp),intent(in),optional::weights(:),tolerance
    integer,intent(in),optional::limit
    real(dp),allocatable::w(:),tmp(:)
    real(dp)::eps,del,nw,a,top,bot
    integer::it,mit
    if(size(y)/=size(mu).or.any(mu<=0.0_dp))then;theta=0.0_dp;status=mass_invalid_argument;return;end if
    allocate(w(size(y)));w=1.0_dp;if(present(weights))w=weights
    mit=20;if(present(limit))mit=limit;eps=epsilon(1.0_dp)**0.25_dp;if(present(tolerance))eps=tolerance;nw=sum(w)
    theta=nw/max(sum(w*(y/mu-1.0_dp)**2),tiny(1.0_dp));a=2.0_dp*sum(w*y*log(max(1.0_dp,y)/mu))-dfr
    do it=1,mit
      theta=abs(theta);tmp=log((y+theta)/(mu+theta));top=a-2.0_dp*sum(w*(y+theta)*tmp)
      bot=2.0_dp*sum(w*((y-mu)/(mu+theta)-tmp));del=top/max(abs(bot),tiny(1.0_dp))*sign(1.0_dp,bot)
      theta=theta-del;if(abs(del)<=eps)exit
    end do
    theta=max(theta,0.0_dp);status=merge(mass_success,mass_no_convergence,it<=mit)
  end subroutine theta_md

  subroutine gamma_shape_estimate(y,mu,alpha,se,status,weights,deviance,df_residual,limit,tolerance)
    real(dp),intent(in)::y(:),mu(:)
    real(dp),intent(out)::alpha,se
    integer,intent(out)::status
    real(dp),intent(in),optional::weights(:),deviance,df_residual,tolerance
    integer,intent(in),optional::limit
    real(dp),allocatable::w(:),fixed(:)
    real(dp)::dbar,eps,step,score,info
    integer::it,mit
    if(size(y)/=size(mu).or.any(mu<=0.0_dp).or.any(y<0.0_dp))then;alpha=0.0_dp;se=0.0_dp;status=mass_invalid_argument;return;end if
    allocate(w(size(y)));w=1.0_dp;if(present(weights))w=weights
    if(present(deviance).and.present(df_residual))then;dbar=deviance/df_residual
    else;dbar=sum(w*(y-mu)**2/max(mu*mu,tiny(1.0_dp)))/max(sum(w)-1.0_dp,1.0_dp);end if
    alpha = (6.0_dp + 2.0_dp * dbar) / &
      max(dbar * (6.0_dp + dbar), tiny(1.0_dp))
    fixed = -y / mu - log(mu) + log(max(w, tiny(1.0_dp))) + 1.0_dp + &
      log(y + merge(1.0_dp, 0.0_dp, y <= tiny(1.0_dp)))
    mit=10;if(present(limit))mit=limit;eps=epsilon(1.0_dp)**0.25_dp;if(present(tolerance))eps=tolerance
    do it=1,mit
      score=sum(w*(fixed+log(alpha)-digamma_mass(w*alpha)));info=sum(w*(w*trigamma_mass(w*alpha)-1.0_dp/alpha))
      step=score/max(info,tiny(1.0_dp));alpha=alpha+step;if(abs(step)<=eps)exit
    end do
    se=sqrt(1.0_dp/max(info,tiny(1.0_dp)));status=merge(mass_success,mass_no_convergence,it<=mit)
  end subroutine gamma_shape_estimate

  pure function lower_string(s) result(out)
    character(len=*),intent(in)::s
    character(len=len(s))::out
    integer::i,c
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        out(i:i) = achar(c + 32)
      else
        out(i:i) = s(i:i)
      end if
    end do
  end function lower_string

end module mass_distribution
