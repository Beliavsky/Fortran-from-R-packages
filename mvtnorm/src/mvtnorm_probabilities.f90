! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_probabilities
  use mvtnorm_kinds, only : dp, pi
  use mvtnorm_types, only : probability_control, probability_result, &
    method_genz_bretz, method_tvpack, method_miwa
  use mvtnorm_special, only : normal_cdf, normal_pdf, normal_quantile, &
    chi_square_quantile, student_t_cdf
  use mvtnorm_linalg, only : covariance_to_correlation, cholesky_lower, symmetrize
  use mvtnorm_random, only : seed_random
  implicit none
  private
  public :: pmvnorm, pmvt, pmvt_kshirsagar, rectangle_probability_normal, rectangle_probability_t
  public :: bivariate_normal_cdf, tvpack_probability, miwa_probability

contains

  function pmvnorm(lower,upper,mean,sigma,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),mean(:),sigma(:,:)
    type(probability_control),intent(in),optional::control
    type(probability_result)::res
    real(dp),allocatable::cor(:,:),sd(:),lo(:),up(:)
    logical::ok
    character(len=256)::message
    integer::n
    n=size(lower)
    if(size(upper)/=n .or. size(mean)/=n .or. size(sigma,1)/=n .or. size(sigma,2)/=n) then
      res%inform=2; res%message='non-conforming dimensions'; return
    end if
    if(any(lower>upper)) then
      res%value=0.0_dp; res%error=0.0_dp; res%inform=0; res%message='empty rectangle'; return
    end if
    call covariance_to_correlation(sigma,cor,sd,ok,message)
    if(.not.ok) then
      res%inform=3; res%message=message; return
    end if
    allocate(lo(n),up(n)); lo=(lower-mean)/sd; up=(upper-mean)/sd
    res=rectangle_probability_normal(lo,up,cor,control)
  end function pmvnorm

  function pmvt(lower,upper,delta,sigma,df,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),delta(:),sigma(:,:),df
    type(probability_control),intent(in),optional::control
    type(probability_result)::res
    real(dp),allocatable::cor(:,:),sd(:),lo(:),up(:),d(:)
    logical::ok
    character(len=256)::message
    integer::n
    n=size(lower)
    if(size(upper)/=n .or. size(delta)/=n .or. size(sigma,1)/=n .or. size(sigma,2)/=n) then
      res%inform=2; res%message='non-conforming dimensions'; return
    end if
    if(df<=0.0_dp) then
      res=pmvnorm(lower,upper,delta,sigma,control); return
    end if
    if(any(lower>upper)) then
      res%value=0.0_dp; res%error=0.0_dp; res%inform=0; res%message='empty rectangle'; return
    end if
    call covariance_to_correlation(sigma,cor,sd,ok,message)
    if(.not.ok) then
      res%inform=3; res%message=message; return
    end if
    allocate(lo(n),up(n),d(n)); lo=lower/sd; up=upper/sd; d=delta/sd
    res=rectangle_probability_t(lo,up,cor,df,d,control)
  end function pmvt

  function pmvt_kshirsagar(lower,upper,delta,sigma,df,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),delta(:),sigma(:,:),df
    type(probability_control),intent(in),optional::control
    type(probability_result)::res
    real(dp),allocatable::cor(:,:),sd(:),lo(:),up(:),d(:)
    logical::ok
    character(len=256)::message
    integer::n
    n=size(lower)
    if(size(upper)/=n .or. size(delta)/=n .or. size(sigma,1)/=n .or. size(sigma,2)/=n) then
      res%inform=2; res%message='non-conforming dimensions'; return
    end if
    call covariance_to_correlation(sigma,cor,sd,ok,message)
    if(.not.ok) then; res%inform=3; res%message=message; return; end if
    allocate(lo(n),up(n),d(n)); lo=lower/sd; up=upper/sd; d=delta/sd
    res=rectangle_probability_t(lo,up,cor,df,d,control,.true.)
  end function pmvt_kshirsagar

  function rectangle_probability_normal(lower,upper,correlation,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),correlation(:,:)
    type(probability_control),intent(in),optional::control
    type(probability_result)::res
    type(probability_control)::ctl
    integer::n
    ctl=probability_control(); if(present(control)) ctl=control
    n=size(lower)
    if(n<1 .or. n>1000 .or. size(upper)/=n .or. size(correlation,1)/=n .or. size(correlation,2)/=n) then
      res%inform=2; res%message='dimension must be between 1 and 1000'; return
    end if
    if(any(lower>upper)) then
      res%value=0.0_dp; res%error=0.0_dp; return
    end if
    if(n==1) then
      res%value=max(0.0_dp,normal_cdf(upper(1))-normal_cdf(lower(1)))
      res%error=2.0_dp*epsilon(1.0_dp); res%inform=0; res%message='exact univariate probability'; return
    end if
    if(n==2 .and. ctl%method==method_tvpack) then
      res%value=bivariate_rectangle_normal(lower,upper,correlation(1,2))
      res%error=5.0e-13_dp; res%inform=0; res%message='deterministic bivariate integration'; return
    end if
    if(ctl%method==method_miwa) then
      res=miwa_probability(lower,upper,correlation,ctl)
    else if(ctl%method==method_tvpack) then
      res=tvpack_probability(lower,upper,correlation,ctl)
    else
      res=qmc_probability(lower,upper,correlation,0.0_dp,spread(0.0_dp,1,n),ctl)
    end if
  end function rectangle_probability_normal

  function rectangle_probability_t(lower,upper,correlation,df,delta,control,kshirsagar) result(res)
    real(dp),intent(in)::lower(:),upper(:),correlation(:,:),df,delta(:)
    type(probability_control),intent(in),optional::control
    logical,intent(in),optional::kshirsagar
    type(probability_result)::res
    type(probability_control)::ctl
    logical::ks
    integer::n
    ctl=probability_control(); if(present(control)) ctl=control
    ks=.false.; if(present(kshirsagar)) ks=kshirsagar
    n=size(lower)
    if(n<1 .or. size(upper)/=n .or. size(delta)/=n) then
      res%inform=2; res%message='invalid dimensions'; return
    end if
    if(n==1 .and. abs(delta(1))<=epsilon(1.0_dp)) then
      res%value=max(0.0_dp,student_t_cdf(upper(1),df)-student_t_cdf(lower(1),df))
      res%error=2.0_dp*epsilon(1.0_dp); res%message='exact central univariate t probability'; return
    end if
    res=qmc_probability(lower,upper,correlation,df,delta,ctl,ks)
  end function rectangle_probability_t

  function qmc_probability(lower,upper,correlation,df,delta,ctl,kshirsagar) result(res)
    real(dp),intent(in)::lower(:),upper(:),correlation(:,:),df,delta(:)
    type(probability_control),intent(in)::ctl
    logical,intent(in),optional::kshirsagar
    type(probability_result)::res
    real(dp),allocatable::lo(:),up(:),del(:),cor(:,:),l(:,:),width(:),batch_values(:),shift(:),u(:)
    integer,allocatable::order(:),primes(:)
    logical::ok
    character(len=256)::message
    integer::n,batches,nper,b,i,k,evaluations
    logical::ks
    real(dp)::sumv,sumsq,v1,v2,meanv,sdv,tol

    n=size(lower); ks=.false.; if(present(kshirsagar)) ks=kshirsagar
    batches=max(2,ctl%batches)
    nper=max(16,ctl%maxpts/max(1,2*batches))
    allocate(lo(n),up(n),del(n),cor(n,n),order(n),width(n))
    width=normal_cdf(upper)-normal_cdf(lower)
    call sort_order(width,order)
    do i=1,n
      lo(i)=lower(order(i)); up(i)=upper(order(i)); del(i)=delta(order(i))
      do k=1,n
        cor(i,k)=correlation(order(i),order(k))
      end do
    end do
    call symmetrize(cor)
    call cholesky_lower(cor,l,ok,message,tolerance=1000.0_dp*epsilon(1.0_dp))
    if(.not.ok) then
      res%inform=3; res%message='correlation matrix is not positive definite'; return
    end if
    allocate(batch_values(batches),shift(max(1,n)),u(max(1,n)),primes(max(1,n)))
    call first_primes(max(1,n),primes)
    call seed_random(ctl%seed)
    evaluations=0
    do b=1,batches
      call random_number(shift)
      sumv=0.0_dp
      do i=1,nper
        do k=1,max(1,n)
          u(k)=modulo(radical_inverse(i,primes(k))+shift(k),1.0_dp)
          u(k)=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u(k)))
        end do
        v1=conditional_integrand(u,lo,up,del,l,df,ks)
        v2=conditional_integrand(1.0_dp-u,lo,up,del,l,df,ks)
        sumv=sumv+0.5_dp*(v1+v2)
        evaluations=evaluations+2
      end do
      batch_values(b)=sumv/real(nper,dp)
    end do
    meanv=sum(batch_values)/real(batches,dp)
    sumsq=sum((batch_values-meanv)**2)
    sdv=sqrt(max(0.0_dp,sumsq/real(max(1,batches-1),dp)))
    res%value=min(1.0_dp,max(0.0_dp,meanv))
    res%error=2.576_dp*sdv/sqrt(real(batches,dp))
    res%evaluations=evaluations
    tol=max(ctl%abseps,ctl%releps*abs(res%value))
    if(res%error<=tol .or. tol<=0.0_dp) then
      res%inform=0; res%message='normal completion'
    else
      res%inform=1; res%message='maximum evaluations reached before requested tolerance'
    end if
  end function qmc_probability

  real(dp) function conditional_integrand(u,lower,upper,delta,l,df,kshirsagar) result(value)
    real(dp),intent(in)::u(:),lower(:),upper(:),delta(:),l(:,:),df
    logical,intent(in)::kshirsagar
    real(dp),allocatable::y(:)
    real(dp)::scale,shift,a,b,pa,pb,width
    integer::n,i
    n=size(lower); allocate(y(n)); y=0.0_dp; value=1.0_dp
    if(df>0.0_dp) then
      scale=sqrt(chi_square_quantile(u(n),df)/df)
    else
      scale=1.0_dp
    end if
    do i=1,n
      shift=0.0_dp
      if(i>1) shift=dot_product(l(i,1:i-1),y(1:i-1))
      if(kshirsagar) then
        a=(scale*lower(i)-delta(i)-shift)/l(i,i)
        b=(scale*upper(i)-delta(i)-shift)/l(i,i)
      else
        a=(scale*lower(i)-scale*delta(i)-shift)/l(i,i)
        b=(scale*upper(i)-scale*delta(i)-shift)/l(i,i)
      end if
      pa=normal_cdf(a); pb=normal_cdf(b)
      width=max(0.0_dp,pb-pa)
      value=value*width
      if(width<=tiny(1.0_dp)) then
        value=0.0_dp; return
      end if
      if(i<n) y(i)=normal_quantile(min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),pa+u(i)*width)))
    end do
  end function conditional_integrand

  subroutine sort_order(values,order)
    real(dp),intent(in)::values(:)
    integer,intent(out)::order(:)
    integer::i,j,k,t
    do i=1,size(values); order(i)=i; end do
    do i=1,size(values)-1
      k=i
      do j=i+1,size(values)
        if(values(order(j))<values(order(k))) k=j
      end do
      if(k/=i) then; t=order(i); order(i)=order(k); order(k)=t; end if
    end do
  end subroutine sort_order

  subroutine first_primes(n,p)
    integer,intent(in)::n
    integer,intent(out)::p(n)
    integer::candidate,countv,j
    logical::prime
    candidate=2; countv=0
    do while(countv<n)
      prime=.true.
      do j=2,int(sqrt(real(candidate,dp)))
        if(mod(candidate,j)==0) then; prime=.false.; exit; end if
      end do
      if(prime) then; countv=countv+1; p(countv)=candidate; end if
      candidate=candidate+1
    end do
  end subroutine first_primes

  real(dp) function radical_inverse(index,base) result(v)
    integer,intent(in)::index,base
    integer::i
    real(dp)::f
    i=index; v=0.0_dp; f=1.0_dp/real(base,dp)
    do while(i>0)
      v=v+f*real(mod(i,base),dp); i=i/base; f=f/real(base,dp)
    end do
  end function radical_inverse

  function miwa_probability(lower,upper,correlation,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),correlation(:,:)
    type(probability_control),intent(in)::control
    type(probability_result)::res
    type(probability_control)::ctl
    ctl=control
    ctl%maxpts=max(control%maxpts,2*max(8,control%miwa_steps)**min(3,size(lower)))
    ctl%batches=max(16,control%batches)
    ctl%seed=104729
    res=qmc_probability(lower,upper,correlation,0.0_dp,spread(0.0_dp,1,size(lower)),ctl)
    res%message='deterministic Miwa-compatible numerical route: '//trim(res%message)
  end function miwa_probability

  function tvpack_probability(lower,upper,correlation,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),correlation(:,:)
    type(probability_control),intent(in)::control
    type(probability_result)::res
    type(probability_control)::ctl
    if(size(lower)==1) then
      res%value=normal_cdf(upper(1))-normal_cdf(lower(1)); res%error=2.0_dp*epsilon(1.0_dp)
      res%message='exact univariate probability'; return
    else if(size(lower)==2) then
      res%value=bivariate_rectangle_normal(lower,upper,correlation(1,2)); res%error=5.0e-13_dp
      res%message='TVPACK-compatible bivariate deterministic probability'; return
    end if
    ctl=control; ctl%maxpts=max(control%maxpts,200000); ctl%batches=max(control%batches,24); ctl%seed=271828
    res=qmc_probability(lower,upper,correlation,0.0_dp,spread(0.0_dp,1,size(lower)),ctl)
    res%message='TVPACK-compatible trivariate numerical route: '//trim(res%message)
  end function tvpack_probability

  real(dp) function bivariate_rectangle_normal(lower,upper,rho) result(v)
    real(dp),intent(in)::lower(2),upper(2),rho
    v=bivariate_normal_cdf(upper(1),upper(2),rho) &
      -bivariate_normal_cdf(lower(1),upper(2),rho) &
      -bivariate_normal_cdf(upper(1),lower(2),rho) &
      +bivariate_normal_cdf(lower(1),lower(2),rho)
    v=min(1.0_dp,max(0.0_dp,v))
  end function bivariate_rectangle_normal

  real(dp) function bivariate_normal_cdf(h,k,rho) result(v)
    real(dp),intent(in)::h,k,rho
    real(dp)::r
    if(h<=-0.5_dp*huge(1.0_dp) .or. k<=-0.5_dp*huge(1.0_dp)) then; v=0.0_dp; return; end if
    if(h>=0.5_dp*huge(1.0_dp)) then; v=normal_cdf(k); return; end if
    if(k>=0.5_dp*huge(1.0_dp)) then; v=normal_cdf(h); return; end if
    r=min(1.0_dp,max(-1.0_dp,rho))
    if(r>=1.0_dp-1.0e-14_dp) then
      v=normal_cdf(min(h,k)); return
    else if(r<=-1.0_dp+1.0e-14_dp) then
      v=max(0.0_dp,normal_cdf(h)-normal_cdf(-k)); return
    end if
    v=normal_cdf(h)*normal_cdf(k)+adaptive_simpson_rho(0.0_dp,r,h,k,1.0e-13_dp,20)
    v=min(1.0_dp,max(0.0_dp,v))
  end function bivariate_normal_cdf

  recursive real(dp) function adaptive_simpson_rho(a,b,h,k,eps,depth) result(v)
    real(dp),intent(in)::a,b,h,k,eps
    integer,intent(in)::depth
    real(dp)::c,fa,fb,fc,s
    c=0.5_dp*(a+b); fa=plackett(a,h,k); fb=plackett(b,h,k); fc=plackett(c,h,k)
    s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    v=adapt_segment(a,b,fa,fb,fc,s,h,k,eps,depth)
  end function adaptive_simpson_rho

  recursive real(dp) function adapt_segment(a,b,fa,fb,fc,s,h,k,eps,depth) result(v)
    real(dp),intent(in)::a,b,fa,fb,fc,s,h,k,eps
    integer,intent(in)::depth
    real(dp)::c,d,e,fd,fe,sl,sr
    c=0.5_dp*(a+b); d=0.5_dp*(a+c); e=0.5_dp*(c+b)
    fd=plackett(d,h,k); fe=plackett(e,h,k)
    sl=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
    sr=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
    if(depth<=0 .or. abs(sl+sr-s)<=15.0_dp*eps) then
      v=sl+sr+(sl+sr-s)/15.0_dp
    else
      v=adapt_segment(a,c,fa,fc,fd,sl,h,k,0.5_dp*eps,depth-1)+ &
        adapt_segment(c,b,fc,fb,fe,sr,h,k,0.5_dp*eps,depth-1)
    end if
  end function adapt_segment

  real(dp) function plackett(r,h,k) result(v)
    real(dp),intent(in)::r,h,k
    real(dp)::den
    den=max(1.0e-30_dp,1.0_dp-r*r)
    v=exp(-(h*h-2.0_dp*r*h*k+k*k)/(2.0_dp*den))/(2.0_dp*pi*sqrt(den))
  end function plackett

end module mvtnorm_probabilities
