! SPDX-License-Identifier: GPL-3.0-only
module anmc_math
  use anmc_kinds, only : dp, i8
  implicit none
  private

  type, public :: probability_control
    integer :: maxpts = 25000
    integer :: batches = 12
    real(dp) :: abseps = 1.0e-3_dp
    real(dp) :: releps = 0.0_dp
    integer :: seed = 12345
  end type probability_control

  type, public :: probability_result
    real(dp) :: value = 0.0_dp
    real(dp) :: error = 1.0_dp
    integer :: inform = 0
    integer :: evaluations = 0
    character(len=256) :: message = ''
  end type probability_result

  public :: genz_bretz, normal_cdf, normal_quantile, random_normals
  public :: cholesky_lower, inverse_spd, pmvnorm

contains

  function genz_bretz(maxpts,abseps,releps,seed,batches) result(control)
    integer,intent(in),optional::maxpts,seed,batches
    real(dp),intent(in),optional::abseps,releps
    type(probability_control)::control
    if(present(maxpts)) control%maxpts=maxpts
    if(present(abseps)) control%abseps=abseps
    if(present(releps)) control%releps=releps
    if(present(seed)) control%seed=seed
    if(present(batches)) control%batches=batches
  end function genz_bretz

  pure real(dp) function normal_cdf(x) result(p)
    real(dp),intent(in)::x
    p=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    ! Peter J. Acklam's rational approximation, expressed independently.
    real(dp),intent(in)::p
    real(dp),parameter::a1=-3.969683028665376e+01_dp,a2=2.209460984245205e+02_dp
    real(dp),parameter::a3=-2.759285104469687e+02_dp,a4=1.383577518672690e+02_dp
    real(dp),parameter::a5=-3.066479806614716e+01_dp,a6=2.506628277459239e+00_dp
    real(dp),parameter::b1=-5.447609879822406e+01_dp,b2=1.615858368580409e+02_dp
    real(dp),parameter::b3=-1.556989798598866e+02_dp,b4=6.680131188771972e+01_dp
    real(dp),parameter::b5=-1.328068155288572e+01_dp
    real(dp),parameter::c1=-7.784894002430293e-03_dp,c2=-3.223964580411365e-01_dp
    real(dp),parameter::c3=-2.400758277161838e+00_dp,c4=-2.549732539343734e+00_dp
    real(dp),parameter::c5=4.374664141464968e+00_dp,c6=2.938163982698783e+00_dp
    real(dp),parameter::d1=7.784695709041462e-03_dp,d2=3.224671290700398e-01_dp
    real(dp),parameter::d3=2.445134137142996e+00_dp,d4=3.754408661907416e+00_dp
    real(dp),parameter::plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp)::q,r
    if(p<=0.0_dp) then
      x=-huge(1.0_dp)
    else if(p>=1.0_dp) then
      x=huge(1.0_dp)
    else if(p<plow) then
      q=sqrt(-2.0_dp*log(p))
      x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if(p<=phigh) then
      q=p-0.5_dp; r=q*q
      x=((((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q)/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q=sqrt(-2.0_dp*log(1.0_dp-p))
      x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function normal_quantile

  subroutine random_normals(z)
    real(dp),intent(out)::z(:)
    real(dp)::u1,u2,r,pi
    integer::i
    pi=acos(-1.0_dp); i=1
    do while(i<=size(z))
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp)); r=sqrt(-2.0_dp*log(u1))
      z(i)=r*cos(2.0_dp*pi*u2)
      if(i+1<=size(z)) z(i+1)=r*sin(2.0_dp*pi*u2)
      i=i+2
    end do
  end subroutine random_normals

  subroutine cholesky_lower(a,l,ok,message,tolerance)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::l(:,:)
    logical,intent(out)::ok
    character(len=*),intent(out),optional::message
    real(dp),intent(in),optional::tolerance
    real(dp)::s,tol
    integer::n,i,j,k
    n=size(a,1); allocate(l(n,n)); l=0.0_dp; ok=.false.
    if(size(a,2)/=n) then
      if(present(message)) message='matrix must be square'; return
    end if
    tol=max(100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))),0.0_dp)
    if(present(tolerance)) tol=max(tolerance,0.0_dp)
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1; s=s-l(i,k)*l(j,k); end do
        if(i==j) then
          if(s<=tol) then
            if(present(message)) message='matrix is not positive definite'; return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
    ok=.true.; if(present(message)) message=''
  end subroutine cholesky_lower

  subroutine inverse_spd(a,ainv,ok,message)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    logical,intent(out)::ok
    character(len=*),intent(out),optional::message
    real(dp),allocatable::l(:,:),y(:),x(:)
    character(len=256)::msg
    integer::n,i,j,k
    call cholesky_lower(a,l,ok,msg)
    if(.not.ok) then
      allocate(ainv(size(a,1),size(a,2))); ainv=0.0_dp
      if(present(message)) message=msg; return
    end if
    n=size(a,1); allocate(ainv(n,n),y(n),x(n)); ainv=0.0_dp
    do j=1,n
      y=0.0_dp
      do i=1,n
        y(i)=merge(1.0_dp,0.0_dp,i==j)
        do k=1,i-1; y(i)=y(i)-l(i,k)*y(k); end do
        y(i)=y(i)/l(i,i)
      end do
      x=0.0_dp
      do i=n,1,-1
        x(i)=y(i)
        do k=i+1,n; x(i)=x(i)-l(k,i)*x(k); end do
        x(i)=x(i)/l(i,i)
      end do
      ainv(:,j)=x
    end do
    ainv=0.5_dp*(ainv+transpose(ainv)); ok=.true.
    if(present(message)) message=''
  end subroutine inverse_spd

  function pmvnorm(lower,upper,mean,sigma,control) result(res)
    real(dp),intent(in)::lower(:),upper(:),mean(:),sigma(:,:)
    type(probability_control),intent(in),optional::control
    type(probability_result)::res
    type(probability_control)::ctl
    real(dp),allocatable::sd(:),lo(:),up(:),cor(:,:),l(:,:),width(:)
    integer,allocatable::ord(:)
    logical::ok
    character(len=256)::msg
    integer::n,i,j
    ctl=probability_control(); if(present(control)) ctl=control
    n=size(lower)
    if(n<1 .or. size(upper)/=n .or. size(mean)/=n .or. size(sigma,1)/=n .or. size(sigma,2)/=n) then
      res%inform=2; res%message='non-conforming dimensions'; return
    end if
    if(any(lower>upper)) then
      res%value=0.0_dp; res%error=0.0_dp; res%message='empty rectangle'; return
    end if
    allocate(sd(n),lo(n),up(n),cor(n,n),width(n),ord(n))
    do i=1,n
      if(sigma(i,i)<=0.0_dp) then; res%inform=3; res%message='non-positive marginal variance'; return; end if
      sd(i)=sqrt(sigma(i,i)); lo(i)=(lower(i)-mean(i))/sd(i); up(i)=(upper(i)-mean(i))/sd(i)
      width(i)=max(0.0_dp,normal_cdf(up(i))-normal_cdf(lo(i))); ord(i)=i
    end do
    if(n==1) then
      res%value=width(1); res%error=2.0_dp*epsilon(1.0_dp); res%evaluations=1; res%message='exact univariate probability'; return
    end if
    call sort_by_width(width,ord)
    lo=lo(ord); up=up(ord); sd=sd(ord)
    do j=1,n
      do i=1,n
        cor(i,j)=sigma(ord(i),ord(j))/(sqrt(sigma(ord(i),ord(i)))*sqrt(sigma(ord(j),ord(j))))
      end do
    end do
    call cholesky_with_jitter(cor,l,ok,msg)
    if(.not.ok) then; res%inform=3; res%message=msg; return; end if
    call qmc_rectangle(lo,up,l,ctl,res)
  end function pmvnorm

  subroutine qmc_rectangle(lo,up,l,ctl,res)
    real(dp),intent(in)::lo(:),up(:),l(:,:)
    type(probability_control),intent(in)::ctl
    type(probability_result),intent(out)::res
    integer::n,nb,per,b,i,j,k,used
    integer(i8)::qstate
    integer,allocatable::primes(:)
    real(dp),allocatable::shift(:),z(:),batch_mean(:)
    real(dp)::s,pl,pu,w,u,pval,m,se,tol
    n=size(lo); nb=max(2,min(ctl%batches,max(2,ctl%maxpts)))
    per=max(1,ctl%maxpts/nb); allocate(primes(max(1,n-1)),shift(max(1,n-1)),z(n),batch_mean(nb))
    call first_primes(max(1,n-1),primes)
    qstate=modulo(int(ctl%seed,kind=i8),2147483646_i8)+1_i8
    used=0
    do b=1,nb
      do i=1,size(shift); shift(i)=local_uniform(qstate); end do
      batch_mean(b)=0.0_dp
      do k=1,per
        pval=1.0_dp; z=0.0_dp
        do i=1,n
          s=0.0_dp
          do j=1,i-1; s=s+l(i,j)*z(j); end do
          pl=normal_cdf((lo(i)-s)/l(i,i)); pu=normal_cdf((up(i)-s)/l(i,i))
          w=max(0.0_dp,pu-pl); pval=pval*w
          if(w<=0.0_dp) exit
          if(i<n) then
            u=modulo(radical_inverse(k,primes(i))+shift(i),1.0_dp)
            u=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),pl+u*w))
            z(i)=normal_quantile(u)
          end if
        end do
        batch_mean(b)=batch_mean(b)+pval
      end do
      batch_mean(b)=batch_mean(b)/real(per,dp); used=used+per
    end do
    m=sum(batch_mean)/real(nb,dp)
    se=sqrt(max(0.0_dp,sum((batch_mean-m)**2)/real(nb-1,dp))/real(nb,dp))
    res%value=min(1.0_dp,max(0.0_dp,m)); res%error=3.0_dp*se; res%evaluations=used
    tol=max(ctl%abseps,ctl%releps*abs(res%value))
    if(res%error<=tol) then; res%inform=0; res%message='randomized Halton-Genz QMC converged'
    else; res%inform=1; res%message='requested tolerance not reached'; end if
  end subroutine qmc_rectangle

  subroutine cholesky_with_jitter(a,l,ok,message)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::l(:,:)
    logical,intent(out)::ok
    character(len=*),intent(out)::message
    real(dp),allocatable::b(:,:)
    real(dp)::jitter
    integer::i,attempt
    b=0.5_dp*(a+transpose(a)); jitter=0.0_dp
    do attempt=1,8
      call cholesky_lower(b,l,ok,message)
      if(ok) return
      if(attempt==1) then; jitter=1.0e-12_dp; else; jitter=10.0_dp*jitter; end if
      b=0.5_dp*(a+transpose(a))
      do i=1,size(a,1); b(i,i)=b(i,i)+jitter; end do
    end do
  end subroutine cholesky_with_jitter

  subroutine sort_by_width(width,ord)
    real(dp),intent(in)::width(:)
    integer,intent(inout)::ord(:)
    integer::i,j,key
    do i=2,size(ord)
      key=ord(i); j=i-1
      do while(j>=1)
        if(width(ord(j))<=width(key)) exit
        ord(j+1)=ord(j); j=j-1
      end do
      ord(j+1)=key
    end do
  end subroutine sort_by_width


  real(dp) function local_uniform(state) result(u)
    integer(i8),intent(inout)::state
    state=modulo(16807_i8*state,2147483647_i8)
    if(state<=0_i8) state=1_i8
    u=real(state,dp)/2147483647.0_dp
  end function local_uniform

  pure real(dp) function radical_inverse(index,base) result(v)
    integer,intent(in)::index,base
    integer::n
    real(dp)::f
    n=index; v=0.0_dp; f=1.0_dp/real(base,dp)
    do while(n>0)
      v=v+f*real(mod(n,base),dp); n=n/base; f=f/real(base,dp)
    end do
  end function radical_inverse

  subroutine first_primes(n,p)
    integer,intent(in)::n
    integer,intent(out)::p(n)
    integer::candidate,k,i
    logical::prime
    if(n<=0) return
    p(1)=2; k=1; candidate=3
    do while(k<n)
      prime=.true.
      do i=1,k
        if(p(i)*p(i)>candidate) exit
        if(mod(candidate,p(i))==0) then; prime=.false.; exit; end if
      end do
      if(prime) then; k=k+1; p(k)=candidate; end if
      candidate=candidate+2
    end do
  end subroutine first_primes

end module anmc_math
