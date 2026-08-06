! SPDX-License-Identifier: GPL-2.0-only
module extra_distr_multivariate
  use extra_distr_kinds, only : dp, pi
  use extra_distr_math
  use extra_distr_rng
  implicit none
  private

  public :: dbvnorm,rbvnorm,dbvpois,rbvpois
  public :: ddirichlet,rdirichlet,ddirmnom,rdirmnom
  public :: dmixnorm,pmixnorm,rmixnorm,dmixpois,pmixpois,rmixpois
  public :: dmnom,rmnom,dmvhyper,rmvhyper

contains

  real(dp) function dbvnorm(x,y,mean1,mean2,sd1,sd2,cor,log_p) result(v)
    real(dp),intent(in)::x,y
    real(dp),intent(in),optional::mean1,mean2,sd1,sd2,cor
    logical,intent(in),optional::log_p
    real(dp)::m1,m2,s1,s2,r,z1,z2,d
    m1=0.0_dp;m2=m1;s1=1.0_dp;s2=s1;r=0.0_dp
    if(present(mean1))m1=mean1;if(present(mean2))m2=mean2
    if(present(sd1))s1=sd1;if(present(sd2))s2=sd2;if(present(cor))r=cor
    if(s1<=0.0_dp.or.s2<=0.0_dp.or.abs(r)>=1.0_dp)then;d=nan_dp();else
      z1=(x-m1)/s1;z2=(y-m2)/s2
      d=exp(-(z1*z1-2.0_dp*r*z1*z2+z2*z2)/(2.0_dp*(1.0_dp-r*r)))/(2.0_dp*pi*s1*s2*sqrt(1.0_dp-r*r))
    end if
    v=apply_density_log(d,log_p)
  end function dbvnorm

  function rbvnorm(n,mean1,mean2,sd1,sd2,cor) result(x)
    integer,intent(in)::n
    real(dp),intent(in),optional::mean1,mean2,sd1,sd2,cor
    real(dp),allocatable::x(:,:)
    real(dp)::m1,m2,s1,s2,r,z1,z2
    integer::i
    m1=0.0_dp;m2=m1;s1=1.0_dp;s2=s1;r=0.0_dp
    if(present(mean1))m1=mean1;if(present(mean2))m2=mean2
    if(present(sd1))s1=sd1;if(present(sd2))s2=sd2;if(present(cor))r=cor
    allocate(x(max(0,n),2))
    do i=1,n
      z1=rnorm_std();z2=rnorm_std()
      x(i,1)=m1+s1*z1
      x(i,2)=m2+s2*(r*z1+sqrt(1.0_dp-r*r)*z2)
    end do
  end function rbvnorm

  real(dp) function dbvpois(x,y,a,b,c,log_p) result(v)
    integer,intent(in)::x,y
    real(dp),intent(in)::a,b,c
    logical,intent(in),optional::log_p
    real(dp)::d
    integer::k
    if(a<0.0_dp.or.b<0.0_dp.or.c<0.0_dp)then
      d=nan_dp()
    else if(x<0.or.y<0)then
      d=0.0_dp
    else
      d=0.0_dp
      do k=0,min(x,y)
        d=d+poisson_pmf(x-k,a)*poisson_pmf(y-k,b)*poisson_pmf(k,c)
      end do
    end if
    v=apply_density_log(d,log_p)
  end function dbvpois

  function rbvpois(n,a,b,c) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::a,b,c
    integer,allocatable::x(:,:)
    integer::i,z
    allocate(x(max(0,n),2))
    do i=1,n
      z=rpoisson_scalar(c)
      x(i,1)=rpoisson_scalar(a)+z
      x(i,2)=rpoisson_scalar(b)+z
    end do
  end function rbvpois

  real(dp) function ddirichlet(x,alpha,log_p,source_compatible) result(v)
    real(dp),intent(in)::x(:),alpha(:)
    logical,intent(in),optional::log_p,source_compatible
    real(dp)::logd,d
    integer::j
    logical::compat
    compat=.true.
    if(present(source_compatible))compat=source_compatible
    if(size(x)/=size(alpha).or.size(x)<2.or.any(alpha<=0.0_dp))then
      d=nan_dp()
    else if(any(x<0.0_dp).or.any(x>1.0_dp))then
      d=0.0_dp
    else if(.not.compat.and.abs(sum(x)-1.0_dp)>1.0e-10_dp)then
      d=0.0_dp
    else
      logd=log_gamma(sum(alpha))-sum(log_gamma(alpha))
      do j=1,size(x)
        if(x(j)==0.0_dp)then
          if(alpha(j)<1.0_dp)then
            logd=pos_inf()
          else
            logd=neg_inf()
          end if
        else
          logd=logd+(alpha(j)-1.0_dp)*log(x(j))
        end if
      end do
      d=exp(logd)
    end if
    v=apply_density_log(d,log_p)
  end function ddirichlet

  function rdirichlet(n,alpha) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::alpha(:)
    real(dp),allocatable::x(:,:)
    real(dp)::s
    integer::i,j
    allocate(x(max(0,n),size(alpha)))
    do i=1,n
      do j=1,size(alpha)
        x(i,j)=rgamma_scalar(alpha(j),1.0_dp)
      end do
      s=sum(x(i,:));x(i,:)=x(i,:)/s
    end do
  end function rdirichlet

  real(dp) function ddirmnom(x,size_value,alpha,log_p) result(v)
    integer,intent(in)::x(:),size_value
    real(dp),intent(in)::alpha(:)
    logical,intent(in),optional::log_p
    real(dp)::logd,d,sa
    integer::j
    if(size(x)/=size(alpha).or.size(x)<2)then
      d=nan_dp()
    else if(size_value<0.or.any(alpha<=0.0_dp))then
      d=nan_dp()
    else if(any(x<0).or.sum(x)/=size_value)then
      d=0.0_dp
    else
      sa=sum(alpha)
      logd=log_factorial(size_value)+log_gamma(sa)-log_gamma(real(size_value,dp)+sa)
      do j=1,size(x)
        logd=logd+log_gamma(real(x(j),dp)+alpha(j))-log_factorial(x(j))-log_gamma(alpha(j))
      end do
      d=exp(logd)
    end if
    v=apply_density_log(d,log_p)
  end function ddirmnom

  function rdirmnom(n,size_value,alpha) result(x)
    integer,intent(in)::n,size_value
    real(dp),intent(in)::alpha(:)
    integer,allocatable::x(:,:)
    real(dp),allocatable::p(:,:)
    integer::i
    p=rdirichlet(n,alpha)
    allocate(x(max(0,n),size(alpha)))
    do i=1,n
      x(i,:)=multinomial_draw(size_value,p(i,:))
    end do
  end function rdirmnom

  real(dp) function dmixnorm(x,mean,sd,alpha,log_p) result(v)
    real(dp),intent(in)::x,mean(:),sd(:),alpha(:)
    logical,intent(in),optional::log_p
    real(dp)::d,total
    integer::j
    total=sum(alpha)
    if(size(mean)/=size(sd).or.size(mean)/=size(alpha).or.size(mean)<1.or. &
       any(sd<=0.0_dp).or.any(alpha<0.0_dp).or.total<=0.0_dp)then
      d=nan_dp()
    else
      d=0.0_dp
      do j=1,size(mean)
        d=d+alpha(j)/total*normal_pdf((x-mean(j))/sd(j))/sd(j)
      end do
    end if
    v=apply_density_log(d,log_p)
  end function dmixnorm

  real(dp) function pmixnorm(q,mean,sd,alpha,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,mean(:),sd(:),alpha(:)
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,total
    integer::j
    total=sum(alpha)
    if(size(mean)/=size(sd).or.size(mean)/=size(alpha).or.size(mean)<1.or. &
       any(sd<=0.0_dp).or.any(alpha<0.0_dp).or.total<=0.0_dp)then
      p=nan_dp()
    else
      p=0.0_dp
      do j=1,size(mean)
        p=p+alpha(j)/total*normal_cdf((q-mean(j))/sd(j))
      end do
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pmixnorm

  function rmixnorm(n,mean,sd,alpha) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sd(:),alpha(:)
    real(dp),allocatable::x(:)
    integer::i,k
    allocate(x(max(0,n)))
    do i=1,n
      k=rcategorical_scalar(alpha)
      x(i)=mean(k)+sd(k)*rnorm_std()
    end do
  end function rmixnorm

  real(dp) function dmixpois(x,lambda,alpha,log_p) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::lambda(:),alpha(:)
    logical,intent(in),optional::log_p
    real(dp)::d,total
    integer::j
    total=sum(alpha)
    if(size(lambda)/=size(alpha).or.size(lambda)<1.or.any(lambda<0.0_dp).or. &
       any(alpha<0.0_dp).or.total<=0.0_dp)then
      d=nan_dp()
    else
      d=0.0_dp
      do j=1,size(lambda)
        d=d+alpha(j)/total*poisson_pmf(x,lambda(j))
      end do
    end if
    v=apply_density_log(d,log_p)
  end function dmixpois

  real(dp) function pmixpois(q,lambda,alpha,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::lambda(:),alpha(:)
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,total
    integer::j
    total=sum(alpha)
    if(size(lambda)/=size(alpha).or.size(lambda)<1.or.any(lambda<0.0_dp).or. &
       any(alpha<0.0_dp).or.total<=0.0_dp)then
      p=nan_dp()
    else
      p=0.0_dp
      do j=1,size(lambda)
        p=p+alpha(j)/total*poisson_cdf(q,lambda(j))
      end do
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pmixpois

  function rmixpois(n,lambda,alpha) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),alpha(:)
    integer,allocatable::x(:)
    integer::i,k
    allocate(x(max(0,n)))
    do i=1,n
      k=rcategorical_scalar(alpha)
      x(i)=rpoisson_scalar(lambda(k))
    end do
  end function rmixpois

  real(dp) function dmnom(x,size_value,prob,log_p) result(v)
    integer,intent(in)::x(:),size_value
    real(dp),intent(in)::prob(:)
    logical,intent(in),optional::log_p
    real(dp)::logd,d,total
    integer::j
    total=sum(prob)
    if(size(x)/=size(prob).or.size(x)<1)then
      d=nan_dp()
    else if(size_value<0.or.any(prob<0.0_dp).or.total<=0.0_dp)then
      d=nan_dp()
    else if(any(x<0).or.sum(x)/=size_value)then
      d=0.0_dp
    else
      logd=log_factorial(size_value)
      do j=1,size(x)
        logd=logd-log_factorial(x(j))
        if(x(j)>0)then
          if(prob(j)<=0.0_dp)then
            logd=neg_inf()
          else
            logd=logd+real(x(j),dp)*log(prob(j)/total)
          end if
        end if
      end do
      d=exp(logd)
    end if
    v=apply_density_log(d,log_p)
  end function dmnom

  function multinomial_draw(size_value,prob) result(x)
    integer,intent(in)::size_value
    real(dp),intent(in)::prob(:)
    integer::x(size(prob))
    integer::j,left
    real(dp)::remaining,pj,total
    x=0;left=size_value;total=sum(prob);remaining=1.0_dp
    do j=1,size(prob)-1
      if(left<=0)exit
      pj=(prob(j)/total)/max(remaining,tiny(1.0_dp))
      pj=min(1.0_dp,max(0.0_dp,pj))
      x(j)=rbinom_scalar(left,pj)
      left=left-x(j);remaining=remaining-prob(j)/total
    end do
    x(size(prob))=left
  end function multinomial_draw

  function rmnom(n,size_value,prob) result(x)
    integer,intent(in)::n,size_value
    real(dp),intent(in)::prob(:)
    integer,allocatable::x(:,:)
    integer::i
    allocate(x(max(0,n),size(prob)))
    do i=1,n
      x(i,:)=multinomial_draw(size_value,prob)
    end do
  end function rmnom

  real(dp) function dmvhyper(x,n_counts,k,log_p) result(v)
    integer,intent(in)::x(:),n_counts(:),k
    logical,intent(in),optional::log_p
    real(dp)::logd,d
    integer::j
    if(size(x)/=size(n_counts).or.size(x)<1)then
      d=nan_dp()
    else if(any(n_counts<0).or.k<0.or.k>sum(n_counts))then
      d=nan_dp()
    else if(any(x<0).or.any(x>n_counts).or.sum(x)/=k)then
      d=0.0_dp
    else
      logd=-log_choose(sum(n_counts),k)
      do j=1,size(x)
        logd=logd+log_choose(n_counts(j),x(j))
      end do
      d=exp(logd)
    end if
    v=apply_density_log(d,log_p)
  end function dmvhyper

  integer function rhyper_scalar(successes,failures,draws) result(x)
    integer,intent(in)::successes,failures,draws
    integer::s,f,d,i
    s=successes;f=failures;x=0
    do i=1,draws
      d=s+f
      if(d<=0)exit
      if(runif_open()<real(s,dp)/real(d,dp))then
        x=x+1;s=s-1
      else
        f=f-1
      end if
    end do
  end function rhyper_scalar

  function rmvhyper(nn,n_counts,k) result(x)
    integer,intent(in)::nn,n_counts(:),k
    integer,allocatable::x(:,:)
    integer::i,j,left,remaining
    allocate(x(max(0,nn),size(n_counts)))
    do i=1,nn
      left=k;remaining=sum(n_counts)
      do j=1,size(n_counts)-1
        x(i,j)=rhyper_scalar(n_counts(j),remaining-n_counts(j),left)
        left=left-x(i,j);remaining=remaining-n_counts(j)
      end do
      x(i,size(n_counts))=left
    end do
  end function rmvhyper

end module extra_distr_multivariate
