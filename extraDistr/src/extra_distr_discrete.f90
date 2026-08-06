! SPDX-License-Identifier: GPL-2.0-only
module extra_distr_discrete
  use extra_distr_kinds, only : dp
  use extra_distr_math
  use extra_distr_rng
  implicit none
  private

  public :: dbern,pbern,qbern,rbern
  public :: dbbinom,pbbinom,rbbinom,dbnbinom,pbnbinom,rbnbinom
  public :: dcat,pcat,qcat,rcat,rcatlp
  public :: ddgamma,pdgamma,rdgamma,ddlaplace,pdlaplace,rdlaplace
  public :: ddnorm,pdnorm,rdnorm,ddunif,pdunif,qdunif,rdunif
  public :: ddweibull,pdweibull,qdweibull,rdweibull
  public :: dgpois,pgpois,rgpois,dlgser,plgser,qlgser,rlgser
  public :: dnhyper,pnhyper,qnhyper,rnhyper
  public :: dskellam,rskellam,rsign
  public :: dtbinom,ptbinom,qtbinom,rtbinom
  public :: dtpois,ptpois,qtpois,rtpois
  public :: dzib,pzib,qzib,rzib,dzinb,pzinb,qzinb,rzinb,dzip,pzip,qzip,rzip

contains

  real(dp) function dbern(x,prob,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in),optional::prob;logical,intent(in),optional::log_p;real(dp)::p,d
    p=0.5_dp;if(present(prob))p=prob
    if(p<0.0_dp.or.p>1.0_dp)then;d=nan_dp();else if(x==0)then;d=1.0_dp-p;else if(x==1)then;d=p;else;d=0.0_dp;end if
    v=apply_density_log(d,log_p)
  end function dbern
  real(dp) function pbern(q,prob,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in),optional::prob
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,pr
    pr=0.5_dp;if(present(prob))pr=prob
    if(pr<0.0_dp.or.pr>1.0_dp)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else if(q<1)then
      p=1.0_dp-pr
    else
      p=1.0_dp
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pbern
  integer function qbern(probability,prob,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability;real(dp),intent(in),optional::prob;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,pr
    p=decode_probability(probability,lower_tail,log_p);pr=0.5_dp;if(present(prob))pr=prob;x=merge(0,1,p<=1.0_dp-pr)
  end function qbern
  function rbern(n,prob) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::prob;integer,allocatable::x(:);real(dp)::p;integer::i
    p=0.5_dp;if(present(prob))p=prob;allocate(x(max(0,n)));do i=1,n;x(i)=merge(1,0,runif_open()<p);end do
  end function rbern

  real(dp) function dbbinom(x,size,alpha,beta,log_p) result(v)
    integer,intent(in)::x,size;real(dp),intent(in),optional::alpha,beta;logical,intent(in),optional::log_p;real(dp)::a,b,d
    a=1.0_dp;b=1.0_dp;if(present(alpha))a=alpha;if(present(beta))b=beta
    if(a<=0.0_dp.or.b<=0.0_dp.or.size<0)then;d=nan_dp();else if(x<0.or.x>size)then;d=0.0_dp;else;d=exp(log_choose(size,x)+log_beta(real(x,dp)+a,real(size-x,dp)+b)-log_beta(a,b));end if;v=apply_density_log(d,log_p)
  end function dbbinom
  real(dp) function pbbinom(q,size,alpha,beta,lower_tail,log_p) result(v)
    integer,intent(in)::q,size
    real(dp),intent(in),optional::alpha,beta
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,a,b
    integer::k
    a=1.0_dp;b=1.0_dp
    if(present(alpha))a=alpha
    if(present(beta))b=beta
    if(a<=0.0_dp.or.b<=0.0_dp.or.size<0)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else if(q>=size)then
      p=1.0_dp
    else
      p=0.0_dp
      do k=0,q
        p=p+dbbinom(k,size,a,b)
      end do
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pbbinom
  function rbbinom(n,size,alpha,beta) result(x)
    integer,intent(in)::n,size;real(dp),intent(in),optional::alpha,beta;integer,allocatable::x(:);real(dp)::a,b,p;integer::i
    a=1.0_dp;b=1.0_dp;if(present(alpha))a=alpha;if(present(beta))b=beta;allocate(x(max(0,n)));do i=1,n;p=rbeta_scalar(a,b);x(i)=rbinom_scalar(size,p);end do
  end function rbbinom

  real(dp) function dbnbinom(x,size,alpha,beta,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in)::size;real(dp),intent(in),optional::alpha,beta;logical,intent(in),optional::log_p;real(dp)::a,b,d
    a=1.0_dp;b=1.0_dp;if(present(alpha))a=alpha;if(present(beta))b=beta
    if(a<=0.0_dp.or.b<=0.0_dp.or.size<=0.0_dp)then;d=nan_dp();else if(x<0)then;d=0.0_dp;else;d=exp(log_gamma(size+real(x,dp))-log_gamma(real(x+1,dp))-log_gamma(size)+log_beta(a+size,b+real(x,dp))-log_beta(a,b));end if;v=apply_density_log(d,log_p)
  end function dbnbinom
  real(dp) function pbnbinom(q,size,alpha,beta,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::size
    real(dp),intent(in),optional::alpha,beta
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,a,b
    integer::k
    a=1.0_dp;b=1.0_dp
    if(present(alpha))a=alpha
    if(present(beta))b=beta
    if(a<=0.0_dp.or.b<=0.0_dp.or.size<=0.0_dp)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else
      p=0.0_dp
      do k=0,q
        p=p+dbnbinom(k,size,a,b)
      end do
      p=min(p,1.0_dp)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pbnbinom
  function rbnbinom(n,size,alpha,beta) result(x)
    integer,intent(in)::n;real(dp),intent(in)::size;real(dp),intent(in),optional::alpha,beta;integer,allocatable::x(:);real(dp)::a,b,p;integer::i
    a=1.0_dp;b=1.0_dp;if(present(alpha))a=alpha;if(present(beta))b=beta;allocate(x(max(0,n)));do i=1,n;p=rbeta_scalar(a,b);x(i)=rnbinom_scalar(size,p);end do
  end function rbnbinom

  real(dp) function dcat(x,prob,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in)::prob(:);logical,intent(in),optional::log_p;real(dp)::d,total
    total=sum(prob);if(any(prob<0.0_dp).or.total<=0.0_dp)then;d=nan_dp();else if(x<1.or.x>size(prob))then;d=0.0_dp;else;d=prob(x)/total;end if;v=apply_density_log(d,log_p)
  end function dcat
  real(dp) function pcat(q,prob,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::prob(:)
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p,total
    total=sum(prob)
    if(any(prob<0.0_dp).or.total<=0.0_dp)then
      p=nan_dp()
    else if(q<1)then
      p=0.0_dp
    else if(q>=size(prob))then
      p=1.0_dp
    else
      p=sum(prob(1:q))/total
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pcat
  integer function qcat(probability,prob,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability,prob(:);logical,intent(in),optional::lower_tail,log_p;real(dp)::p,c,total;integer::i
    p=decode_probability(probability,lower_tail,log_p);total=sum(prob);c=0.0_dp;x=size(prob);do i=1,size(prob);c=c+prob(i)/total;if(p<=c)then;x=i;return;end if;end do
  end function qcat
  function rcat(n,prob) result(x)
    integer,intent(in)::n;real(dp),intent(in)::prob(:);integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=rcategorical_scalar(prob);end do
  end function rcat
  function rcatlp(n,log_prob) result(x)
    integer,intent(in)::n;real(dp),intent(in)::log_prob(:);integer,allocatable::x(:);real(dp),allocatable::p(:);real(dp)::m
    m=maxval(log_prob);allocate(p(size(log_prob)));p=exp(log_prob-m);x=rcat(n,p)
  end function rcatlp

  real(dp) function ddgamma(x,shape,rate,scale,log_p) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::log_p
    real(dp)::r,s,d
    r=1.0_dp;if(present(rate))r=rate
    if(r>0.0_dp)then;s=1.0_dp/r;else;s=nan_dp();end if
    if(present(scale))s=scale
    if(shape<=0.0_dp.or.s<=0.0_dp)then
      d=nan_dp()
    else if(x<0)then
      d=0.0_dp
    else
      d=regularized_gamma_p(shape,real(x+1,dp)/s)-regularized_gamma_p(shape,real(x,dp)/s)
    end if
    v=apply_density_log(d,log_p)
  end function ddgamma
  real(dp) function pdgamma(q,shape,rate,scale,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::r,s,pr
    r=1.0_dp;if(present(rate))r=rate
    if(r>0.0_dp)then;s=1.0_dp/r;else;s=nan_dp();end if
    if(present(scale))s=scale
    if(shape<=0.0_dp.or.s<=0.0_dp)then
      pr=nan_dp()
    else if(q<0)then
      pr=0.0_dp
    else
      pr=regularized_gamma_p(shape,real(q+1,dp)/s)
    end if
    v=apply_tail(pr,lower_tail,log_p)
  end function pdgamma
  function rdgamma(n,shape,rate,scale) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::rate,scale
    integer,allocatable::x(:)
    real(dp)::r,s
    integer::i
    r=1.0_dp;if(present(rate))r=rate
    s=1.0_dp/r;if(present(scale))s=scale
    allocate(x(max(0,n)))
    do i=1,n
      x(i)=floor(rgamma_scalar(shape,s))
    end do
  end function rdgamma

  real(dp) function ddlaplace(x,location,scale,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in)::location,scale;logical,intent(in),optional::log_p;real(dp)::d
    if(scale<=0.0_dp.or.scale>=1.0_dp)then;d=nan_dp();else;d=(1.0_dp-scale)/(1.0_dp+scale)*scale**abs(real(x,dp)-location);end if;v=apply_density_log(d,log_p)
  end function ddlaplace
  real(dp) function pdlaplace(q,location,scale,lower_tail,log_p,source_compatible) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::location,scale
    logical,intent(in),optional::lower_tail,log_p,source_compatible
    real(dp)::p,z
    logical::compat,left_branch
    z=real(q,dp)-location;compat=.true.
    if(present(source_compatible))compat=source_compatible
    left_branch=merge(q<0,z<0.0_dp,compat)
    if(scale<=0.0_dp.or.scale>=1.0_dp)then
      p=nan_dp()
    else if(left_branch)then
      p=scale**(-floor(z))/(1.0_dp+scale)
    else
      p=1.0_dp-scale**(floor(z)+1.0_dp)/(1.0_dp+scale)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pdlaplace
  function rdlaplace(n,location,scale) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::location,scale
    real(dp),allocatable::x(:)
    integer::i,u,w
    allocate(x(max(0,n)))
    do i=1,n
      u=int(log(runif_open())/log(scale))
      w=int(log(runif_open())/log(scale))
      x(i)=location+real(u-w,dp)
    end do
  end function rdlaplace

  real(dp) function ddnorm(x,mean,sd,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in),optional::mean,sd;logical,intent(in),optional::log_p;real(dp)::m,s,d
    m=0.0_dp;s=1.0_dp;if(present(mean))m=mean;if(present(sd))s=sd;if(s<=0.0_dp)then;d=nan_dp();else;d=normal_cdf((real(x+1,dp)-m)/s)-normal_cdf((real(x,dp)-m)/s);end if;v=apply_density_log(d,log_p)
  end function ddnorm
  real(dp) function pdnorm(q,mean,sd,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in),optional::mean,sd
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::m,s,p
    m=0.0_dp;s=1.0_dp
    if(present(mean))m=mean
    if(present(sd))s=sd
    if(s<=0.0_dp)then
      p=nan_dp()
    else
      p=normal_cdf((real(q+1,dp)-m)/s)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pdnorm
  function rdnorm(n,mean,sd) result(x)
    integer,intent(in)::n;real(dp),intent(in),optional::mean,sd;integer,allocatable::x(:);real(dp)::m,s;integer::i
    m=0.0_dp;s=1.0_dp;if(present(mean))m=mean;if(present(sd))s=sd;allocate(x(max(0,n)));do i=1,n;x(i)=floor(m+s*rnorm_std());end do
  end function rdnorm

  real(dp) function ddunif(x,min_value,max_value,log_p) result(v)
    integer,intent(in)::x,min_value,max_value;logical,intent(in),optional::log_p;real(dp)::d
    if(min_value>max_value)then;d=nan_dp();else if(x<min_value.or.x>max_value)then;d=0.0_dp;else;d=1.0_dp/real(max_value-min_value+1,dp);end if;v=apply_density_log(d,log_p)
  end function ddunif
  real(dp) function pdunif(q,min_value,max_value,lower_tail,log_p) result(v)
    integer,intent(in)::q,min_value,max_value
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    if(max_value<min_value)then
      p=nan_dp()
    else if(q<min_value)then
      p=0.0_dp
    else if(q>=max_value)then
      p=1.0_dp
    else
      p=real(q-min_value+1,dp)/real(max_value-min_value+1,dp)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pdunif
  integer function qdunif(prob,min_value,max_value,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;integer,intent(in)::min_value,max_value;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(prob,lower_tail,log_p);if(p<=0.0_dp)then;x=min_value;else;x=ceiling(p*real(max_value-min_value+1,dp)+real(min_value-1,dp));end if
  end function qdunif
  function rdunif(n,min_value,max_value) result(x)
    integer,intent(in)::n,min_value,max_value;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=min_value+int(runif_open()*real(max_value-min_value+1,dp));end do
  end function rdunif

  real(dp) function ddweibull(x,shape1,shape2,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in)::shape1,shape2;logical,intent(in),optional::log_p;real(dp)::d
    if(shape1<=0.0_dp.or.shape1>=1.0_dp.or.shape2<=0.0_dp)then;d=nan_dp();else if(x<0)then;d=0.0_dp;else;d=shape1**(real(x,dp)**shape2)-shape1**(real(x+1,dp)**shape2);end if;v=apply_density_log(d,log_p)
  end function ddweibull
  real(dp) function pdweibull(q,shape1,shape2,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::shape1,shape2
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    if(shape1<=0.0_dp.or.shape1>=1.0_dp.or.shape2<=0.0_dp)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else
      p=1.0_dp-shape1**(real(q+1,dp)**shape2)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pdweibull
  integer function qdweibull(prob,shape1,shape2,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,shape1,shape2;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(prob,lower_tail,log_p);if(p<=0.0_dp)then;x=0;else;x=ceiling((log(1.0_dp-p)/log(shape1))**(1.0_dp/shape2)-1.0_dp);end if
  end function qdweibull
  function rdweibull(n,shape1,shape2) result(x)
    integer,intent(in)::n;real(dp),intent(in)::shape1,shape2;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=ceiling((log(runif_open())/log(shape1))**(1.0_dp/shape2)-1.0_dp);end do
  end function rdweibull

  real(dp) function dgpois(x,shape,rate,scale,log_p) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::log_p
    real(dp)::r,s,d,p
    r=1.0_dp;if(present(rate))r=rate
    if(r>0.0_dp)then;s=1.0_dp/r;else;s=nan_dp();end if
    if(present(scale))s=scale
    if(shape<=0.0_dp.or.s<=0.0_dp)then
      d=nan_dp()
    else
      p=1.0_dp/(1.0_dp+s)
      d=nbinom_pmf(x,shape,p)
    end if
    v=apply_density_log(d,log_p)
  end function dgpois
  real(dp) function pgpois(q,shape,rate,scale,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::rate,scale
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::r,s,p,pr
    r=1.0_dp;if(present(rate))r=rate
    if(r>0.0_dp)then;s=1.0_dp/r;else;s=nan_dp();end if
    if(present(scale))s=scale
    if(shape<=0.0_dp.or.s<=0.0_dp)then
      pr=nan_dp()
    else
      p=1.0_dp/(1.0_dp+s)
      pr=nbinom_cdf(q,shape,p)
    end if
    v=apply_tail(pr,lower_tail,log_p)
  end function pgpois
  function rgpois(n,shape,rate,scale) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::rate,scale
    integer,allocatable::x(:)
    real(dp)::r,s,p
    integer::i
    r=1.0_dp;if(present(rate))r=rate
    s=1.0_dp/r;if(present(scale))s=scale
    p=1.0_dp/(1.0_dp+s)
    allocate(x(max(0,n)))
    do i=1,n
      x(i)=rnbinom_scalar(shape,p)
    end do
  end function rgpois

  real(dp) function dlgser(x,theta,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in)::theta;logical,intent(in),optional::log_p;real(dp)::d
    if(theta<=0.0_dp.or.theta>=1.0_dp)then;d=nan_dp();else if(x<1)then;d=0.0_dp;else;d=-theta**x/(real(x,dp)*log(1.0_dp-theta));end if;v=apply_density_log(d,log_p)
  end function dlgser
  real(dp) function plgser(q,theta,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::theta
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    integer::k
    if(theta<=0.0_dp.or.theta>=1.0_dp)then
      p=nan_dp()
    else if(q<1)then
      p=0.0_dp
    else
      p=0.0_dp
      do k=1,q
        p=p+dlgser(k,theta)
      end do
      p=min(p,1.0_dp)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function plgser
  integer function qlgser(prob,theta,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob,theta;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,pk
    p=decode_probability(prob,lower_tail,log_p);x=1;pk=dlgser(1,theta);do while(p>pk);p=p-pk;pk=pk*theta*real(x,dp)/real(x+1,dp);x=x+1;end do
  end function qlgser
  function rlgser(n,theta) result(x)
    integer,intent(in)::n;real(dp),intent(in)::theta;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=qlgser(runif_open(),theta);end do
  end function rlgser

  real(dp) function dnhyper(x,n,m,r,log_p) result(v)
    integer,intent(in)::x,n,m,r;logical,intent(in),optional::log_p;real(dp)::d
    if(n<0.or.m<0.or.r<0.or.r>m)then;d=nan_dp();else if(x<r.or.x>n+r)then;d=0.0_dp;else;d=exp(log_choose(x-1,r-1)+log_choose(m+n-x,m-r)-log_choose(m+n,m));end if;v=apply_density_log(d,log_p)
  end function dnhyper
  real(dp) function pnhyper(q,n,m,r,lower_tail,log_p) result(v)
    integer,intent(in)::q,n,m,r
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    integer::k
    if(n<0.or.m<0.or.r<0.or.r>m)then
      p=nan_dp()
    else if(q<r)then
      p=0.0_dp
    else if(q>=n+r)then
      p=1.0_dp
    else
      p=0.0_dp
      do k=r,q
        p=p+dnhyper(k,n,m,r)
      end do
      p=min(p,1.0_dp)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pnhyper
  integer function qnhyper(prob,n,m,r,lower_tail,log_p) result(x)
    real(dp),intent(in)::prob;integer,intent(in)::n,m,r;logical,intent(in),optional::lower_tail,log_p;real(dp)::p,c
    p=decode_probability(prob,lower_tail,log_p);c=0.0_dp;do x=r,n+r;c=c+dnhyper(x,n,m,r);if(c>=p)return;end do
  end function qnhyper
  function rnhyper(nn,n,m,r) result(x)
    integer,intent(in)::nn,n,m,r;integer,allocatable::x(:);integer::i;allocate(x(max(0,nn)));do i=1,nn;x(i)=qnhyper(runif_open(),n,m,r);end do
  end function rnhyper

  real(dp) function dskellam(x,mu1,mu2,log_p) result(v)
    integer,intent(in)::x;real(dp),intent(in)::mu1,mu2;logical,intent(in),optional::log_p;real(dp)::d,z
    if(mu1<0.0_dp.or.mu2<0.0_dp)then;d=nan_dp();else if(mu1==0.0_dp.and.mu2==0.0_dp)then;d=merge(1.0_dp,0.0_dp,x==0);else if(mu1==0.0_dp)then;d=poisson_pmf(-x,mu2);else if(mu2==0.0_dp)then;d=poisson_pmf(x,mu1);else;z=2.0_dp*sqrt(mu1*mu2);d=exp(-(mu1+mu2))*(mu1/mu2)**(0.5_dp*real(x,dp))*bessel_i_integer(abs(x),z);end if;v=apply_density_log(d,log_p)
  end function dskellam
  function rskellam(n,mu1,mu2) result(x)
    integer,intent(in)::n;real(dp),intent(in)::mu1,mu2;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=rpoisson_scalar(mu1)-rpoisson_scalar(mu2);end do
  end function rskellam
  function rsign(n) result(x)
    integer,intent(in)::n;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;x(i)=merge(1,-1,runif_open()>0.5_dp);end do
  end function rsign

  real(dp) function dtbinom(x,size,prob,a,b,log_p) result(v)
    integer,intent(in)::x,size
    real(dp),intent(in)::prob
    integer,intent(in),optional::a,b
    logical,intent(in),optional::log_p
    integer::lo,hi
    real(dp)::d,den
    lo=0;hi=size
    if(present(a))lo=max(0,a)
    if(present(b))hi=min(size,b)
    if(size<0.or.prob<0.0_dp.or.prob>1.0_dp.or.lo>hi)then
      d=nan_dp()
    else
      den=binom_cdf(hi,size,prob)-binom_cdf(lo-1,size,prob)
      if(den<=0.0_dp)then
        d=nan_dp()
      else if(x<lo.or.x>hi)then
        d=0.0_dp
      else
        d=binom_pmf(x,size,prob)/den
      end if
    end if
    v=apply_density_log(d,log_p)
  end function dtbinom
  real(dp) function ptbinom(q,size,prob,a,b,lower_tail,log_p) result(v)
    integer,intent(in)::q,size
    real(dp),intent(in)::prob
    integer,intent(in),optional::a,b
    logical,intent(in),optional::lower_tail,log_p
    integer::lo,hi
    real(dp)::p,den
    lo=0;hi=size
    if(present(a))lo=max(0,a)
    if(present(b))hi=min(size,b)
    if(size<0.or.prob<0.0_dp.or.prob>1.0_dp.or.lo>hi)then
      p=nan_dp()
    else
      den=binom_cdf(hi,size,prob)-binom_cdf(lo-1,size,prob)
      if(den<=0.0_dp)then
        p=nan_dp()
      else if(q<lo)then
        p=0.0_dp
      else if(q>=hi)then
        p=1.0_dp
      else
        p=(binom_cdf(q,size,prob)-binom_cdf(lo-1,size,prob))/den
      end if
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function ptbinom
  integer function qtbinom(probability,size,prob,a,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability,prob;integer,intent(in)::size;integer,intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;integer::lo,hi;real(dp)::p,c
    p=decode_probability(probability,lower_tail,log_p);lo=0;hi=size;if(present(a))lo=max(0,a);if(present(b))hi=min(size,b);c=0.0_dp;do x=lo,hi;c=c+dtbinom(x,size,prob,lo,hi);if(c>=p)return;end do
  end function qtbinom
  function rtbinom(n,size,prob,a,b) result(x)
    integer,intent(in)::n,size;real(dp),intent(in)::prob;integer,intent(in),optional::a,b;integer,allocatable::x(:);integer::lo,hi,i
    lo=0;hi=size;if(present(a))lo=max(0,a);if(present(b))hi=min(size,b);allocate(x(max(0,n)));do i=1,n;x(i)=qtbinom(runif_open(),size,prob,lo,hi);end do
  end function rtbinom

  real(dp) function dtpois(x,lambda,a,b,log_p) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::lambda
    integer,intent(in),optional::a,b
    logical,intent(in),optional::log_p
    integer::lo,hi
    real(dp)::d,den
    lo=0;hi=1000
    if(lambda>=0.0_dp)hi=max(1000,poisson_quantile(1.0_dp-1.0e-14_dp,lambda))
    if(present(a))lo=max(0,a)
    if(present(b))hi=b
    if(lambda<0.0_dp.or.lo>hi)then
      d=nan_dp()
    else
      den=poisson_cdf(hi,lambda)-poisson_cdf(lo-1,lambda)
      if(den<=0.0_dp)then
        d=nan_dp()
      else if(x<lo.or.x>hi)then
        d=0.0_dp
      else
        d=poisson_pmf(x,lambda)/den
      end if
    end if
    v=apply_density_log(d,log_p)
  end function dtpois
  real(dp) function ptpois(q,lambda,a,b,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::lambda
    integer,intent(in),optional::a,b
    logical,intent(in),optional::lower_tail,log_p
    integer::lo,hi
    real(dp)::p,den
    lo=0;hi=1000
    if(lambda>=0.0_dp)hi=max(1000,poisson_quantile(1.0_dp-1.0e-14_dp,lambda))
    if(present(a))lo=max(0,a)
    if(present(b))hi=b
    if(lambda<0.0_dp.or.lo>hi)then
      p=nan_dp()
    else
      den=poisson_cdf(hi,lambda)-poisson_cdf(lo-1,lambda)
      if(den<=0.0_dp)then
        p=nan_dp()
      else if(q<lo)then
        p=0.0_dp
      else if(q>=hi)then
        p=1.0_dp
      else
        p=(poisson_cdf(q,lambda)-poisson_cdf(lo-1,lambda))/den
      end if
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function ptpois
  integer function qtpois(probability,lambda,a,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability,lambda;integer,intent(in),optional::a,b;logical,intent(in),optional::lower_tail,log_p;integer::lo,hi;real(dp)::p,c
    p=decode_probability(probability,lower_tail,log_p);lo=0;hi=max(1000,poisson_quantile(1.0_dp-1.0e-14_dp,lambda));if(present(a))lo=max(0,a);if(present(b))hi=b;c=0.0_dp;do x=lo,hi;c=c+dtpois(x,lambda,lo,hi);if(c>=p)return;end do
  end function qtpois
  function rtpois(n,lambda,a,b) result(x)
    integer,intent(in)::n;real(dp),intent(in)::lambda;integer,intent(in),optional::a,b;integer,allocatable::x(:);integer::lo,hi,i
    lo=0;hi=max(1000,poisson_quantile(1.0_dp-1.0e-14_dp,lambda));if(present(a))lo=max(0,a);if(present(b))hi=b;allocate(x(max(0,n)));do i=1,n;x(i)=qtpois(runif_open(),lambda,lo,hi);end do
  end function rtpois

  real(dp) function dzib(x,size,prob,pi0,log_p) result(v)
    integer,intent(in)::x,size
    real(dp),intent(in)::prob,pi0
    logical,intent(in),optional::log_p
    real(dp)::d
    if(size<0.or.prob<0.0_dp.or.prob>1.0_dp.or.pi0<0.0_dp.or.pi0>1.0_dp)then
      d=nan_dp()
    else if(x==0)then
      d=pi0+(1.0_dp-pi0)*binom_pmf(0,size,prob)
    else
      d=(1.0_dp-pi0)*binom_pmf(x,size,prob)
    end if
    v=apply_density_log(d,log_p)
  end function dzib
  real(dp) function pzib(q,size,prob,pi0,lower_tail,log_p) result(v)
    integer,intent(in)::q,size
    real(dp),intent(in)::prob,pi0
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    if(size<0.or.prob<0.0_dp.or.prob>1.0_dp.or.pi0<0.0_dp.or.pi0>1.0_dp)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else
      p=pi0+(1.0_dp-pi0)*binom_cdf(q,size,prob)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pzib
  integer function qzib(probability,size,prob,pi0,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability,prob,pi0;integer,intent(in)::size;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(probability,lower_tail,log_p);if(p<pi0)then;x=0;else;x=binom_quantile((p-pi0)/(1.0_dp-pi0),size,prob);end if
  end function qzib
  function rzib(n,size,prob,pi0) result(x)
    integer,intent(in)::n,size;real(dp),intent(in)::prob,pi0;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;if(runif_open()<pi0)then;x(i)=0;else;x(i)=rbinom_scalar(size,prob);end if;end do
  end function rzib

  real(dp) function dzinb(x,size,prob,pi0,log_p) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::size,prob,pi0
    logical,intent(in),optional::log_p
    real(dp)::d
    if(size<=0.0_dp.or.prob<=0.0_dp.or.prob>1.0_dp.or.pi0<0.0_dp.or.pi0>1.0_dp)then
      d=nan_dp()
    else if(x==0)then
      d=pi0+(1.0_dp-pi0)*nbinom_pmf(0,size,prob)
    else
      d=(1.0_dp-pi0)*nbinom_pmf(x,size,prob)
    end if
    v=apply_density_log(d,log_p)
  end function dzinb
  real(dp) function pzinb(q,size,prob,pi0,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::size,prob,pi0
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    if(size<=0.0_dp.or.prob<=0.0_dp.or.prob>1.0_dp.or.pi0<0.0_dp.or.pi0>1.0_dp)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else
      p=pi0+(1.0_dp-pi0)*nbinom_cdf(q,size,prob)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pzinb
  integer function qzinb(probability,size,prob,pi0,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability,size,prob,pi0;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(probability,lower_tail,log_p);if(p<pi0)then;x=0;else;x=nbinom_quantile((p-pi0)/(1.0_dp-pi0),size,prob);end if
  end function qzinb
  function rzinb(n,size,prob,pi0) result(x)
    integer,intent(in)::n;real(dp),intent(in)::size,prob,pi0;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;if(runif_open()<pi0)then;x(i)=0;else;x(i)=rnbinom_scalar(size,prob);end if;end do
  end function rzinb

  real(dp) function dzip(x,lambda,pi0,log_p) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::lambda,pi0
    logical,intent(in),optional::log_p
    real(dp)::d
    if(lambda<0.0_dp.or.pi0<0.0_dp.or.pi0>1.0_dp)then
      d=nan_dp()
    else if(x==0)then
      d=pi0+(1.0_dp-pi0)*poisson_pmf(0,lambda)
    else
      d=(1.0_dp-pi0)*poisson_pmf(x,lambda)
    end if
    v=apply_density_log(d,log_p)
  end function dzip
  real(dp) function pzip(q,lambda,pi0,lower_tail,log_p) result(v)
    integer,intent(in)::q
    real(dp),intent(in)::lambda,pi0
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::p
    if(lambda<0.0_dp.or.pi0<0.0_dp.or.pi0>1.0_dp)then
      p=nan_dp()
    else if(q<0)then
      p=0.0_dp
    else
      p=pi0+(1.0_dp-pi0)*poisson_cdf(q,lambda)
    end if
    v=apply_tail(p,lower_tail,log_p)
  end function pzip
  integer function qzip(probability,lambda,pi0,lower_tail,log_p) result(x)
    real(dp),intent(in)::probability,lambda,pi0;logical,intent(in),optional::lower_tail,log_p;real(dp)::p
    p=decode_probability(probability,lower_tail,log_p);if(p<pi0)then;x=0;else;x=poisson_quantile((p-pi0)/(1.0_dp-pi0),lambda);end if
  end function qzip
  function rzip(n,lambda,pi0) result(x)
    integer,intent(in)::n;real(dp),intent(in)::lambda,pi0;integer,allocatable::x(:);integer::i;allocate(x(max(0,n)));do i=1,n;if(runif_open()<pi0)then;x(i)=0;else;x(i)=rpoisson_scalar(lambda);end if;end do
  end function rzip

end module extra_distr_discrete
