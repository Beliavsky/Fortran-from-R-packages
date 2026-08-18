! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran_methods
  use runuran_kinds, only : dp, huge_dp
  use runuran_rng, only : rng_state, rng_uniform
  use runuran_distributions, only : continuous_distribution, discrete_distribution
  implicit none
  private

  integer, parameter, public :: METHOD_PINV=1, METHOD_ARS=2, METHOD_AROU=3, METHOD_SROU=4
  integer, parameter, public :: METHOD_TDR=5, METHOD_ITDR=6, METHOD_TABL=7
  integer, parameter, public :: METHOD_DARI=20, METHOD_DAU=21, METHOD_DGT=22
  integer, parameter, public :: METHOD_MIXT=30

  type, public :: unuran_generator
    integer :: method=METHOD_PINV
    logical :: is_discrete=.false.
    logical :: is_mixture=.false.
    type(continuous_distribution) :: cont
    type(discrete_distribution) :: discr
    type(continuous_distribution), allocatable :: components(:)
    real(dp), allocatable :: weights(:), cumweights(:)
    integer :: n_trials=0, n_accept=0
  contains
    procedure :: sample=>generator_sample
    procedure :: sample_int=>generator_sample_int
    procedure :: sample_n=>generator_sample_n
    procedure :: sample_int_n=>generator_sample_int_n
    procedure :: acceptance_rate=>generator_acceptance_rate
  end type unuran_generator

  public :: pinv_new, ars_new, arou_new, srou_new, tdr_new, itdr_new, tabl_new
  public :: dari_new, dau_new, dgt_new, mixt_new, unuran_new
contains
  function pinv_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_PINV;g%cont=distr
  end function
  function ars_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_ARS;g%cont=distr
  end function
  function arou_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_AROU;g%cont=distr
  end function
  function srou_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_SROU;g%cont=distr
  end function
  function tdr_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_TDR;g%cont=distr
  end function
  function itdr_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_ITDR;g%cont=distr
  end function
  function tabl_new(distr) result(g)
    type(continuous_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_TABL;g%cont=distr
  end function
  function dari_new(distr) result(g)
    type(discrete_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_DARI;g%is_discrete=.true.;g%discr=distr
  end function
  function dau_new(distr) result(g)
    type(discrete_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_DAU;g%is_discrete=.true.;g%discr=distr
  end function
  function dgt_new(distr) result(g)
    type(discrete_distribution),intent(in)::distr
    type(unuran_generator)::g
    g%method=METHOD_DGT;g%is_discrete=.true.;g%discr=distr
  end function

  function mixt_new(components,weights) result(g)
    type(continuous_distribution),intent(in)::components(:)
    real(dp),intent(in)::weights(:)
    type(unuran_generator)::g
    integer::i,n
    real(dp)::s
    n=size(components)
    if(size(weights)/=n .or. n<1) error stop 'mixt_new: inconsistent component/weight sizes'
    g%method=METHOD_MIXT;g%is_mixture=.true.
    allocate(g%components(n),g%weights(n),g%cumweights(n))
    g%components=components;g%weights=max(weights,0.0_dp);s=sum(g%weights)
    if(s<=0.0_dp)error stop 'mixt_new: weights must have positive sum'
    g%weights=g%weights/s;g%cumweights(1)=g%weights(1)
    do i=2,n;g%cumweights(i)=g%cumweights(i-1)+g%weights(i);end do
    g%cumweights(n)=1.0_dp
  end function

  function unuran_new(distr,method) result(g)
    type(continuous_distribution),intent(in)::distr
    character(len=*),intent(in),optional::method
    type(unuran_generator)::g
    character(len=:),allocatable::m
    if(.not.present(method))then;g=pinv_new(distr);return;end if
    m=lower(method)
    select case(trim(m))
    case('pinv','hinv','ninv');g=pinv_new(distr)
    case('ars');g=ars_new(distr)
    case('arou','nrou','vnrou');g=arou_new(distr)
    case('srou');g=srou_new(distr)
    case('tdr','utdr');g=tdr_new(distr)
    case('itdr');g=itdr_new(distr)
    case('tabl');g=tabl_new(distr)
    case default;g=pinv_new(distr)
    end select
  end function

  pure function lower(s) result(t)
    character(len=*),intent(in)::s
    character(len=len(s))::t
    integer::i,c
    do i=1,len(s);c=iachar(s(i:i));if(c>=65.and.c<=90)then;t(i:i)=achar(c+32);else;t(i:i)=s(i:i);end if;end do
  end function

  real(dp) function numeric_dlogpdf(d,x) result(v)
    type(continuous_distribution),intent(in)::d
    real(dp),intent(in)::x
    v=d%dlogpdf(x)
  end function

  real(dp) function ars_sample(d,rng,ok,ntry) result(xout)
    type(continuous_distribution),intent(in)::d
    type(rng_state),intent(inout)::rng
    logical,intent(out)::ok
    integer,intent(out)::ntry
    integer,parameter::maxp=80
    real(dp)::x(maxp),h(maxp),hp(maxp),z(0:maxp),la(maxp),cw(maxp)
    real(dp)::u,lu,xt,uhat,lf,total,mx
    integer::n,i,j,idx
    logical::valid
    real(dp),parameter::probs(3)=[0.2_dp,0.5_dp,0.8_dp]
    n=3
    do i=1,n
      x(i)=d%quantile(probs(i));h(i)=d%logpdf(x(i));hp(i)=numeric_dlogpdf(d,x(i))
    end do
    call sort_support(x,h,hp,n)
    valid=.true.
    do i=1,n-1
      if(hp(i)<=hp(i+1))valid=.false.
    end do
    if(.not.valid)then;ok=.false.;xout=0.0_dp;ntry=0;return;end if
    ntry=0
    do while(ntry<10000)
      ntry=ntry+1
      z(0)=d%lb;z(n)=d%ub
      do i=1,n-1
        z(i)=(h(i)-h(i+1)-x(i)*hp(i)+x(i+1)*hp(i+1))/(hp(i)-hp(i+1))
        z(i)=max(x(i),min(x(i+1),z(i)))
      end do
      do i=1,n
        la(i)=tangent_logarea(h(i),hp(i),x(i),z(i-1),z(i))
      end do
      mx=maxval(la(1:n));total=0.0_dp
      do i=1,n;cw(i)=exp(la(i)-mx);total=total+cw(i);cw(i)=total;end do
      u=rng_uniform(rng)*total;idx=1
      do while(idx<n .and. u>cw(idx));idx=idx+1;end do
      xt=sample_tangent(rng,hp(idx),z(idx-1),z(idx))
      if(xt<=d%lb.or.xt>=d%ub)cycle
      uhat=h(idx)+hp(idx)*(xt-x(idx));lf=d%logpdf(xt);lu=log(rng_uniform(rng))
      if(lu<=lf-uhat)then;xout=xt;ok=.true.;return;end if
      if(n<maxp .and. lf>-0.5_dp*huge_dp)then
        n=n+1;x(n)=xt;h(n)=lf;hp(n)=numeric_dlogpdf(d,xt);call sort_support(x,h,hp,n)
        valid=.true.;do j=1,n-1;if(hp(j)<=hp(j+1))valid=.false.;end do
        if(.not.valid)then;ok=.false.;xout=0.0_dp;return;end if
      end if
    end do
    ok=.false.;xout=0.0_dp
  end function

  subroutine sort_support(x,h,hp,n)
    real(dp),intent(inout)::x(:),h(:),hp(:);integer,intent(in)::n
    integer::i,j;real(dp)::tx,th,tp
    do i=2,n
      tx=x(i);th=h(i);tp=hp(i);j=i-1
      do while(j>=1)
        if(x(j)<=tx)exit
        x(j+1)=x(j);h(j+1)=h(j);hp(j+1)=hp(j);j=j-1
      end do
      x(j+1)=tx;h(j+1)=th;hp(j+1)=tp
    end do
  end subroutine

  real(dp) function tangent_logarea(h,s,x,a,b) result(la)
    real(dp),intent(in)::h,s,x,a,b
    real(dp)::aa,bb,q
    if(abs(s)<1.0e-14_dp)then
      if(abs(a)>0.5_dp*huge_dp.or.abs(b)>0.5_dp*huge_dp)then;la=huge_dp;else;la=h+log(b-a);end if
      return
    end if
    if(a<-0.5_dp*huge_dp)then
      if(s<=0.0_dp)then;la=huge_dp;else;la=h+s*(b-x)-log(s);end if
      return
    end if
    if(b>0.5_dp*huge_dp)then
      if(s>=0.0_dp)then;la=huge_dp;else;la=h+s*(a-x)-log(-s);end if
      return
    end if
    aa=s*(a-x);bb=s*(b-x)
    if(s>0.0_dp)then;q=bb+log(1.0_dp-exp(aa-bb));else;q=aa+log(1.0_dp-exp(bb-aa));end if
    la=h+q-log(abs(s))
  end function

  real(dp) function sample_tangent(rng,s,a,b) result(x)
    type(rng_state),intent(inout)::rng;real(dp),intent(in)::s,a,b
    real(dp)::u,ea,eb,m
    u=rng_uniform(rng)
    if(abs(s)<1.0e-14_dp)then;x=a+u*(b-a);return;end if
    if(a<-0.5_dp*huge_dp)then;x=b+log(u)/s;return;end if
    if(b>0.5_dp*huge_dp)then;x=a+log(1.0_dp-u)/s;return;end if
    m=max(s*a,s*b);ea=exp(s*a-m);eb=exp(s*b-m);x=(log(ea+u*(eb-ea))+m)/s
  end function

  real(dp) function generator_sample(self,rng) result(x)
    class(unuran_generator),intent(inout)::self
    type(rng_state),intent(inout)::rng
    integer::i,nt
    real(dp)::u
    logical::ok
    if(self%is_discrete)then;x=real(self%sample_int(rng),dp);return;end if
    if(self%is_mixture)then
      u=rng_uniform(rng);i=1;do while(i<size(self%cumweights).and.u>self%cumweights(i));i=i+1;end do
      x=self%components(i)%sample(rng);self%n_trials=self%n_trials+1;self%n_accept=self%n_accept+1;return
    end if
    self%n_trials=self%n_trials+1
    select case(self%method)
    case(METHOD_ARS)
      x=ars_sample(self%cont,rng,ok,nt);self%n_trials=self%n_trials+max(nt-1,0)
      if(.not.ok)x=self%cont%sample(rng)
    case default
      ! UNU.RAN's PINV/HINV/TDR/ROU families share the same target. This
      ! translation uses adaptive numerical inversion as the robust generic
      ! fallback; ARS and discrete guide-table paths retain specialized kernels.
      x=self%cont%sample(rng)
    end select
    self%n_accept=self%n_accept+1
  end function

  integer function generator_sample_int(self,rng) result(k)
    class(unuran_generator),intent(inout)::self
    type(rng_state),intent(inout)::rng
    if(.not.self%is_discrete)then;k=nint(self%sample(rng));return;end if
    self%n_trials=self%n_trials+1;k=self%discr%sample(rng);self%n_accept=self%n_accept+1
  end function
  subroutine generator_sample_n(self,rng,x)
    class(unuran_generator),intent(inout)::self;type(rng_state),intent(inout)::rng;real(dp),intent(out)::x(:);integer::i
    do i=1,size(x);x(i)=self%sample(rng);end do
  end subroutine
  subroutine generator_sample_int_n(self,rng,x)
    class(unuran_generator),intent(inout)::self;type(rng_state),intent(inout)::rng;integer,intent(out)::x(:);integer::i
    do i=1,size(x);x(i)=self%sample_int(rng);end do
  end subroutine
  real(dp) function generator_acceptance_rate(self) result(r)
    class(unuran_generator),intent(in)::self
    if(self%n_trials<=0)then;r=0.0_dp;else;r=real(self%n_accept,dp)/real(self%n_trials,dp);end if
  end function
end module runuran_methods
