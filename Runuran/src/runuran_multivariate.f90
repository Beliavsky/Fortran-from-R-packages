! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran_multivariate
  use runuran_kinds, only : dp, pi, huge_dp
  use runuran_rng, only : rng_state, rng_uniform, rng_normal, rng_chisq, rng_exponential
  implicit none
  private
  integer, parameter :: MV_CUSTOM=0,MV_NORMAL=1,MV_T=2,MV_CAUCHY=3,MV_EXP=4

  abstract interface
    function mv_pdf_callback(x,params) result(y)
      import dp
      real(dp),intent(in)::x(:),params(:)
      real(dp)::y
    end function
  end interface

  type, public :: multivariate_distribution
    integer :: id=MV_CUSTOM, dim=1
    character(len=32)::name='custom'
    real(dp),allocatable::mean(:),cov(:,:),chol(:,:),params(:),lower(:),upper(:)
    procedure(mv_pdf_callback),pointer,nopass::pdf_cb=>null(),logpdf_cb=>null()
  contains
    procedure :: pdf=>mv_pdf
    procedure :: logpdf=>mv_logpdf
    procedure :: sample=>mv_sample
  end type

  type, public :: mv_generator
    integer :: method=1
    type(multivariate_distribution)::distr
    real(dp),allocatable::state(:)
    integer::n_trials=0,n_accept=0
  contains
    procedure::sample=>mvg_sample
    procedure::sample_n=>mvg_sample_n
    procedure::acceptance_rate=>mvg_acceptance_rate
  end type

  public :: udmultinormal, udmultistudent, udmulticauchy, udmultiexponential, udmultivariate
  public :: hitro_new, vnrou_new, gibbs_new
contains
  subroutine cholesky_spd(a,l,ok)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(size(a,1),size(a,2))
    logical,intent(out)::ok
    integer::i,j,k,n
    real(dp)::s
    n=size(a,1)
    l=0.0_dp
    ok=.true.
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1
        s=s-l(i,k)*l(j,k)
        end do
        if(i==j)then
          if(s<=0.0_dp)then
          ok=.false.
          return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine

  real(dp) function logdet_chol(l) result(v)
    real(dp),intent(in)::l(:,:)
    integer::i
    v=0.0_dp
    do i=1,size(l,1)
    v=v+2.0_dp*log(l(i,i))
    end do
  end function

  subroutine solve_lower(l,b,x)
    real(dp),intent(in)::l(:,:),b(:)
    real(dp),intent(out)::x(size(b))
    integer::i,j
    do i=1,size(b)
    x(i)=b(i)
    do j=1,i-1
    x(i)=x(i)-l(i,j)*x(j)
    end do
    x(i)=x(i)/l(i,i)
    end do
  end subroutine

  function udmultinormal(mean,cov,lower,upper) result(d)
    real(dp),intent(in)::mean(:),cov(:,:)
    real(dp),intent(in),optional::lower(:),upper(:)
    type(multivariate_distribution)::d
    logical::ok
    d%id=MV_NORMAL
    d%name='multinormal'
    d%dim=size(mean)
    allocate(d%mean(d%dim),d%cov(d%dim,d%dim),d%chol(d%dim,d%dim),d%lower(d%dim),d%upper(d%dim))
    d%mean=mean
    d%cov=cov
    d%lower=-huge_dp
    d%upper=huge_dp
    if(present(lower))d%lower=lower
    if(present(upper))d%upper=upper
    call cholesky_spd(cov,d%chol,ok)
    if(.not.ok)error stop 'udmultinormal: covariance is not SPD'
  end function
  function udmultistudent(df,mean,cov,lower,upper) result(d)
    real(dp),intent(in)::df,mean(:),cov(:,:)
    real(dp),intent(in),optional::lower(:),upper(:)
    type(multivariate_distribution)::d
    d=udmultinormal(mean,cov,lower,upper)
    d%id=MV_T
    d%name='multistudent'
    allocate(d%params(1))
    d%params=[df]
  end function
  function udmulticauchy(mean,cov,lower,upper) result(d)
    real(dp),intent(in)::mean(:),cov(:,:)
    real(dp),intent(in),optional::lower(:),upper(:)
    type(multivariate_distribution)::d
    d=udmultistudent(1.0_dp,mean,cov,lower,upper)
    d%id=MV_CAUCHY
    d%name='multicauchy'
  end function
  function udmultiexponential(dim) result(d)
    integer,intent(in)::dim
    type(multivariate_distribution)::d
    d%id=MV_EXP
    d%name='multiexponential'
    d%dim=dim
    allocate(d%mean(dim),d%lower(dim),d%upper(dim))
    d%mean=1.0_dp
    d%lower=0.0_dp
    d%upper=huge_dp
  end function
  function udmultivariate(dim,pdf,logpdf,params,lower,upper) result(d)
    integer,intent(in)::dim
    procedure(mv_pdf_callback)::pdf
    procedure(mv_pdf_callback),optional::logpdf
    real(dp),intent(in),optional::params(:),lower(:),upper(:)
    type(multivariate_distribution)::d
    d%id=MV_CUSTOM
    d%dim=dim
    d%pdf_cb=>pdf
    allocate(d%mean(dim),d%lower(dim),d%upper(dim))
    d%mean=0.0_dp
    d%lower=-huge_dp
    d%upper=huge_dp
    if(present(logpdf))d%logpdf_cb=>logpdf
    if(present(params))then
    allocate(d%params(size(params)))
    d%params=params
    else
    allocate(d%params(0))
    end if
    if(present(lower))d%lower=lower
    if(present(upper))d%upper=upper
  end function

  real(dp) function mv_logpdf(self,x) result(lp)
    class(multivariate_distribution),intent(in)::self
    real(dp),intent(in)::x(:)
    real(dp),allocatable::y(:)
    real(dp)::q,df
    if(any(x<self%lower).or.any(x>self%upper))then
    lp=-huge_dp
    return
    end if
    select case(self%id)
    case(MV_CUSTOM)
      if(associated(self%logpdf_cb))then
      lp=self%logpdf_cb(x,self%params)
      else
      lp=log(max(self%pdf_cb(x,self%params),tiny(1.0_dp)))
      end if
    case(MV_NORMAL)
      allocate(y(self%dim))
      call solve_lower(self%chol,x-self%mean,y)
      q=sum(y*y)
      lp=-0.5_dp*(real(self%dim,dp)*log(2.0_dp*pi)+logdet_chol(self%chol)+q)
    case(MV_T,MV_CAUCHY)
      df=self%params(1)
      allocate(y(self%dim))
      call solve_lower(self%chol,x-self%mean,y)
      q=sum(y*y)
      lp=log_gamma(0.5_dp*(df+self%dim))-log_gamma(0.5_dp*df) &
        -0.5_dp*(real(self%dim,dp)*log(df*pi)+logdet_chol(self%chol)) &
        -0.5_dp*(df+self%dim)*log(1.0_dp+q/df)
    case(MV_EXP)
      if(any(x<0.0_dp))then
      lp=-huge_dp
      else
      lp=-sum(x)
      end if
    case default
    lp=-huge_dp
    end select
  end function
  real(dp) function mv_pdf(self,x) result(p)
  class(multivariate_distribution),intent(in)::self
  real(dp),intent(in)::x(:)
  real(dp)::lp
  lp=self%logpdf(x)
  if(lp<log(tiny(1.0_dp)))then
  p=0.0_dp
  else
  p=exp(lp)
  end if
  end function

  subroutine mv_sample(self,rng,x)
    class(multivariate_distribution),intent(in)::self
    type(rng_state),intent(inout)::rng
    real(dp),intent(out)::x(:)
    real(dp),allocatable::z(:)
    real(dp)::s
    integer::i
    if(size(x)/=self%dim)error stop 'mv_sample: wrong output dimension'
    select case(self%id)
    case(MV_NORMAL,MV_T,MV_CAUCHY)
      allocate(z(self%dim))
      do
      i=0
      z=0.0_dp
      do i=1,self%dim
      z(i)=rng_normal(rng)
      end do
        x=self%mean+matmul(self%chol,z)
        if(self%id==MV_T.or.self%id==MV_CAUCHY)then
        s=sqrt(rng_chisq(rng,self%params(1))/self%params(1))
        x=self%mean+(x-self%mean)/s
        end if
        if(all(x>=self%lower).and.all(x<=self%upper))exit
      end do
    case(MV_EXP)
      do i=1,self%dim
      x(i)=rng_exponential(rng)
      end do
    case default
      error stop 'mv_sample: custom distributions require HITRO/Gibbs generator'
    end select
  end subroutine

  function hitro_new(distr,start) result(g)
    type(multivariate_distribution),intent(in)::distr
    real(dp),intent(in),optional::start(:)
    type(mv_generator)::g
    g%method=1
    g%distr=distr
    allocate(g%state(distr%dim))
    g%state=distr%mean
    if(present(start))g%state=start
  end function
  function vnrou_new(distr,start) result(g)
    type(multivariate_distribution),intent(in)::distr
    real(dp),intent(in),optional::start(:)
    type(mv_generator)::g
    g=hitro_new(distr,start)
    g%method=2
  end function
  function gibbs_new(distr,start) result(g)
    type(multivariate_distribution),intent(in)::distr
    real(dp),intent(in),optional::start(:)
    type(mv_generator)::g
    g=hitro_new(distr,start)
    g%method=3
  end function

  subroutine mvg_sample(self,rng,x)
    class(mv_generator),intent(inout)::self
    type(rng_state),intent(inout)::rng
    real(dp),intent(out)::x(:)
    real(dp),allocatable::dir(:),prop(:)
    real(dp)::normd,step,lp0,lp1
    integer::i
    if(self%distr%id/=MV_CUSTOM)then
    call self%distr%sample(rng,x)
    self%state=x
    self%n_trials=self%n_trials+1
    self%n_accept=self%n_accept+1
    return
    end if
    allocate(dir(self%distr%dim),prop(self%distr%dim))
    lp0=self%distr%logpdf(self%state)
    ! Hit-and-run random-walk Metropolis kernel. It preserves the requested
    ! multivariate target and is the native fallback for HITRO/VNROU on
    ! arbitrary user callbacks.
    do i=1,self%distr%dim
    dir(i)=rng_normal(rng)
    end do
    normd=sqrt(sum(dir*dir))
    dir=dir/normd
    step=rng_normal(rng)
    prop=self%state+step*dir
    lp1=self%distr%logpdf(prop)
    self%n_trials=self%n_trials+1
    if(log(rng_uniform(rng))<=lp1-lp0)then
    self%state=prop
    self%n_accept=self%n_accept+1
    end if
    x=self%state
  end subroutine
  subroutine mvg_sample_n(self,rng,x)
    class(mv_generator),intent(inout)::self
    type(rng_state),intent(inout)::rng
    real(dp),intent(out)::x(:,:)
    integer::i
    if(size(x,1)/=self%distr%dim)error stop 'mvg_sample_n: first dimension mismatch'
    do i=1,size(x,2)
    call self%sample(rng,x(:,i))
    end do
  end subroutine
  real(dp) function mvg_acceptance_rate(self) result(r)
  class(mv_generator),intent(in)::self
  if(self%n_trials==0)then
  r=0.0_dp
  else
  r=real(self%n_accept,dp)/real(self%n_trials,dp)
  end if
  end function
end module runuran_multivariate
