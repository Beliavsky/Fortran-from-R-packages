! SPDX-License-Identifier: Artistic-2.0
module lamp_process
  use ecd_kinds, only : dp, pi, ecd_ok, ecd_invalid
  use ecd_rng, only : rng_state, rng_uniform, rng_normal, rng_exponential
  use ecld_models, only : ecld_new, ecld_sd
  implicit none
  private

  type, public :: lamp_model
    real(dp) :: lambda=4.0_dp, alpha=0.5_dp, beta=1.0_dp
    integer :: parameterization=1, random_walk=1, sd_method=0
    real(dp) :: target_sd=-1.0_dp, t_infinity=86400000.0_dp
    integer :: random_count=1000000
    real(dp) :: n_lower=0.0_dp, n_upper=1000.0_dp
  end type lamp_model

  type, public :: lamp_result
    real(dp), allocatable :: z(:),b(:),n(:),tau(:)
    integer :: used_tau=0, status=ecd_ok
  end type lamp_result

  public :: lamp_new, lamp_sd_factor, stable_random
  public :: lamp_generate_tau, lamp_stable_random_walk, lamp_simulate_once, lamp_simulate

contains

  function lamp_new(lambda,alpha,beta,random_walk,t_infinity,random_count,target_sd,sd_method,n_lower,n_upper,status) result(m)
    real(dp), intent(in), optional :: lambda,alpha,beta,t_infinity,target_sd,n_lower,n_upper
    integer, intent(in), optional :: random_walk,random_count,sd_method
    integer, intent(out), optional :: status
    type(lamp_model) :: m
    if(present(lambda))m%lambda=lambda
    if(present(alpha))m%lambda=2.0_dp/alpha
    m%alpha=2.0_dp/m%lambda
    if(present(beta))m%beta=beta
    if(present(random_walk))m%random_walk=random_walk
    if(present(t_infinity))m%t_infinity=t_infinity
    if(present(random_count))m%random_count=random_count
    if(present(target_sd))m%target_sd=target_sd
    if(present(sd_method))m%sd_method=sd_method
    if(present(n_lower))m%n_lower=n_lower
    if(present(n_upper))m%n_upper=n_upper
    if(present(status))status=ecd_ok
    if(m%lambda<=0.0_dp .or. abs(m%beta)>1.0_dp .or. m%t_infinity<=0.0_dp .or. &
       m%random_count<1 .or. m%n_lower<0.0_dp .or. m%n_upper<m%n_lower .or. &
       all(m%random_walk/=[1,11,2,22])) then
      if(present(status))status=ecd_invalid
    end if
  end function lamp_new

  function lamp_sd_factor(m) result(f)
    type(lamp_model), intent(in) :: m
    real(dp) :: f,sd0
    if(m%target_sd<0.0_dp) then; f=1.0_dp; return; end if
    sd0=ecld_sd(ecld_new(lambda=m%lambda))
    if(m%random_walk==1 .or. m%random_walk==11) then
      f=m%target_sd/sd0
    else
      f=(m%target_sd/sd0)**(1.0_dp/(1.0_dp+m%alpha/2.0_dp))
    end if
  end function lamp_sd_factor

  function stable_random(rng,alpha,beta,scale,location) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: alpha,beta
    real(dp), intent(in), optional :: scale,location
    real(dp) :: x,s,mu,v,w,tanpa,b,fac
    s=1.0_dp; mu=0.0_dp
    if(present(scale))s=scale
    if(present(location))mu=location
    v=pi*(rng_uniform(rng)-0.5_dp); w=rng_exponential(rng)
    if(abs(alpha-1.0_dp)>1.0e-10_dp) then
      tanpa=tan(0.5_dp*pi*alpha)
      b=atan(beta*tanpa)/alpha
      fac=(1.0_dp+(beta*tanpa)**2)**(0.5_dp/alpha)
      x=fac*sin(alpha*(v+b))/cos(v)**(1.0_dp/alpha)* &
        (cos(v-alpha*(v+b))/w)**((1.0_dp-alpha)/alpha)
    else
      x=(2.0_dp/pi)*((0.5_dp*pi+beta*v)*tan(v)- &
        beta*log((0.5_dp*pi*w*cos(v))/(0.5_dp*pi+beta*v)))
    end if
    x=mu+s*x
  end function stable_random

  subroutine lamp_generate_tau(m,rng,tau,status)
    type(lamp_model), intent(in) :: m
    type(rng_state), intent(inout) :: rng
    real(dp), allocatable, intent(out) :: tau(:)
    integer, intent(out), optional :: status
    real(dp) :: g,scale,sd_factor,signv
    integer :: i
    if(present(status))status=ecd_ok
    if(m%beta==1.0_dp .and. m%lambda<=2.0_dp) then
      allocate(tau(0)); if(present(status))status=ecd_invalid; return
    end if
    allocate(tau(m%random_count))
    sd_factor=merge(lamp_sd_factor(m),1.0_dp,m%sd_method==0)
    if(m%alpha<1.0_dp) then; g=cos(pi*m%alpha/2.0_dp)**(1.0_dp/m%alpha)
    else if(abs(m%alpha-1.0_dp)<1e-12_dp) then; g=0.125_dp
    else; g=1.0_dp; end if
    scale=g/sd_factor
    do i=1,size(tau)
      tau(i)=stable_random(rng,m%alpha,m%beta,scale)
      if(m%beta==1.0_dp) then
        signv=merge(-1.0_dp,1.0_dp,rng_uniform(rng)<0.5_dp)
        tau(i)=tau(i)*signv
      end if
    end do
  end subroutine lamp_generate_tau

  function lamp_stable_random_walk(m,rng,n,b) result(x)
    type(lamp_model), intent(in) :: m
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: n,b
    real(dp) :: x,epsilon_v,scale
    if(n<=0) then; x=0.0_dp; return; end if
    epsilon_v=m%t_infinity**(-1.0_dp/m%lambda)
    select case(m%random_walk)
    case(1)
      x=merge(-1.0_dp,1.0_dp,rng_uniform(rng)<0.5_dp)*rng_exponential(rng)
    case(2)
      x=epsilon_v*sqrt(real(n,dp))*rng_normal(rng)
    case(11)
      scale=sqrt(2.0_dp*rng_exponential(rng))
      x=real(b,dp)/sqrt(real(n,dp))*scale
    case(22)
      x=real(b,dp)*epsilon_v
    case default
      x=0.0_dp
    end select
  end function lamp_stable_random_walk

  subroutine lamp_simulate_once(m,rng,result,drop,keep_tau)
    type(lamp_model), intent(in) :: m
    type(rng_state), intent(inout) :: rng
    type(lamp_result), intent(out) :: result
    integer, intent(in), optional :: drop,keep_tau
    real(dp), allocatable :: tau(:),zt(:),bt(:),nt(:)
    real(dp) :: t,tn,cnt,tinf2,sdf
    integer :: i,d,keep,n,b,nz,maxz
    d=10; keep=1
    if(present(drop))d=drop
    if(present(keep_tau))keep=keep_tau
    call lamp_generate_tau(m,rng,tau,result%status)
    if(result%status/=ecd_ok)return
    sdf=merge(lamp_sd_factor(m),1.0_dp,m%sd_method==1)
    tinf2=m%t_infinity*sdf
    maxz=max(1,size(tau)); allocate(zt(maxz),bt(maxz),nt(maxz))
    tn=0.0_dp; n=0; b=0; nz=0; i=1
    do while(i<=size(tau)-max(0,d))
      t=tau(i); i=i+1
      if(abs(t)>tinf2)cycle
      tn=tn+abs(t); n=n+1; b=b+merge(-1,1,t<0.0_dp)
      do while(tn>=tinf2)
        n=n-1; b=b-merge(-1,1,t<0.0_dp)
        cnt=merge(real(n,dp)**(m%lambda/2.0_dp)/m%t_infinity,0.0_dp,n>0)
        if(n>0 .and. cnt>=m%n_lower .and. cnt<=m%n_upper) then
          nz=nz+1; nt(nz)=cnt; bt(nz)=lamp_stable_random_walk(m,rng,n,b); zt(nz)=bt(nz)*cnt
        end if
        tn=tn-tinf2; b=merge(-1,1,t<0.0_dp); n=merge(1,0,tn>0.0_dp)
      end do
    end do
    allocate(result%z(nz),result%b(nz),result%n(nz))
    result%z=zt(:nz); result%b=bt(:nz); result%n=nt(:nz); result%used_tau=i-1
    select case(keep)
    case(0); allocate(result%tau(0))
    case(1); allocate(result%tau(max(0,size(tau)-i+1))); if(size(result%tau)>0)result%tau=tau(i:)
    case default; allocate(result%tau(size(tau))); result%tau=tau
    end select
  end subroutine lamp_simulate_once

  subroutine lamp_simulate(m,rng,target_length,result,drop)
    type(lamp_model), intent(in) :: m
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: target_length
    type(lamp_result), intent(out) :: result
    integer, intent(in), optional :: drop
    type(lamp_result) :: one
    real(dp), allocatable :: z(:),b(:),n(:)
    integer :: used,take,old
    allocate(z(0),b(0),n(0)); used=0
    do while(size(z)<target_length)
      call lamp_simulate_once(m,rng,one,drop,0)
      if(one%status/=ecd_ok) then; result%status=one%status; return; end if
      take=min(size(one%z),target_length-size(z))
      old=size(z)
      z=[z,one%z(:take)]; b=[b,one%b(:take)]; n=[n,one%n(:take)]
      used=used+one%used_tau
      if(take==0 .and. old==size(z))exit
    end do
    call move_alloc(z,result%z); call move_alloc(b,result%b); call move_alloc(n,result%n)
    allocate(result%tau(0)); result%used_tau=used
  end subroutine lamp_simulate

end module lamp_process
