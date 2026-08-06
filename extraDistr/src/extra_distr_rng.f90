! SPDX-License-Identifier: GPL-2.0-only
module extra_distr_rng
  use extra_distr_kinds, only : dp, pi
  use extra_distr_math, only : poisson_quantile, binom_quantile, nan_dp
  implicit none
  private
  public :: seed_rng, runif_open, rnorm_std, rgamma_scalar, rbeta_scalar
  public :: rchisq_scalar, rt_scalar, rpoisson_scalar, rbinom_scalar
  public :: rnbinom_scalar, rcategorical_scalar

contains

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=modulo(seed+104729*i+37*i*i, huge(1)-1)+1
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  real(dp) function runif_open() result(u)
    do
      call random_number(u)
      if(u>0.0_dp .and. u<1.0_dp) exit
    end do
  end function runif_open

  real(dp) function rnorm_std() result(z)
    real(dp) :: u1,u2
    u1=runif_open(); u2=runif_open()
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rnorm_std

  recursive real(dp) function rgamma_scalar(shape,scale) result(x)
    real(dp), intent(in)::shape,scale
    real(dp)::d,c,z,u
    if(shape<=0.0_dp .or. scale<=0.0_dp) then
      x=nan_dp()
    else if(shape<1.0_dp) then
      x=rgamma_scalar(shape+1.0_dp,scale)*runif_open()**(1.0_dp/shape)
    else
      d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
      do
        do
          z=rnorm_std()
          if(1.0_dp+c*z>0.0_dp) exit
        end do
        x=d*(1.0_dp+c*z)**3
        u=runif_open()
        if(u<1.0_dp-0.0331_dp*z**4) exit
        if(log(u)<0.5_dp*z*z+d*(1.0_dp-x/d+log(x/d))) exit
      end do
      x=scale*x
    end if
  end function rgamma_scalar

  real(dp) function rbeta_scalar(a,b) result(x)
    real(dp), intent(in)::a,b
    real(dp)::g1,g2
    g1=rgamma_scalar(a,1.0_dp); g2=rgamma_scalar(b,1.0_dp)
    x=g1/(g1+g2)
  end function rbeta_scalar

  real(dp) function rchisq_scalar(df) result(x)
    real(dp), intent(in)::df
    x=rgamma_scalar(0.5_dp*df,2.0_dp)
  end function rchisq_scalar

  real(dp) function rt_scalar(df) result(x)
    real(dp), intent(in)::df
    x=rnorm_std()/sqrt(rchisq_scalar(df)/df)
  end function rt_scalar

  integer function rpoisson_scalar(lambda) result(k)
    real(dp), intent(in)::lambda
    real(dp)::l,p,u,z
    if(lambda<30.0_dp) then
      l=exp(-lambda); p=1.0_dp; k=0
      do
        k=k+1; p=p*runif_open()
        if(p<=l) exit
      end do
      k=k-1
    else
      ! Inversion using a uniform is reliable for all lambda and preserves discreteness.
      u=runif_open(); k=poisson_quantile(u,lambda)
      if(k==huge(1)) then
        z=rnorm_std(); k=max(0,nint(lambda+sqrt(lambda)*z))
      end if
    end if
  end function rpoisson_scalar

  integer function rbinom_scalar(n,p) result(k)
    integer, intent(in)::n
    real(dp), intent(in)::p
    integer :: i
    if(n<=200) then
      k=0
      do i=1,n
        if(runif_open()<p) k=k+1
      end do
    else
      k=binom_quantile(runif_open(),n,p)
    end if
  end function rbinom_scalar

  integer function rnbinom_scalar(size,p) result(k)
    real(dp), intent(in)::size,p
    real(dp)::lambda
    lambda=rgamma_scalar(size,(1.0_dp-p)/p)
    k=rpoisson_scalar(lambda)
  end function rnbinom_scalar

  integer function rcategorical_scalar(prob) result(k)
    real(dp), intent(in)::prob(:)
    integer :: i
    real(dp)::u,s,total
    total=sum(max(prob,0.0_dp)); u=runif_open()*total; s=0.0_dp
    k=size(prob)
    do i=1,size(prob)
      s=s+max(prob(i),0.0_dp)
      if(u<=s) then; k=i; return; end if
    end do
  end function rcategorical_scalar

end module extra_distr_rng
