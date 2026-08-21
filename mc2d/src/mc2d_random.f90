! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_random
  use mc2d_kinds, only : dp, pi, nan_dp
  implicit none
  private
  public :: seed_random, random_normal, random_gamma, random_beta, random_poisson
  public :: random_binomial, random_uniform_open
contains
  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n,i
    integer, allocatable :: put(:)
    call random_seed(size=n); allocate(put(n))
    do i=1,n
      put(i)=modulo(seed+104729*i+37*i*i,huge(1)-1)+1
    end do
    call random_seed(put=put)
  end subroutine seed_random

  real(dp) function random_uniform_open() result(u)
    call random_number(u)
    u=min(1.0_dp-epsilon(1.0_dp),max(tiny(1.0_dp),u))
  end function random_uniform_open

  real(dp) function random_normal() result(z)
    real(dp) :: u1,u2
    u1=random_uniform_open(); u2=random_uniform_open()
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function random_normal

  recursive real(dp) function random_gamma(shape,scale) result(x)
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    real(dp) :: d,c,z,u,v,s
    s=1.0_dp; if(present(scale)) s=scale
    if(shape<=0.0_dp .or. s<=0.0_dp) then; x=nan_dp(); return; end if
    if(shape<1.0_dp) then
      u=random_uniform_open(); x=s*random_gamma(shape+1.0_dp)*u**(1.0_dp/shape); return
    end if
    d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z=random_normal(); v=1.0_dp+c*z
        if(v>0.0_dp) exit
      end do
      v=v*v*v; u=random_uniform_open()
      if(u<1.0_dp-0.0331_dp*z**4) exit
      if(log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x=s*d*v
  end function random_gamma

  real(dp) function random_beta(a,b) result(x)
    real(dp),intent(in)::a,b
    real(dp)::g1,g2
    g1=random_gamma(a); g2=random_gamma(b); x=g1/(g1+g2)
  end function random_beta

  integer function random_poisson(lambda) result(k)
    real(dp), intent(in) :: lambda
    real(dp) :: l,p,u,y,em,t,sq,alxm,g
    integer :: kk
    if(lambda<0.0_dp) then; k=-1; return; end if
    if(lambda==0.0_dp) then; k=0; return; end if
    if(lambda<30.0_dp) then
      l=exp(-lambda); p=1.0_dp; k=-1
      do
        k=k+1; p=p*random_uniform_open()
        if(p<=l) exit
      end do
      return
    end if
    sq=sqrt(2.0_dp*lambda); alxm=log(lambda); g=lambda*alxm-log_gamma(lambda+1.0_dp)
    do
      do
        y=tan(pi*random_uniform_open()); em=sq*y+lambda
        if(em>=0.0_dp) exit
      end do
      em=floor(em); t=0.9_dp*(1.0_dp+y*y)*exp(em*alxm-log_gamma(em+1.0_dp)-g)
      u=random_uniform_open()
      if(u<=t) exit
    end do
    kk=int(em); k=kk
  end function random_poisson

  integer function random_binomial(n,p) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::p
    integer::i
    real(dp)::u
    x=0
    if(n<=0) return
    do i=1,n
      call random_number(u); if(u<p) x=x+1
    end do
  end function random_binomial
end module mc2d_random
