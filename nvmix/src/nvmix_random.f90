! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_random
  use nvmix_kinds, only : dp, i8, pi
  implicit none
  private
  integer(i8) :: state=88172645463325252_i8
  public :: seed_random, random_uniform, random_normal, random_gamma, halton
contains
  subroutine seed_random(seed)
    integer(i8), intent(in) :: seed
    state=ieor(abs(seed)+104729_i8,4101842887655102017_i8)
    if(state==0_i8)state=88172645463325252_i8
  end subroutine
  real(dp) function random_uniform() result(u)
    integer(i8) :: x
    x=state
    x=ieor(x,shiftl(x,13)); x=ieor(x,shiftr(x,7)); x=ieor(x,shiftl(x,17))
    state=x
    u=real(iand(x,int(z'7FFFFFFFFFFFFFFF',i8)),dp)/real(huge(1_i8),dp)
    u=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u))
  end function
  real(dp) function random_normal() result(z)
    real(dp) :: u1,u2
    u1=random_uniform(); u2=random_uniform()
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function
  recursive real(dp) function random_gamma(shape,scale) result(x)
    real(dp), intent(in) :: shape,scale
    real(dp) :: d,c,z,u,v
    if(shape<=0.0_dp .or. scale<=0.0_dp) then; x=0.0_dp; return; end if
    if(shape<1.0_dp) then
      x=random_gamma(shape+1.0_dp,scale)*random_uniform()**(1.0_dp/shape)
      return
    end if
    d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
    do
      z=random_normal(); v=(1.0_dp+c*z)**3
      if(v<=0.0_dp)cycle
      u=random_uniform()
      if(u<1.0_dp-0.0331_dp*z**4)exit
      if(log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v)))exit
    end do
    x=scale*d*v
  end function
  pure real(dp) function halton(index,base,shift) result(u)
    integer, intent(in) :: index,base
    real(dp), intent(in), optional :: shift
    integer :: i
    real(dp) :: f,s
    i=index; f=1.0_dp; u=0.0_dp
    do while(i>0)
      f=f/real(base,dp); u=u+f*real(mod(i,base),dp); i=i/base
    end do
    s=0.0_dp; if(present(shift))s=shift
    u=modulo(u+s,1.0_dp)
    u=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u))
  end function
end module nvmix_random
