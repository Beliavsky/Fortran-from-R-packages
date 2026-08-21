! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_random
  use mvtnorm_kinds, only : dp, pi
  implicit none
  private
  public :: seed_random, random_normal, random_normals

contains

  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=modulo(seed+104729*i+37*i*i,huge(1)-1)+1
    end do
    call random_seed(put=put)
  end subroutine seed_random

  real(dp) function random_normal() result(z)
    real(dp) :: u1,u2
    call random_number(u1); call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function random_normal

  subroutine random_normals(z)
    real(dp), intent(out) :: z(:)
    integer :: i,n
    real(dp) :: u1,u2,r
    n=size(z); i=1
    do while(i<=n)
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp)); r=sqrt(-2.0_dp*log(u1))
      z(i)=r*cos(2.0_dp*pi*u2)
      if(i+1<=n) z(i+1)=r*sin(2.0_dp*pi*u2)
      i=i+2
    end do
  end subroutine random_normals
end module mvtnorm_random
