! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_random
  use tvmvp_kinds, only : dp
  implicit none
  private
  public :: seed_random, random_normal, fill_random_normal
contains
  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=modulo(seed+104729*i,2147483646)+1
    end do
    call random_seed(put=put)
  end subroutine seed_random

  real(dp) function random_normal()
    real(dp) :: u1,u2
    real(dp), parameter :: twopi=6.2831853071795864769252867665590058_dp
    call random_number(u1); call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    random_normal=sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
  end function random_normal

  subroutine fill_random_normal(x)
    real(dp), intent(out) :: x(:,:)
    integer :: i,j
    do j=1,size(x,2)
      do i=1,size(x,1)
        x(i,j)=random_normal()
      end do
    end do
  end subroutine fill_random_normal
end module tvmvp_random
