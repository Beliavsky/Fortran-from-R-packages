! SPDX-License-Identifier: GPL-3.0-only
program moebius_example
  use mixedindtests
  implicit none
  real(dp) :: x(20,3)
  type(moebius_result) :: out
  integer :: i,k

  do i=1,20
    x(i,1)=real(mod(i,4),dp)
    x(i,2)=sin(real(i,dp))
    x(i,3)=real(mod(3*i+1,7),dp)
  end do
  out=EstDepMoebius(x,3)
  do k=1,size(out%spearman)
    print '(i3,1x,i2,3(1x,f10.5))',k,out%cardinality(k), &
      out%spearman(k),out%vdw(k),out%savage(k)
  end do
end program moebius_example
