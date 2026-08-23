! SPDX-License-Identifier: GPL-3.0-or-later
program test_core
   use pracma
   implicit none
   real(dp),allocatable :: x(:),ai(:,:),rr(:)
   complex(dp),allocatable :: z(:)
   integer :: st

   x=linspace(0.0_dp,1.0_dp,5)
   call check(maxval(abs(x-[0.0_dp,0.25_dp,0.5_dp,0.75_dp,1.0_dp]))<1e-14_dp,'linspace')
   call check(abs(determinant(reshape([4.0_dp,2.0_dp,7.0_dp,6.0_dp],[2,2]))-10.0_dp)<1e-12_dp,'determinant')
   ai=inv(reshape([4.0_dp,2.0_dp,7.0_dp,6.0_dp],[2,2]),st)
   call check(st==pracma_ok,'inverse status')
   call check(maxval(abs(matmul(reshape([4.0_dp,2.0_dp,7.0_dp,6.0_dp],[2,2]),ai)-eye(2)))<1e-11_dp,'inverse')
   call check(abs(polyval([1.0_dp,-3.0_dp,2.0_dp],2.0_dp))<1e-14_dp,'polyval')
   z=roots([1.0_dp,-3.0_dp,2.0_dp],status=st)
   call check(st==pracma_ok,'roots status')
   rr=sort_real(real(z,dp))
   call check(maxval(abs(rr-[1.0_dp,2.0_dp]))<1e-8_dp,'roots')
   call check(abs(gammainc(1.0_dp,1.0_dp)-(1.0_dp-exp(-1.0_dp)))<1e-12_dp,'gammainc')
   call check(abs(zeta(2.0_dp)-pi_dp*pi_dp/6.0_dp)<1e-8_dp,'zeta')
   call check(abs(agmean_value(1.0_dp,2.0_dp)-1.4567910310469068_dp)<1e-12_dp,'agmean')
   print '(a)','test_core: PASS'
contains
   real(dp) function agmean_value(a,b) result(v)
      real(dp),intent(in)::a,b
      call agmean(a,b,v)
   end function agmean_value
   subroutine check(ok,name)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      if(.not.ok)then; write(*,'(a,1x,a)')'FAIL:',trim(name); error stop 1; end if
   end subroutine check
end program test_core
