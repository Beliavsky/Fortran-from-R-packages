! SPDX-License-Identifier: GPL-3.0-or-later
module test_ode_functions
   use pracma_kinds, only : dp
   implicit none
contains
   subroutine decay(t,y,dydt)
      real(dp),intent(in)::t,y(:); real(dp),intent(out)::dydt(:)
      dydt=-y+0.0_dp*t
   end subroutine decay
end module test_ode_functions

program test_interpolation_ode_signal
   use pracma
   use test_ode_functions
   implicit none
   real(dp),allocatable::yi(:),cv(:),sm(:),p(:)
   complex(dp),allocatable::f(:),back(:)
   type(ode_result)::odes
   type(peak_result)::peaks
   integer :: i
   yi=interp1([0.0_dp,1.0_dp,2.0_dp],[0.0_dp,1.0_dp,4.0_dp],[0.5_dp,1.5_dp],'linear')
   call check(maxval(abs(yi-[0.5_dp,2.5_dp]))<1e-14_dp,'interp1 linear')
   yi=interp1_pchip([0.0_dp,1.0_dp,2.0_dp],[0.0_dp,1.0_dp,4.0_dp],[0.0_dp,1.0_dp,2.0_dp])
   call check(maxval(abs(yi-[0.0_dp,1.0_dp,4.0_dp]))<1e-13_dp,'pchip knots')
   odes=ode45(decay,[0.0_dp,1.0_dp],[1.0_dp])
   call check(odes%converged.and.abs(odes%y(1,size(odes%t))-exp(-1.0_dp))<2e-6_dp,'ode45')
   cv=conv([1.0_dp,2.0_dp],[1.0_dp,1.0_dp])
   call check(maxval(abs(cv-[1.0_dp,3.0_dp,2.0_dp]))<1e-14_dp,'conv')
   f=fft(cmplx([1.0_dp,2.0_dp,3.0_dp,4.0_dp],0.0_dp,dp)); back=ifft(f)
   call check(maxval(abs(real(back,dp)-[1.0_dp,2.0_dp,3.0_dp,4.0_dp]))<1e-12_dp,'fft inverse')
   sm=movavg([1.0_dp,2.0_dp,100.0_dp,4.0_dp,5.0_dp],3)
   call check(size(sm)==5,'movavg')
   peaks=findpeaks([0.0_dp,1.0_dp,0.0_dp,2.0_dp,0.0_dp])
   call check(size(peaks%indices)==2.and.all(peaks%indices==[2,4]),'findpeaks')
   p=periodogram(sin(2*pi_dp*[(real(i,dp),i=0,15)]/16.0_dp))
   call check(maxloc(p,dim=1)==2,'periodogram')
   print '(a)','test_interpolation_ode_signal: PASS'
contains
   subroutine check(ok,name)
      logical,intent(in)::ok; character(len=*),intent(in)::name
      if(.not.ok)then; write(*,'(a,1x,a)')'FAIL:',trim(name); error stop 1; end if
   end subroutine check
end program test_interpolation_ode_signal
