! SPDX-License-Identifier: GPL-3.0-only
program test_estimation
   use smoots
   implicit none
   integer,parameter::n=220
   real(dp)::y(n),x
   integer::i
   type(smooth_result)::a,b,o,na,om,nm,oam,nam,d1,d2
   type(confidence_result)::ci
   do i=1,n
      x=real(i,dp)/real(n,dp)
      y(i)=1.0_dp+2.0_dp*x+0.4_dp*sin(6.0_dp*x)+0.08_dp*sin(47.0_dp*x)
   end do
   call msmooth(y,a,p=1,algorithm='A')
   call msmooth(y,b,p=3,algorithm='B')
   call msmooth(y,o,p=1,algorithm='O')
   call msmooth(y,na,p=1,algorithm='NA')
   call msmooth(y,om,p=1,algorithm='OM')
   call msmooth(y,nm,p=1,algorithm='NM')
   call msmooth(y,oam,p=1,algorithm='OAM')
   call msmooth(y,nam,p=1,algorithm='NAM')
   call check(a,'A');call check(b,'B');call check(o,'O');call check(na,'NA')
   call check(om,'OM');call check(nm,'NM');call check(oam,'OAM');call check(nam,'NAM')
   call dsmooth(y,d1,d=1);call dsmooth(y,d2,d=2)
   call check(d1,'d1');call check(d2,'d2')
   call conf_bounds(a,0.95_dp,ci,1)
   if(ci%status/=sm_ok.or.any(ci%lower>ci%upper))error stop 'confidence bounds'
   if(maxval(abs(a%estimate-(y-a%residuals)))>1.0e-12_dp)error stop 'residual identity'
   print '(a)','test_estimation: PASS'
contains
   subroutine check(r,name)
      type(smooth_result),intent(in)::r
      character(len=*),intent(in)::name
      if(r%status/=sm_ok.and.r%status/=sm_iteration_limit)then;print *,name,r%status;error stop 1;end if
      if(r%b0<=0.0_dp.or.r%b0>=0.5_dp)error stop 'bandwidth range'
      if(.not.allocated(r%estimate))error stop 'missing estimate'
   end subroutine check
end program test_estimation
