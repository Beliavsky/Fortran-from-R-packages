! SPDX-License-Identifier: GPL-3.0-only
program test_sample_design
   use highorderportfolios
   implicit none
   integer,parameter :: t=160,n=4
   real(dp) :: x(t,n),lmd(4),w0(n),m0(4),f0
   type(sample_moments) :: s
   type(portfolio_result) :: r
   integer :: i,j
   character(len=6),parameter :: methods(3)=[character(len=6)::'Q-MVSK','MM','DC']

   do i=1,t
      x(i,1)=0.0012_dp+0.010_dp*sin(0.13_dp*i)+0.002_dp*sin(0.037_dp*i)**2
      x(i,2)=0.0007_dp+0.007_dp*cos(0.11_dp*i)
      x(i,3)=0.0004_dp+0.012_dp*sin(0.07_dp*i+0.8_dp)
      x(i,4)=0.0009_dp+0.006_dp*cos(0.05_dp*i)+0.001_dp*sin(0.17_dp*i)**3
   end do
   call estimate_sample_moments(x,s)
   call check(s%status==hop_success,'moment estimation')
   w0=1.0_dp/n
   lmd=[1.0_dp,5.0_dp,18.333333333333333_dp,55.0_dp]
   m0=eval_portfolio_moments(w0,s)
   f0=-lmd(1)*m0(1)+lmd(2)*m0(2)-lmd(3)*m0(3)+lmd(4)*m0(4)
   do j=1,3
      call design_mvsk_portfolio_via_sample_moments(lmd,s,r,w_init=w0, &
           method=trim(methods(j)),maxiter=120,ftol=1.0e-10_dp,wtol=1.0e-9_dp)
      call check(r%status==hop_success .or. r%status==hop_not_converged,'solver status')
      call check(abs(sum(r%w)-1.0_dp)<1.0e-10_dp,'sum weights')
      call check(minval(r%w)>=-1.0e-12_dp,'nonnegative weights')
      call check(r%objective_history(size(r%objective_history))<=f0+1.0e-8_dp,'objective improvement')
   end do
   print '(a)', 'test_sample_design: PASS'
contains
   subroutine check(ok,label)
      logical,intent(in)::ok
      character(len=*),intent(in)::label
      if(.not.ok) then
         print '(a)', 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine check
end program test_sample_design
