! SPDX-License-Identifier: GPL-3.0-only
program test_tilting_validation
   use highorderportfolios
   implicit none
   integer,parameter :: t=140,n=4
   real(dp) :: x(t,n),w0(n),m0(4),d(4),kappa,diff(n),tracking
   type(sample_moments) :: s
   type(portfolio_result) :: r
   integer :: i,j
   character(len=7),parameter :: methods(2)=[character(len=7)::'Q-MVSKT','L-MVSKT']

   do i=1,t
      x(i,1)=0.0015_dp+0.012_dp*sin(0.09_dp*i)+0.004_dp*sin(0.021_dp*i)**2
      x(i,2)=0.0008_dp+0.008_dp*cos(0.13_dp*i)
      x(i,3)=0.0003_dp+0.014_dp*sin(0.06_dp*i+0.9_dp)
      x(i,4)=0.0011_dp+0.006_dp*cos(0.04_dp*i)+0.002_dp*sin(0.15_dp*i)**3
   end do
   call estimate_sample_moments(x,s,adjust_magnitude=.true.)
   w0=1.0_dp/n
   m0=eval_portfolio_moments(w0,s)
   d=abs(m0)
   kappa=0.3_dp*sqrt(dot_product(w0,matmul(s%covariance,w0)))
   do j=1,2
      call design_mvsktilting_portfolio_via_sample_moments(d,s,r,w_init=w0,w0=w0, &
           w0_moments=m0,kappa=kappa,method=trim(methods(j)),maxiter=150)
      call check(r%status==hop_success .or. r%status==hop_not_converged,'tilting status')
      call check(abs(sum(r%w)-1.0_dp)<1.0e-10_dp,'sum weights')
      call check(minval(r%w)>=-1.0e-12_dp,'nonnegative weights')
      diff=r%w-w0
      tracking=dot_product(diff,matmul(s%covariance,diff))
      call check(tracking<=kappa*kappa*(1.0_dp+1.0e-9_dp)+1.0e-14_dp,'tracking constraint')
      call check(abs(r%delta-minval(r%improvement))<1.0e-12_dp,'delta definition')
   end do

   call design_mvsk_portfolio_via_sample_moments([1.0_dp,1.0_dp,1.0_dp,1.0_dp],s,r,leverage=1.2_dp)
   call check(r%status==hop_invalid_argument,'invalid leverage')
   print '(a)', 'test_tilting_validation: PASS'
contains
   subroutine check(ok,label)
      logical,intent(in)::ok
      character(len=*),intent(in)::label
      if(.not.ok) then
         print '(a)', 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine check
end program test_tilting_validation
