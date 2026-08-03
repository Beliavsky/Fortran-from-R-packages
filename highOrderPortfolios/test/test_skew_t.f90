! SPDX-License-Identifier: GPL-3.0-only
program test_skew_t
   use highorderportfolios
   use highorder_linalg, only: cholesky_upper
   implicit none
   integer,parameter :: n=3
   type(skew_t_parameters) :: p,fitp
   type(portfolio_result) :: r
   real(dp) :: lambda(4),w0(n),m0(4),f0,x(120,n)
   integer :: i,j,status
   character(len=7),parameter :: methods(6)=[character(len=7):: &
      'L-MVSK','DC','Q-MVSK','SQUAREM','RFPA','PGD']

   allocate(p%mu(n),p%gamma(n),p%scatter(n,n),p%chol_scatter(n,n))
   p%mu=[0.0015_dp,0.0008_dp,0.0011_dp]
   p%gamma=[0.0010_dp,-0.0004_dp,0.0007_dp]
   p%scatter=reshape([0.00010_dp,0.00002_dp,0.00001_dp, &
                      0.00002_dp,0.00008_dp,0.000015_dp, &
                      0.00001_dp,0.000015_dp,0.00012_dp],[n,n])
   call cholesky_upper(p%scatter,p%chol_scatter,status)
   call check(status==0,'cholesky')
   p%nu=12.0_dp
   call set_coefficients(p)
   p%status=hop_success
   w0=1.0_dp/n
   lambda=[1.0_dp,4.0_dp,10.0_dp,20.0_dp]
   m0=eval_portfolio_moments(w0,p)
   f0=-lambda(1)*m0(1)+lambda(2)*m0(2)-lambda(3)*m0(3)+lambda(4)*m0(4)
   do j=1,6
      call design_mvsk_portfolio_via_skew_t(lambda,p,r,w_init=w0, &
           method=trim(methods(j)),maxiter=300,ftol=1.0e-10_dp,wtol=1.0e-8_dp, &
           tau=1000.0_dp,initial_eta=10.0_dp)
      call check(r%status==hop_success .or. r%status==hop_not_converged,'solver status')
      call check(abs(sum(r%w)-1.0_dp)<1.0e-9_dp,'sum weights')
      call check(minval(r%w)>=-1.0e-12_dp,'nonnegative weights')
      call check(r%objective_history(size(r%objective_history))<=f0+1.0e-7_dp,'objective improvement')
   end do

   do i=1,120
      x(i,1)=0.02_dp*sin(0.13_dp*i)+0.004_dp*sin(0.031_dp*i)**2
      x(i,2)=0.015_dp*cos(0.09_dp*i)+0.002_dp*sin(0.041_dp*i)**3
      x(i,3)=0.018_dp*sin(0.07_dp*i+0.4_dp)+0.003_dp*cos(0.023_dp*i)**2
   end do
   call estimate_skew_t(x,fitp,nu_lb=9.0_dp,max_iter=30,pxem=.false.)
   call check(fitp%status==hop_success .or. fitp%status==hop_not_converged, &
      'estimate_skew_t status')
   call check(fitp%nu>8.0_dp,'finite fourth moment')
   print '(a)', 'test_skew_t: PASS'
contains
   subroutine set_coefficients(q)
      type(skew_t_parameters),intent(inout)::q
      real(dp)::nu
      nu=q%nu
      q%a11=nu/(nu-2.0_dp); q%a21=q%a11
      q%a22=2.0_dp*nu**2/((nu-2.0_dp)**2*(nu-4.0_dp))
      q%a31=16.0_dp*nu**3/((nu-2.0_dp)**3*(nu-4.0_dp)*(nu-6.0_dp))
      q%a32=6.0_dp*nu**2/((nu-2.0_dp)**2*(nu-4.0_dp))
      q%a41=(12.0_dp*nu+120.0_dp)*nu**4/ &
         ((nu-2.0_dp)**4*(nu-4.0_dp)*(nu-6.0_dp)*(nu-8.0_dp))
      q%a42=6.0_dp*(2.0_dp*nu+4.0_dp)*nu**3/ &
         ((nu-2.0_dp)**3*(nu-4.0_dp)*(nu-6.0_dp))
      q%a43=3.0_dp*nu**2/((nu-2.0_dp)*(nu-4.0_dp))
   end subroutine set_coefficients
   subroutine check(ok,label)
      logical,intent(in)::ok
      character(len=*),intent(in)::label
      if(.not.ok) then
         print '(a)', 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine check
end program test_skew_t
