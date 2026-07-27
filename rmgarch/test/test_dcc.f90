! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
program test_dcc
   use rmgarch
   implicit none
   integer, parameter :: n = 700, m = 2
   real(dp) :: qbar(m,m), z(n,m), q(m,m,n), r(m,m,n)
   real(dp) :: ll(n), qb(m,m), nb(m,m), qf(m,m,3), rf(m,m,3)
   type(dcc_spec) :: spec, adcc_spec
   type(dcc_fit_result) :: fit
   logical :: ok
   integer :: t

   call seed_rng(12345)
   qbar = reshape([1.0_dp,0.50_dp,0.50_dp,1.0_dp],[m,m])
   spec = make_dcc_spec([0.05_dp],[0.92_dp])
   call simulate_dcc(n,spec,qbar,z,q,r,burn=400)
   call dcc_filter(z,spec,q,r,ll,qb,nb,ok)
   call assert_true(ok,'DCC filter validity')
   do t = 1, n
      call assert_close(r(1,1,t),1.0_dp,1.0e-10_dp,'DCC diagonal')
      call assert_true(abs(r(1,2,t)) < 1.0_dp,'DCC correlation bound')
   end do
   fit = fit_dcc11(z,max_iterations=700)
   call assert_true(fit%status == 0 .or. fit%status == 1,'DCC fit status')
   call assert_true(fit%spec%alpha(1) >= 0.0_dp,'DCC alpha nonnegative')
   call assert_true(fit%spec%beta(1) >= 0.0_dp,'DCC beta nonnegative')
   call assert_true(sum(fit%spec%alpha)+sum(fit%spec%beta) < 1.0_dp,'DCC stationarity')
   call dcc_forecast(fit%spec,fit%qbar,fit%nbar,fit%q(:,:,n),z(n,:),3,qf,rf)
   call assert_close(rf(1,1,1),1.0_dp,1.0e-10_dp,'forecast diagonal')

   adcc_spec = make_dcc_spec([0.03_dp],[0.91_dp],[0.02_dp],dist_student,8.0_dp)
   call simulate_dcc(n,adcc_spec,qbar,z,q,r,burn=250)
   call dcc_filter(z,adcc_spec,q,r,ll,qb,nb,ok)
   call assert_true(ok,'Student ADCC filter validity')
   call assert_true(all([(r(1,1,t) > 0.999999_dp,t=1,n)]),'Student ADCC diagonal')
   print '(a)', 'DCC tests passed.'
contains
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) error stop message
   end subroutine assert_true
   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected) > tolerance) error stop message
   end subroutine assert_close
end program test_dcc
