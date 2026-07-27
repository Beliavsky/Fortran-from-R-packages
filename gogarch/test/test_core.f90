! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_core
   use gogarch
   use gogarch_linalg, only : identity_matrix, is_orthogonal, is_symmetric, outer_product
   use test_helpers
   implicit none
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp) :: r(2,2), u(3,3), matched(3,3), v(6), a(3,3)
   real(dp) :: x(6,2), cov(2,2), evec(2,2), eval(2), root(2,2), invroot(2,2)
   real(dp) :: ssi(2,2,6), gamma(2,2)
   real(dp) :: y(8), residuals(8), variance(8), standardized(8), ll
   real(dp) :: expected_h2, sim_y(500), sim_h(500), mf(4), vf(4)
   type(garch11_fit) :: fit
   logical :: ok
   integer :: t

   r = rd2(pi/4.0_dp)
   call assert_true(is_orthogonal(r,1.0e-12_dp),'Rd2 must be orthogonal')
   call assert_close(r(2,1),sqrt(0.5_dp),1.0e-12_dp,'Rd2 sine placement')

   u = uprod_r([0.2_dp,0.4_dp,0.3_dp],ok)
   call assert_true(ok .and. is_orthogonal(u,1.0e-11_dp),'UprodR must form an orthogonal matrix')
   matched = umatch(identity_matrix(3),u,ok)
   call assert_true(ok .and. is_orthogonal(matched,1.0e-11_dp),'Umatch must preserve orthogonality')

   v = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   a = unvech(v,ok)
   call assert_true(ok .and. is_symmetric(a,1.0e-14_dp),'unvech must return a symmetric matrix')
   call assert_close(a(3,2),5.0_dp,0.0_dp,'unvech lower-triangle ordering')
   call assert_true(maxval(abs(vech(a)-v)) < 1.0e-14_dp,'vech and unvech must invert each other')

   x(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   x(:,2) = [2.0_dp,1.0_dp,4.0_dp,2.0_dp,5.0_dp,3.0_dp]
   call initialize_gogarch(x,cov,evec,eval,root,invroot,ok)
   call assert_true(ok,'GO-GARCH initialization must succeed')
   call assert_true(maxval(abs(matmul(root,root)-cov)) < 1.0e-10_dp,'covariance square root reconstruction')
   call assert_true(maxval(abs(matmul(invroot,matmul(cov,invroot))-identity_matrix(2))) < 1.0e-10_dp, &
      'inverse covariance square root')

   do t = 1, 6
      ssi(:,:,t) = outer_product(matmul(x(t,:),invroot),matmul(x(t,:),invroot))-identity_matrix(2)
   end do
   gamma = cora(ssi,lag=1,standardize=.true.,ok=ok)
   call assert_true(ok .and. is_symmetric(gamma,1.0e-12_dp),'cora must return a symmetric standardized matrix')

   y = [0.2_dp,-0.1_dp,0.4_dp,-0.3_dp,0.1_dp,0.0_dp,0.3_dp,-0.2_dp]
   call filter_garch11(y,0.0_dp,0.05_dp,0.10_dp,0.80_dp,residuals,variance,standardized,ll,ok)
   call assert_true(ok,'GARCH filter must accept valid parameters')
   expected_h2 = 0.05_dp+0.10_dp*y(1)**2+0.80_dp*variance(1)
   call assert_close(variance(2),expected_h2,1.0e-13_dp,'GARCH recursion')
   call assert_all_finite(standardized,'standardized GARCH residuals')

   call seed_rng(24680)
   call simulate_garch11(500,0.02_dp,0.05_dp,0.08_dp,0.87_dp,sim_y,sim_h,burnin=300)
   fit = fit_garch11(sim_y,max_iterations=220)
   call assert_true(fit%status <= 1,'GARCH fit must return a usable optimum')
   call assert_true(fit%omega > 0.0_dp .and. fit%alpha >= 0.0_dp .and. fit%beta >= 0.0_dp, &
      'GARCH fitted parameter constraints')
   call assert_true(fit%alpha+fit%beta < 1.0_dp,'GARCH covariance stationarity constraint')
   call forecast_garch11(fit,4,mf,vf)
   call assert_true(all(vf > 0.0_dp),'GARCH variance forecasts must be positive')

   write(*,'(a)') 'Core numerical tests passed.'
end program test_core
