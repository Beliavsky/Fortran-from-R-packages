! SPDX-License-Identifier: GPL-3.0-only
program test_elliptical
   use fitheavytail
   use test_support, only: check, make_data
   implicit none
   integer :: i
   real(dp) :: x(240,3)
   type(heavy_tail_fit) :: tyler, cauchy

   call make_data(x)
   call fit_tyler(x,tyler,max_iter=300,ptol=2.0e-4_dp)
   call check(tyler%status==ht_success.or. &
      tyler%status==ht_no_convergence,'Tyler status')
   call check(allocated(tyler%covariance),'Tyler covariance')
   call check(all([(tyler%covariance(i,i)>0.0_dp,i=1,3)]), &
      'Tyler positive diagonal')
   call check(abs(sum(tyler%mu)-0.4_dp)<1.0_dp,'Tyler location')

   call fit_cauchy(x,cauchy,max_iter=300,ptol=2.0e-4_dp)
   call check(cauchy%status==ht_success.or. &
      cauchy%status==ht_no_convergence,'Cauchy status')
   call check(allocated(cauchy%scatter),'Cauchy scatter')
   call check(all([(cauchy%scatter(i,i)>0.0_dp,i=1,3)]), &
      'Cauchy positive diagonal')
   call check(maxval(abs(tyler%covariance-cauchy%covariance))<3.0_dp, &
      'robust covariance agreement')
   write(*,'(a)') 'test_elliptical: PASS'
end program test_elliptical
