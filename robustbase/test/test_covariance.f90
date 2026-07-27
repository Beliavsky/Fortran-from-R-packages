! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program test_covariance
   use robustbase
   use test_support
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer,parameter::n=80,p=2
   real(dp)::x(n,p),class_center(p),class_cov(p,p),ao(n),aoc(p),vals(p),vecs(p,p)
   type(robust_cov_result)::ogk,mcd,comed
   integer::i,info
   call seed_rng(12345)
   do i=1,n
      x(i,1)=sin(0.37_dp*real(i,dp))+0.01_dp*real(i,dp)
      x(i,2)=0.6_dp*x(i,1)+cos(0.21_dp*real(i,dp))
   end do
   x(n-7:n,1)=x(n-7:n,1)+20.0_dp
   x(n-7:n,2)=x(n-7:n,2)-15.0_dp
   class_center=[sum(x(:,1))/real(n,dp),sum(x(:,2))/real(n,dp)]
   call covariance_matrix(x,class_center,class_cov)
   call cov_ogk(x,ogk)
   call cov_comed(x,comed,n_iter=2)
   call cov_mcd(x,mcd,alpha=0.75_dp,n_starts=150,max_csteps=20)
   call assert_true(abs(ogk%center(1))<abs(class_center(1)),'OGK robust center')
   call assert_true(abs(mcd%center(1))<abs(class_center(1)),'MCD robust center')
   call assert_true(abs(comed%center(1))<abs(class_center(1)),'comedian robust center')
   call assert_true(maxval(abs(ogk%covariance-transpose(ogk%covariance)))<1.0e-10_dp,'OGK symmetry')
   call assert_true(maxval(abs(mcd%covariance-transpose(mcd%covariance)))<1.0e-10_dp,'MCD symmetry')
   call assert_true(maxval(abs(comed%covariance-transpose(comed%covariance)))<1.0e-10_dp,'comedian symmetry')
   call symmetric_eigen(ogk%covariance,vals,vecs,info)
   call assert_true(info==0 .and. all(vals>=-1.0e-10_dp),'OGK PSD')
   call symmetric_eigen(mcd%covariance,vals,vecs,info)
   call assert_true(info==0 .and. all(vals>=-1.0e-10_dp),'MCD PSD')
   call adjusted_outlyingness(x,ao,aoc,n_directions=120)
   call assert_true(sum(ao(n-7:n))/8.0_dp>sum(ao(1:n-8))/real(n-8,dp),'outliers have larger adjusted outlyingness')
   call assert_true(ieee_is_finite(cov_gk(x(:,1),x(:,2))),'GK finite')
   write(*,'(a)')'Robust covariance and outlyingness tests passed.'
end program test_covariance
