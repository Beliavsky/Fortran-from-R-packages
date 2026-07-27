! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program test_next_batch
   use robustbase
   use test_support
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer,parameter::n=80
   real(dp)::a(6,3),x(n,2),z6(6,3),y(n),u,prob,h,fd1,fd2,eps,vals(2),vecs(2,2)
   type(pca_result)::pca
   type(detmcd_result)::dmcd1,dmcd2
   type(fast_lts_result)::flts_exact,flts_det
   type(by_logistic_result)::byfit
   integer::i,info

   a(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   a(:,2)=[2.0_dp,1.0_dp,0.0_dp,-1.0_dp,-2.0_dp,-3.0_dp]
   a(:,3)=a(:,1)+2.0_dp*a(:,2)
   call assert_true(rank_mm(a)==2,'rankMM dependent-column rank')
   call classical_pca(a,pca,center_data=.true.,scale_data=.false.,return_scores=.true.)
   call assert_true(pca%rank==1,'classPC centered rank')
   z6=matmul(pca%scores,transpose(pca%loadings))
   do i=1,3
      z6(:,i)=z6(:,i)+pca%center(i)
   end do
   call assert_true(maxval(abs(z6-a))<1.0e-10_dp,'classPC reconstruction')
   if(pca%rank>1) call assert_true(all(pca%eigenvalues(1:pca%rank-1)>=pca%eigenvalues(2:pca%rank)),'classPC descending eigenvalues')

   do i=1,n
      x(i,1)=sin(0.17_dp*real(i,dp))+0.01_dp*real(i,dp)
      x(i,2)=0.7_dp*x(i,1)+cos(0.11_dp*real(i,dp))
   end do
   x(n-9:n,1)=x(n-9:n,1)+18.0_dp
   x(n-9:n,2)=x(n-9:n,2)-14.0_dp
   call cov_detmcd(x,dmcd1,alpha=0.75_dp)
   call cov_detmcd(x,dmcd2,alpha=0.75_dp)
   call assert_true(dmcd1%best_start>=1 .and. dmcd1%best_start<=6,'detMCD selected six-pack start')
   call assert_true(all(shape(dmcd1%initial_orderings)==[n,6]),'detMCD saved six orderings')
   call assert_true(all(dmcd1%best_subset==dmcd2%best_subset),'detMCD deterministic subset')
   call assert_true(maxval(abs(dmcd1%estimate%raw_center-dmcd2%estimate%raw_center))<1.0e-13_dp,'detMCD deterministic center')
   call assert_true(abs(dmcd1%estimate%center(1))<abs(sum(x(:,1))/real(n,dp)),'detMCD robust center')
   call symmetric_eigen(dmcd1%estimate%covariance,vals,vecs,info)
   call assert_true(info==0 .and. all(vals>=-1.0e-10_dp),'detMCD covariance PSD')

   do i=1,12
      x(i,1)=1.0_dp
      x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/11.0_dp
      y(i)=0.5_dp+2.0_dp*x(i,2)+0.01_dp*sin(real(i,dp))
   end do
   y(10:12)=y(10:12)+12.0_dp
   call fast_lts_regression(x(1:12,:),y(1:12),flts_exact,alpha=0.5_dp,sampling='exact',max_csteps=30)
   call fast_lts_regression(x(1:12,:),y(1:12),flts_det,alpha=0.5_dp,sampling='deterministic',n_starts=40,max_csteps=30)
   call assert_true(flts_exact%exhaustive .and. flts_exact%trials==66,'FAST-LTS exact enumeration')
   call assert_true(abs(flts_exact%coefficients(2)-2.0_dp)<0.08_dp,'FAST-LTS exact slope')
   call assert_true(abs(flts_det%coefficients(2)-2.0_dp)<0.12_dp,'FAST-LTS deterministic slope')
   call assert_true(flts_exact%objective<=flts_det%objective+1.0e-10_dp,'FAST-LTS exact objective dominates')

   eps=1.0e-5_dp
   h=0.37_dp
   fd1=(by_phi(h+eps,1.0_dp,0.5_dp)-by_phi(h-eps,1.0_dp,0.5_dp))/(2.0_dp*eps)
   fd2=(by_phi_derivative(h+eps,1.0_dp,0.5_dp)-by_phi_derivative(h-eps,1.0_dp,0.5_dp))/(2.0_dp*eps)
   call assert_close(by_phi_derivative(h,1.0_dp,0.5_dp),fd1,2.0e-7_dp,'BY first derivative')
   call assert_close(by_phi_second_derivative(h,1.0_dp,0.5_dp),fd2,3.0e-6_dp,'BY second derivative')

   call seed_rng(97531)
   do i=1,n
      x(i,1)=1.0_dp
      x(i,2)=-2.5_dp+5.0_dp*real(i-1,dp)/real(n-1,dp)
      prob=1.0_dp/(1.0_dp+exp(-(-0.35_dp+1.25_dp*x(i,2))))
      call random_number(u)
      y(i)=merge(1.0_dp,0.0_dp,u<prob)
   end do
   y(n-7:n)=0.0_dp
   call by_logistic_fit(x,y,byfit,const=0.5_dp,max_iter=400,tol=1.0e-7_dp)
   call assert_true(byfit%converged,'BY logistic convergence')
   call assert_true(byfit%coefficients(2)>0.35_dp,'BY logistic positive slope')
   call assert_true(ieee_is_finite(byfit%objective) .and. all(ieee_is_finite(byfit%standard_errors)),'BY finite inference')
   call assert_true(maxval(abs(byfit%covariance-transpose(byfit%covariance)))<1.0e-9_dp,'BY covariance symmetry')
   write(*,'(a)')'PCA/rank, deterministic MCD, advanced LTS, and BY logistic tests passed.'
end program test_next_batch
