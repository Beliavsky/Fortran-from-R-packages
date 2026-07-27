! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_distributions
  use fmultivar, only : dp, i8, pi, dnorm2d, pnorm2d, dt2d, pt2d, dcauchy2d, pcauchy2d, &
    elliptical2d_density, mvnorm_pdf, mvt_pdf, mvnorm_rng, mvt_rng, mvnorm_rect_prob, &
    mvnorm_equicoordinate_quantile, sample_mean_cov, mvsnorm_pdf, mvst_pdf, &
    mvsnorm_rng, mvst_rng, mvsnorm_rect_prob
  implicit none
  real(dp) :: rho,p,se,q,f1,f2,p_exact
  real(dp), allocatable :: x(:,:),mean(:),cov(:,:),lower(:),upper(:)
  real(dp) :: mu2(2),om2(2,2),alpha2(2),point2(2)
  logical :: ok
  integer :: failures
  failures=0;rho=0.6_dp
  call check_close(pnorm2d(0.0_dp,0.0_dp,rho),0.25_dp+asin(rho)/(2.0_dp*pi),2.0e-14_dp,'bvn origin',failures)
  call check_close(pnorm2d(0.7_dp,-0.2_dp,0.0_dp), &
    0.5_dp*erfc(-0.7_dp/sqrt(2.0_dp))*0.5_dp*erfc(0.2_dp/sqrt(2.0_dp)),1.0e-12_dp,'bvn independent',failures)
  call check_close(pt2d(0.7_dp,-0.2_dp,0.0_dp,5.0_dp), &
    student_cdf_local(0.7_dp,5.0_dp)*student_cdf_local(-0.2_dp,5.0_dp),2.0e-11_dp,'bvt independent',failures)
  call check_close(dcauchy2d(0.4_dp,-0.3_dp,rho),dt2d(0.4_dp,-0.3_dp,rho,1.0_dp),1.0e-14_dp,'cauchy density',failures)
  call check_close(pcauchy2d(0.4_dp,-0.3_dp,rho),pt2d(0.4_dp,-0.3_dp,rho,1.0_dp),1.0e-12_dp,'cauchy cdf',failures)
  call check_close(elliptical2d_density(0.4_dp,-0.3_dp,rho,'norm'), &
    dnorm2d(0.4_dp,-0.3_dp,rho),1.0e-14_dp,'elliptical normal',failures)
  call check_close(elliptical2d_density(0.4_dp,-0.3_dp,rho,'t',5.0_dp), &
    dt2d(0.4_dp,-0.3_dp,rho,5.0_dp),1.0e-14_dp,'elliptical t',failures)
  call check_positive(elliptical2d_density(0.4_dp,-0.3_dp,rho,'logistic'),'elliptical logistic',failures)
  call check_positive(elliptical2d_density(0.4_dp,-0.3_dp,rho,'laplace'),'elliptical laplace',failures)
  call check_positive(elliptical2d_density(0.4_dp,-0.3_dp,rho,'kotz',1.2_dp),'elliptical kotz',failures)
  call check_positive(elliptical2d_density(0.4_dp,-0.3_dp,rho,'epower',1.1_dp,0.8_dp),'elliptical epower',failures)

  mu2=[0.2_dp,-0.1_dp];om2=reshape([1.0_dp,0.4_dp,0.4_dp,2.0_dp],[2,2]);point2=[0.3_dp,-0.4_dp]
  f1=mvnorm_pdf(point2,mu2,om2,ok);call check_true(ok.and.f1>0.0_dp,'mvnorm density',failures)
  f2=mvt_pdf(point2,mu2,om2,7.0_dp,ok);call check_true(ok.and.f2>0.0_dp,'mvt density',failures)
  call mvnorm_rng(40000,mu2,om2,x,12345_i8,ok);call check_true(ok,'mvnorm rng status',failures)
  call sample_mean_cov(x,mean,cov)
  call check_close(maxval(abs(mean-mu2)),0.0_dp,0.03_dp,'mvnorm rng mean',failures)
  call check_close(maxval(abs(cov-om2)),0.0_dp,0.06_dp,'mvnorm rng covariance',failures)
  call mvt_rng(50000,mu2,om2,8.0_dp,x,23456_i8,ok);call sample_mean_cov(x,mean,cov)
  call check_close(maxval(abs(mean-mu2)),0.0_dp,0.04_dp,'mvt rng mean',failures)
  call check_close(maxval(abs(cov-om2*8.0_dp/6.0_dp)),0.0_dp,0.12_dp,'mvt rng covariance',failures)

  allocate(lower(2),upper(2));lower=-huge(1.0_dp);upper=[0.2_dp,0.3_dp]
  call mvnorm_rect_prob(lower,upper,[0.0_dp,0.0_dp], &
    reshape([1.0_dp,rho,rho,1.0_dp],[2,2]),p,se,180000,999_i8,ok)
  p_exact=pnorm2d(0.2_dp,0.3_dp,rho)
  call check_true(ok.and.abs(p-p_exact)<max(0.012_dp,5.0_dp*se), &
    'mvnorm rectangle',failures)
  q=mvnorm_equicoordinate_quantile(0.25_dp,[0.0_dp,0.0_dp],reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]),70000)
  call check_close(q,0.0_dp,0.045_dp,'equicoordinate quantile',failures)

  alpha2=[0.0_dp,0.0_dp]
  call check_close(mvsnorm_pdf(point2,mu2,om2,alpha2),mvnorm_pdf(point2,mu2,om2),1.0e-13_dp,'skew normal alpha zero',failures)
  call check_close(mvst_pdf(point2,mu2,om2,alpha2,7.0_dp),mvt_pdf(point2,mu2,om2,7.0_dp),1.0e-13_dp,'skew t alpha zero',failures)
  alpha2=[4.0_dp,-1.0_dp]
  call mvsnorm_rng(35000,mu2,om2,alpha2,x,34567_i8,ok);call sample_mean_cov(x,mean,cov)
  call check_true(ok.and.mean(1)>mu2(1)+0.35_dp,'skew normal rng direction',failures)
  call mvst_rng(35000,mu2,om2,alpha2,8.0_dp,x,45678_i8,ok);call sample_mean_cov(x,mean,cov)
  call check_true(ok.and.mean(1)>mu2(1)+0.35_dp,'skew t rng direction',failures)
  lower=-huge(1.0_dp);upper=[0.5_dp,0.5_dp]
  call mvsnorm_rect_prob(lower,upper,mu2,om2,alpha2,p,se,60000,888_i8,ok)
  call check_true(ok.and.p>=0.0_dp.and.p<=1.0_dp.and.se<0.01_dp,'skew normal probability',failures)

  if(failures>0)then
    write(*,'(a,i0)')'Distribution test failures: ',failures;error stop 1
  end if
  write(*,'(a)')'Distribution tests passed.'
contains
  function student_cdf_local(z,nu) result(v)
    use fmultivar, only : student_t_cdf
    real(dp),intent(in)::z,nu
    real(dp)::v
    v=student_t_cdf(z,nu)
  end function student_cdf_local
  subroutine check_close(actual,expected,tol,name,nfail)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.(abs(actual-expected)<=tol))then
      write(*,'(a,2es18.8)')trim(name)//' failed: ',actual,expected;nfail=nfail+1
    end if
  end subroutine check_close
  subroutine check_positive(value,name,nfail)
    real(dp),intent(in)::value
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.(value>0.0_dp))then;write(*,'(a)')trim(name)//' failed';nfail=nfail+1;end if
  end subroutine check_positive
  subroutine check_true(cond,name,nfail)
    logical,intent(in)::cond
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.cond)then;write(*,'(a)')trim(name)//' failed';nfail=nfail+1;end if
  end subroutine check_true
end program test_distributions
