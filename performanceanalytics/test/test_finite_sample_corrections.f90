! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program test_finite_sample_corrections
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use performanceanalytics_mod, only: dp, mean_value, exact_shrinkage_result, &
    exact_vm2_terms, exact_vm3_terms, exact_vm3_kstat_terms, exact_vm4_terms, &
    exact_m2_shrinkage, exact_m3_shrinkage, exact_m4_shrinkage, structured_coskewness
  use test_support_mod, only: assert_close, assert_true, assert_vector_close
  implicit none
  integer,parameter::n=12,p=4
  real(dp)::r(n,p),x(n,p),f(n,2)
  real(dp)::m11(p,p),m21(p,p),m22(p,p),m31(p,p),m32(p,p),m41(p,p),m42(p,p),m33(p,p)
  real(dp)::terms2(3),terms3(3),terms3u(3),terms4(3),ref2(3),ref3(3),ref4(3),sm3(p,p*p)
  type(exact_shrinkage_result)::fit
  integer::t,i,j

  do t=1,n
    r(t,1)=sin(0.37_dp*real(t,dp))+0.02_dp*real(t,dp)
    r(t,2)=cos(0.23_dp*real(t,dp))-0.03_dp*real(t,dp)+0.2_dp*r(t,1)
    r(t,3)=0.4_dp*real(mod(t,5)-2,dp)+0.1_dp*sin(0.71_dp*real(t,dp))
    r(t,4)=0.3_dp*r(t,1)-0.2_dp*r(t,2)+sin(0.11_dp*real(t*t,dp))
    f(t,1)=0.5_dp*r(t,1)-0.3_dp*r(t,2)+0.2_dp*r(t,3)
    f(t,2)=0.4_dp*r(t,4)+0.1_dp*cos(0.13_dp*real(t,dp))
  end do
  do j=1,p;x(:,j)=r(:,j)-mean_value(r(:,j));end do
  do i=1,p;do j=1,p
    m11(i,j)=sum(x(:,i)*x(:,j))/real(n,dp)
    m21(i,j)=sum(x(:,i)**2*x(:,j))/real(n,dp)
    m22(i,j)=sum(x(:,i)**2*x(:,j)**2)/real(n,dp)
    m31(i,j)=sum(x(:,i)**3*x(:,j))/real(n,dp)
    m32(i,j)=sum(x(:,i)**3*x(:,j)**2)/real(n,dp)
    m41(i,j)=sum(x(:,i)**4*x(:,j))/real(n,dp)
    m42(i,j)=sum(x(:,i)**4*x(:,j)**2)/real(n,dp)
    m33(i,j)=sum(x(:,i)**3*x(:,j)**3)/real(n,dp)
  end do;end do

  call exact_vm2_terms(m11,m22,n,terms2)
  call direct_vm2(x,m11,ref2)
  call assert_vector_close(terms2,ref2,2.0e-12_dp,'exact VM2 influence identities')

  call exact_vm3_terms(x,m11,m21,m22,m31,m42,m33,terms3)
  call direct_vm3(x,m11,ref3)
  call assert_vector_close(terms3,ref3,2.0e-11_dp,'exact VM3 influence identities')

  call exact_vm3_kstat_terms(x,real(n,dp)*m11,real(n,dp)*m21,real(n,dp)*m22, &
    real(n,dp)*m31,real(n,dp)*m42,real(n,dp)*m33,terms3u)
  call assert_vector_close(terms3u,[0.65728955243026177_dp,0.043782526665478896_dp, &
    0.14117214288281874_dp],2.0e-12_dp,'unbiased VM3 fixed reference')
  call structured_coskewness(r,'indep',sm3,unbiased_marg=.true.)
  call assert_close(sm3(1,1),real(n*n,dp)/real((n-1)*(n-2),dp)*sum(x(:,1)**3)/real(n,dp), &
    2.0e-12_dp,'unbiased structured marginal coskewness')

  call exact_vm4_terms(x,m11,m21,m22,m31,m32,m41,m42,terms4)
  call direct_vm4(x,m11,ref4)
  call assert_vector_close(terms4,ref4,3.0e-10_dp,'exact VM4 influence identities')

  call exact_m2_shrinkage(r,[1,2,3,4],fit,f)
  call check_fit(fit,2,5,'M2 multi-target shrinkage')
  call assert_vector_close(fit%b,[0.15561624598817927_dp,0.19508173594644207_dp, &
    0.084095526913055801_dp,0.13058707618305826_dp,0.050445935187985042_dp], &
    3.0e-11_dp,'M2 exact target-covariance reference')
  call exact_m2_shrinkage(r,[1],fit)
  call assert_close(fit%lambda(1),max(0.0_dp,min(1.0_dp,fit%b(1)/fit%a(1,1))),1.0e-12_dp,'M2 scalar QP')

  call exact_m3_shrinkage(r,[1,2,3,4,5,6],fit,f)
  call check_fit(fit,3,7,'M3 multi-target shrinkage')
  call assert_vector_close(fit%b,[0.36119072603242464_dp,0.42995357140345181_dp, &
    0.36083102579074133_dp,0.35936395485612815_dp,0.13460187600334894_dp, &
    0.46900746585785286_dp,0.29700671042630306_dp],4.0e-10_dp, &
    'M3 exact target-covariance reference')
  call exact_m3_shrinkage(r,[1,2,6],fit,unbiased_mse=.true.)
  call check_fit(fit,3,3,'M3 unbiased shrinkage')
  call assert_vector_close(fit%b,[0.51611740954744323_dp,0.61350702576478311_dp, &
    0.65728955243026199_dp],3.0e-11_dp,'M3 unbiased target reference')
  call assert_true(fit%unbiased_mse,'M3 unbiased flag')
  call exact_m3_shrinkage(r,[6],fit)
  call assert_close(fit%lambda(1),max(0.0_dp,min(1.0_dp,fit%b(1)/fit%a(1,1))),1.0e-12_dp,'M3 scalar QP')

  call exact_m4_shrinkage(r,[1,2,3,4],fit,f)
  call check_fit(fit,4,5,'M4 multi-target shrinkage')
  call assert_vector_close(fit%b,[1.3403724588604473_dp,1.4926017432519485_dp, &
    0.72851239871241280_dp,1.1297299026601697_dp,0.23812117485351991_dp], &
    8.0e-9_dp,'M4 exact target-covariance reference')
  call exact_m4_shrinkage(r,[1],fit)
  call assert_close(fit%lambda(1),max(0.0_dp,min(1.0_dp,fit%b(1)/fit%a(1,1))),1.0e-12_dp,'M4 scalar QP')

  print '(a)','Exact finite-sample M2/M3/M4 shrinkage-correction tests passed.'
contains
  subroutine check_fit(z,order,n_targets,label)
    type(exact_shrinkage_result),intent(in)::z
    integer,intent(in)::order,n_targets
    character(len=*),intent(in)::label
    call assert_true(z%order==order,trim(label)//' order')
    call assert_true(z%n_targets==n_targets,trim(label)//' target count')
    call assert_true(z%converged,trim(label)//' convergence')
    call assert_true(all(z%lambda>=-1.0e-12_dp),trim(label)//' nonnegative weights')
    call assert_true(sum(z%lambda)<=1.0_dp+1.0e-10_dp,trim(label)//' simplex weight')
    call assert_true(all(ieee_is_finite(z%lambda)) .and. all(ieee_is_finite(z%a)) .and. &
      all(ieee_is_finite(z%b)) .and. all(ieee_is_finite(z%estimate)),trim(label)//' finite outputs')
    call assert_close(maxval(abs(z%a-transpose(z%a))),0.0_dp,2.0e-12_dp,trim(label)//' symmetric A')
  end subroutine check_fit

  pure integer function mult3(i,j,k) result(v)
    integer,intent(in)::i,j,k
    if(i==k)then;v=1;else if(i==j .or. j==k)then;v=3;else;v=6;end if
  end function mult3

  pure integer function mult4(i,j,k,l) result(v)
    integer,intent(in)::i,j,k,l
    if(i==l)then
      v=1
    else if((i==k .and. k<l) .or. (i<j .and. j==l))then
      v=4
    else if(i==j .and. k==l .and. j<k)then
      v=6
    else if(i==j .or. j==k .or. k==l)then
      v=12
    else
      v=24
    end if
  end function mult4

  subroutine direct_vm2(xc,cov,ans)
    real(dp),intent(in)::xc(:,:),cov(:,:)
    real(dp),intent(out)::ans(3)
    real(dp)::ifs(size(xc,1)),ift(size(xc,1)),avg(size(xc,1))
    integer::ii,jj,tt,pp,nn
    nn=size(xc,1);pp=size(xc,2);ans=0.0_dp;avg=0.0_dp
    do tt=1,nn
      do ii=1,pp;avg(tt)=avg(tt)+(xc(tt,ii)**2-cov(ii,ii))/real(pp,dp);end do
    end do
    do ii=1,pp;do jj=1,pp
      ifs=xc(:,ii)*xc(:,jj)-cov(ii,jj)
      ans(1)=ans(1)+sum(ifs*ifs)/real(nn*nn,dp)
      if(ii==jj)then
        ift=xc(:,ii)**2-cov(ii,ii)
        ans(3)=ans(3)+sum(ifs*ift)/real(nn*nn,dp)
        ans(2)=ans(2)+sum(ifs*avg)/real(nn*nn,dp)
      end if
    end do;end do
  end subroutine direct_vm2

  subroutine direct_vm3(xc,cov,ans)
    real(dp),intent(in)::xc(:,:),cov(:,:)
    real(dp),intent(out)::ans(3)
    real(dp),allocatable::m3(:,:,:),if_diag(:,:),avg(:),ifs(:),ift(:)
    integer::ii,jj,kk,pp,nn,m
    nn=size(xc,1);pp=size(xc,2);allocate(m3(pp,pp,pp),if_diag(nn,pp),avg(nn),ifs(nn),ift(nn))
    do ii=1,pp;do jj=1,pp;do kk=1,pp
      m3(ii,jj,kk)=sum(xc(:,ii)*xc(:,jj)*xc(:,kk))/real(nn,dp)
    end do;end do;end do
    do ii=1,pp
      if_diag(:,ii)=xc(:,ii)**3-m3(ii,ii,ii)-3.0_dp*xc(:,ii)*cov(ii,ii)
    end do
    avg=sum(if_diag,dim=2)/real(pp,dp);ans=0.0_dp
    do ii=1,pp;do jj=ii,pp;do kk=jj,pp
      ifs=xc(:,ii)*xc(:,jj)*xc(:,kk)-m3(ii,jj,kk)-xc(:,ii)*cov(jj,kk)- &
        xc(:,jj)*cov(ii,kk)-xc(:,kk)*cov(ii,jj)
      m=mult3(ii,jj,kk);ans(1)=ans(1)+real(m,dp)*sum(ifs*ifs)/real(nn*nn,dp)
      if(ii==kk)then
        ans(3)=ans(3)+sum(ifs*if_diag(:,ii))/real(nn*nn,dp)
        ans(2)=ans(2)+sum(ifs*avg)/real(nn*nn,dp)
      end if
    end do;end do;end do
  end subroutine direct_vm3

  subroutine direct_vm4(xc,cov,ans)
    real(dp),intent(in)::xc(:,:),cov(:,:)
    real(dp),intent(out)::ans(3)
    real(dp),allocatable::m3(:,:,:),m4(:,:,:,:),if4diag(:,:),if2(:,:),avg4(:),avg22(:),ifs(:),ift(:)
    integer::ii,jj,kk,ll,pp,nn,m
    nn=size(xc,1);pp=size(xc,2)
    allocate(m3(pp,pp,pp),m4(pp,pp,pp,pp),if4diag(nn,pp),if2(nn,pp),avg4(nn),avg22(nn),ifs(nn),ift(nn))
    do ii=1,pp;do jj=1,pp;do kk=1,pp
      m3(ii,jj,kk)=sum(xc(:,ii)*xc(:,jj)*xc(:,kk))/real(nn,dp)
      do ll=1,pp;m4(ii,jj,kk,ll)=sum(xc(:,ii)*xc(:,jj)*xc(:,kk)*xc(:,ll))/real(nn,dp);end do
    end do;end do;end do
    do ii=1,pp
      if2(:,ii)=xc(:,ii)**2-cov(ii,ii)
      if4diag(:,ii)=xc(:,ii)**4-m4(ii,ii,ii,ii)-4.0_dp*xc(:,ii)*m3(ii,ii,ii)
    end do
    avg4=sum(if4diag,dim=2)/real(pp,dp);avg22=0.0_dp
    do ii=1,pp;avg22=avg22+2.0_dp*cov(ii,ii)*if2(:,ii)/real(pp,dp);end do
    ans=0.0_dp
    do ii=1,pp;do jj=ii,pp;do kk=jj,pp;do ll=kk,pp
      ifs=xc(:,ii)*xc(:,jj)*xc(:,kk)*xc(:,ll)-m4(ii,jj,kk,ll)- &
        xc(:,ii)*m3(jj,kk,ll)-xc(:,jj)*m3(ii,kk,ll)-xc(:,kk)*m3(ii,jj,ll)-xc(:,ll)*m3(ii,jj,kk)
      m=mult4(ii,jj,kk,ll);ans(1)=ans(1)+real(m,dp)*sum(ifs*ifs)/real(nn*nn,dp)
      if(ii==ll)then
        ift=if4diag(:,ii)
        ans(3)=ans(3)+sum(ifs*ift)/real(nn*nn,dp)
        ans(2)=ans(2)+sum(ifs*avg4)/real(nn*nn,dp)
      else if(ii==jj .and. kk==ll .and. jj<kk)then
        ift=if2(:,ii)*cov(kk,kk)+cov(ii,ii)*if2(:,kk)
        ans(3)=ans(3)+real(m,dp)*sum(ifs*ift)/real(nn*nn,dp)
        ans(2)=ans(2)+real(m,dp)*sum(ifs*avg22)/real(nn*nn,dp)
      end if
    end do;end do;end do;end do
  end subroutine direct_vm4
end program test_finite_sample_corrections
