! SPDX-License-Identifier: GPL-3.0-only
program test_discriminant_robust
  use mass
  use test_support
  implicit none
  integer, parameter :: n=40
  real(dp) :: x(n,2)
  integer :: group(n), i, status
  integer, allocatable :: predicted(:)
  type(lda_model) :: lda
  type(qda_model) :: qda
  type(covariance_result) :: robust_cov

  do i=1,n/2
    x(i,1)=-2.0_dp+0.15_dp*sin(real(i,dp))
    x(i,2)=-1.0_dp+0.15_dp*cos(real(i,dp))
    group(i)=1
  end do
  do i=n/2+1,n
    x(i,1)=2.0_dp+0.15_dp*sin(real(i,dp))
    x(i,2)=1.0_dp+0.15_dp*cos(real(i,dp))
    group(i)=2
  end do
  call lda_fit(x,group,lda)
  call lda_predict(lda,x,predicted,status=status)
  call assert_true(status==mass_success,'LDA prediction status')
  call assert_true(count(predicted==group)>=38,'LDA classification')

  call qda_fit(x,group,qda)
  call qda_predict(qda,x,predicted,status=status)
  call assert_true(status==mass_success,'QDA prediction status')
  call assert_true(count(predicted==group)>=38,'QDA classification')

  x(n,1)=25.0_dp
  call cov_trob(x,robust_cov,nu=5.0_dp)
  call assert_true(robust_cov%status==mass_success,'cov.trob status')
  call assert_all_finite(reshape(robust_cov%covariance,[size(robust_cov%covariance)]), &
    'cov.trob finite')
  call cov_mcd(x,robust_cov,nsamp=200,seed=2026)
  call assert_true(robust_cov%status==mass_success,'cov.mcd status')
  call assert_true(all(diagonal(robust_cov%covariance)>0.0_dp),'cov.mcd positive diagonal')
  write(*,'(a)') 'test_discriminant_robust: PASS'
contains
  function diagonal(a) result(d)
    real(dp),intent(in)::a(:,:)
    real(dp)::d(min(size(a,1),size(a,2)))
    integer::j
    do j=1,size(d);d(j)=a(j,j);end do
  end function diagonal
end program test_discriminant_robust
