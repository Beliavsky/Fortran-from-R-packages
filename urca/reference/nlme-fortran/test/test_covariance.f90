! SPDX-License-Identifier: GPL-2.0-or-later
program test_covariance
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS
  use nlme_types
  use nlme_correlation, only : correlation_matrix, arma_autocorrelation
  use nlme_variance, only : variance_sd
  use nlme_pdmat, only : pd_matrix
  use test_support
  implicit none
  type(correlation_spec)::c
  type(variance_spec)::v
  type(pd_spec)::p
  real(dp),allocatable::m(:,:),sd(:),acf(:)
  real(dp)::t(4)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp]
  integer::status,g(4)=[1,1,2,2]

  c%kind=COR_AR1
  allocate(c%par(1))
  c%par=0.5_dp
  call correlation_matrix(c,t,m,status)
  call check(status==NLME_SUCCESS,'AR1 status')
  call check_close(m(1,4),0.125_dp,1.0e-12_dp,'AR1 lag 3')

  c%kind=COR_COMPOUND_SYMM
  c%par=0.2_dp
  call correlation_matrix(c,t,m,status)
  call check_close(m(1,4),0.2_dp,1.0e-12_dp,'compound symmetry')

  c%kind=COR_EXPONENTIAL
  c%par=2.0_dp
  c%nugget=.false.
  call correlation_matrix(c,t,m,status)
  call check_close(m(1,3),exp(-1.0_dp),1.0e-12_dp,'exponential correlation')

  call arma_autocorrelation([0.6_dp],[0.2_dp],4,acf,status)
  call check(status==NLME_SUCCESS,'ARMA ACF status')
  call check_close(acf(1),1.0_dp,1.0e-12_dp,'ARMA ACF zero')
  call check(acf(2)>0.6_dp .and. acf(2)<0.8_dp,'ARMA ACF lag one')

  v%kind=VAR_POWER
  allocate(v%par(1))
  v%par=0.5_dp
  call variance_sd(v,[1.0_dp,4.0_dp,9.0_dp,16.0_dp],g,sd,status)
  call check_vector_close(sd,[1.0_dp,2.0_dp,3.0_dp,4.0_dp],1.0e-12_dp,'varPower')

  p%kind=PD_LOG_CHOL
  p%dim=2
  allocate(p%par(3))
  p%par=[log(2.0_dp),log(1.5_dp),0.4_dp]
  call pd_matrix(p,m,status)
  call check(status==NLME_SUCCESS,'pdLogChol status')
  call check_close(m(1,1),4.0_dp,1.0e-12_dp,'pd diagonal 1')
  call check_close(m(2,1),0.8_dp,1.0e-12_dp,'pd off diagonal')
  call check_close(m(2,2),2.41_dp,1.0e-12_dp,'pd diagonal 2')
  write(*,'(a)')'test_covariance: PASS'
end program test_covariance
