! SPDX-License-Identifier: GPL-3.0-only
program test_distribution_glm
  use mass
  use test_support
  implicit none
  integer, parameter :: n=20
  real(dp) :: values(n), design(n,1), counts(n), mean_count
  integer, allocatable :: nb_draws(:)
  type(density_fit_result) :: distribution
  type(regression_result) :: poisson_fit, nb_fit
  type(loglinear_result) :: loglin
  integer :: i, status

  do i=1,n
    values(i)=2.0_dp+0.5_dp*sin(real(i,dp))+0.2_dp*cos(0.4_dp*real(i,dp))
    counts(i)=real(mod(i,5)+1,dp)
  end do
  call fit_distribution(values,'normal',distribution)
  call assert_true(distribution%status == mass_success, 'normal fit status')
  call assert_close(distribution%estimates(1),sum(values)/real(n,dp),1.0e-12_dp,'normal mean')
  call assert_true(distribution%estimates(2)>0.0_dp,'normal scale')

  design=1.0_dp
  call poisson_glm_fit(design,counts,poisson_fit)
  mean_count=sum(counts)/real(n,dp)
  call assert_true(poisson_fit%status == mass_success,'Poisson GLM status')
  call assert_close(exp(poisson_fit%coefficients(1)),mean_count,1.0e-8_dp,'Poisson intercept')

  call glm_nb_fit(design,counts,nb_fit,init_theta=5.0_dp)
  call assert_true(nb_fit%theta>0.0_dp,'negative-binomial theta positive')
  call assert_all_finite(nb_fit%coefficients,'negative-binomial coefficients finite')

  call loglinear_fit(design,counts,loglin)
  call assert_true(loglin%status == mass_success,'loglinear status')
  call assert_true(loglin%deviance>=0.0_dp,'loglinear deviance')

  call rnegbin(100,[2.5_dp],4.0_dp,nb_draws,seed=991,status=status)
  call assert_true(status == mass_success,'rnegbin status')
  call assert_true(all(nb_draws>=0.0_dp),'rnegbin support')
  write(*,'(a)') 'test_distribution_glm: PASS'
end program test_distribution_glm
