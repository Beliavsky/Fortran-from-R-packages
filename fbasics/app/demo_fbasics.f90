! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
program demo_fbasics
  use fbasics
  implicit none
  integer, parameter :: n=300
  real(dp) :: x(n), a(2,2)
  integer :: i, info
  real(dp), allocatable :: ainv(:,:)
  type(basic_stats_result) :: stats
  type(distribution_fit) :: nfit, tfit
  type(test_result) :: jb
  call set_lcg_seed(20260724_8)
  do i=1,n
    x(i)=0.001_dp+0.015_dp*rt_lcg(7.0_dp)
  end do
  stats=basic_stats(x)
  call fit_normal(x,nfit)
  call fit_student(x,tfit)
  jb=jarque_bera_test(x)
  a=reshape([4.0_dp,2.0_dp,2.0_dp,3.0_dp],[2,2])
  call matrix_inverse(a,ainv,info)
  write(*,'(a,i0)')'observations: ',stats%nobs
  write(*,'(a,es14.6)')'mean: ',stats%mean
  write(*,'(a,es14.6)')'standard deviation: ',stats%stdev
  write(*,'(a,es14.6)')'skewness: ',stats%skewness
  write(*,'(a,es14.6)')'excess kurtosis: ',stats%kurtosis
  write(*,'(a,2(1x,es14.6))')'normal parameters:',nfit%parameters
  write(*,'(a,3(1x,es14.6))')'Student parameters:',tfit%parameters
  write(*,'(a,es14.6,a,es14.6)')'Jarque-Bera: ',jb%statistic,' p-value: ',jb%p_value
  write(*,'(a,es14.6)')'expected zero-drift max drawdown, T=100: ',maxdd_expectation(0.0_dp,1.0_dp,100.0_dp)
  write(*,'(a,4(1x,f9.5))')'inverse matrix:',ainv
end program demo_fbasics
