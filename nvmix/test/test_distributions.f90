! SPDX-License-Identifier: GPL-3.0-or-later
program test_distributions
  use nvmix
  implicit none
  real(dp) :: x(2),loc(2),scale(2,2),v
  type(probability_result) :: pr
  x=[0.2_dp,-0.1_dp]; loc=0.0_dp
  scale=reshape([1.0_dp,0.3_dp,0.3_dp,2.0_dp],[2,2])
  v=dnorm_mv(x,loc,scale)
  call assert_close(v,0.11212619979818034_dp,2.0e-13_dp,'normal density')
  v=dstudent_mv(x,7.0_dp,loc,scale)
  call assert_close(v,0.11128854053255031_dp,2.0e-13_dp,'Student density')
  call assert_close(normal_cdf(0.0_dp),0.5_dp,1.0e-15_dp,'normal cdf')
  call assert_close(student_cdf(student_quantile(0.93_dp,5.5_dp),5.5_dp),0.93_dp,2.0e-12_dp,'t inversion')
  call assert_close(gamma_cdf(gamma_quantile(0.8_dp,2.5_dp,1.2_dp),2.5_dp,1.2_dp),0.8_dp,2.0e-12_dp,'gamma inversion')
  pr=pnorm_mv([-huge(1.0_dp)],[1.0_dp],[0.0_dp],reshape([1.0_dp],[1,1]))
  call assert_close(pr%value,normal_cdf(1.0_dp),2.0e-14_dp,'normal probability')
  print '(a)','test_distributions: PASS'
contains
  subroutine assert_close(a,b,tol,label)
    real(dp), intent(in) :: a,b,tol
    character(*), intent(in) :: label
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      write(*,'(a,3es24.15)')trim(label)//' mismatch: ',a,b,abs(a-b)
      error stop 1
    end if
  end subroutine
end program
