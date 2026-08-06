! SPDX-License-Identifier: GPL-3.0-only
program test_probability
  use spantest_kinds, only : dp
  use spantest_probability, only : normal_cdf,normal_quantile,student_t_cdf,f_upper_tail
  implicit none
  call check(abs(normal_cdf(0.0_dp)-0.5_dp)<1.0e-15_dp,'normal CDF at zero')
  call check(abs(normal_quantile(0.975_dp)-1.95996398454005_dp)<5.0e-8_dp,'normal quantile')
  call check(abs(student_t_cdf(0.0_dp,10.0_dp)-0.5_dp)<1.0e-15_dp,'Student t CDF at zero')
  call check(abs(f_upper_tail(4.964602743730714_dp,1.0_dp,10.0_dp)-0.05_dp)<2.0e-8_dp,'F upper tail')
  print '(a)','test_probability: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if (.not.ok) error stop msg
  end subroutine
end program test_probability
