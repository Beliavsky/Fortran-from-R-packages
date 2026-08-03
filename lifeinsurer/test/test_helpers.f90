! SPDX-License-Identifier: GPL-2.0-or-later
program test_helpers
  use lifeinsurer
  implicit none
  real(dp),allocatable::b(:),x(:)
  type(frequency_correction)::c
  call check(is_single_premium_contract(1),'single premium')
  call check(is_regular_premium_contract(20),'regular premium')
  call check(premium_refund_period_default(20,0)==20,'refund whole period')
  call check(premium_refund_period_default(80,25)==25,'refund deferral')
  b=death_benefit_linear_decreasing(20,0,21)
  call close(b(1),1.0_dp,1e-14_dp,'linear first')
  call close(b(21),0.0_dp,1e-14_dp,'linear last')
  b=death_benefit_annuity_decreasing(0.0_dp,20,0,21)
  call close(b(11),0.5_dp,1e-14_dp,'annuity zero rate')
  c=correction_payment_frequency(0.03_dp,12,0.0_dp)
  call close(c%beta,11.0_dp/24.0_dp,1e-14_dp,'frequency correction')
  call close(frequency_charge(12,3.0_dp,2.0_dp,1.0_dp,0.0_dp),0.03_dp,1e-14_dp,'frequency charge')
  x=[2.0e-8_dp,2.0e-2_dp,2.0_dp,987654321.987654321_dp]
  x=round_value(x,2)
  call close(x(1),0.0_dp,1e-14_dp,'round small')
  call close(x(2),0.02_dp,1e-14_dp,'round cents')
  call close(x(4),987654321.99_dp,1e-8_dp,'round large')
  print '(a)','test_helpers: PASS'
contains
  subroutine check(ok,msg); logical,intent(in)::ok; character(*),intent(in)::msg
    if(.not.ok) then; print '(a,1x,a)','FAIL:',msg; error stop 1; end if
  end subroutine
  subroutine close(a,b,tol,msg); real(dp),intent(in)::a,b,tol; character(*),intent(in)::msg
    call check(abs(a-b)<=tol*(1.0_dp+abs(b)),msg)
  end subroutine
end program
