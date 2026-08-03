! SPDX-License-Identifier: GPL-2.0-or-later
program test_cashflows
  use lifeinsurer
  implicit none
  type(insurance_tariff)::t
  type(cash_flow_set)::cf
  integer::st
  t%tariff_type=tariff_endowment; t%policy_period=20; t%premium_period=20
  t%premium_refund=0.0_dp
  call build_cash_flows(t,cf,st); call check(st==0,'endowment status')
  call check(maxval(abs(cf%premiums_advance(1:20)-1.0_dp))<1e-14_dp .and. &
    abs(cf%premiums_advance(21))<1e-14_dp,'endowment premiums')
  call check(abs(cf%survival_advance(21)-1.0_dp)<1e-14_dp,'endowment survival')
  call check(maxval(abs(cf%death_sum_insured(1:20)-1.0_dp))<1e-14_dp .and. abs(cf%death_sum_insured(21))<1e-14_dp,'endowment death')
  t%tariff_type=tariff_annuity; t%policy_period=55; t%premium_period=1
  call build_cash_flows(t,cf,st); call check(abs(sum(cf%premiums_advance)-1.0_dp)<1e-14_dp,'annuity single premium')
  call check(abs(sum(cf%survival_advance)-55.0_dp)<1e-14_dp,'annuity benefits')
  t%policy_period=80; t%premium_period=25; t%deferral_period=25; t%premium_refund=1.0_dp
  call build_cash_flows(t,cf,st)
  call check(maxval(abs(cf%premiums_advance(1:25)-1.0_dp))<1e-14_dp,'deferred premiums')
  call check(maxval(abs(cf%survival_advance(1:25)))<1e-14_dp .and. &
    abs(sum(cf%survival_advance)-55.0_dp)<1e-14_dp,'deferred annuity')
  call check(maxval(abs(cf%death_gross_premium(1:25)-[(real(st,dp),st=1,25)]))<1e-14_dp,'refund cumulative')
  call check(maxval(abs(cf%death_refund_past(1:25)-1.0_dp))<1e-14_dp .and. &
    abs(sum(cf%death_refund_past)-25.0_dp)<1e-14_dp,'refund flag')
  t%tariff_type=tariff_endowment_dread; t%policy_period=20; t%premium_period=20; t%deferral_period=0; t%premium_refund=0
  call build_cash_flows(t,cf,st)
  call check(abs(sum(cf%disease_sum_insured)-20.0_dp)<1e-14_dp,'dread disease')
  print '(a)','test_cashflows: PASS'
contains
  subroutine check(ok,msg); logical,intent(in)::ok; character(*),intent(in)::msg
    if(.not.ok) then; print '(a,1x,a)','FAIL:',msg; error stop 1; end if
  end subroutine
end program
