! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_contract
  use lifeinsurer_kinds, only : dp, lir_success
  use lifeinsurer_types, only : insurance_tariff, mortality_table, contract_result, &
    frequency_correction, profit_rate_table
  use lifeinsurer_helpers, only : correction_payment_frequency
  use lifeinsurer_cashflows, only : build_transition_probabilities, build_cash_flows
  use lifeinsurer_pv, only : calculate_present_values
  use lifeinsurer_actuarial, only : calculate_premiums, calculate_reserves
  use lifeinsurer_profit, only : calculate_profit_participation
  implicit none
  private
  public :: calculate_contract, premium_waiver, extend_contract, contract_grid_premium
contains
  subroutine calculate_contract(tariff,mortality,result,rates)
    type(insurance_tariff),intent(in)::tariff
    type(mortality_table),intent(in)::mortality
    type(contract_result),intent(out)::result
    type(profit_rate_table),intent(in),optional::rates
    type(frequency_correction)::pc,bc
    real(dp)::v
    integer::st
    result%tariff=tariff
    call build_transition_probabilities(mortality,tariff%policy_period,tariff%invalidity_ends_contract,result%transition,st)
    if(st/=lir_success) then; result%status=st; result%message='invalid mortality table or contract period'; return; end if
    call build_cash_flows(tariff,result%cash_flows,st)
    if(st/=lir_success) then; result%status=st; result%message='invalid tariff cash-flow parameters'; return; end if
    pc=correction_payment_frequency(tariff%interest,max(1,tariff%premium_frequency),0.0_dp)
    bc=correction_payment_frequency(tariff%interest,max(1,tariff%benefit_frequency),0.0_dp)
    v=1.0_dp/(1.0_dp+tariff%interest)
    call calculate_present_values(result%cash_flows,result%transition,v,max(1,tariff%premium_frequency),pc,&
      max(1,tariff%benefit_frequency),bc,result%present_values,st)
    if(st/=lir_success) then; result%status=st; result%message='present-value calculation failed'; return; end if
    call calculate_premiums(tariff,result%present_values,result%premiums,st)
    if(st/=lir_success) then; result%status=st; result%message='premium calculation failed'; return; end if
    call calculate_reserves(tariff,result%present_values,result%premiums,result%reserves,st)
    if(st/=lir_success) then; result%status=st; result%message='reserve calculation failed'; return; end if
    if(present(rates)) then
      call calculate_profit_participation(rates,result%reserves,result%premiums,tariff%sum_insured,out=result%profits,status=st)
      if(st/=lir_success) then; result%status=st; result%message='profit participation failed'; return; end if
    end if
    result%status=lir_success; result%message='success'
  end subroutine

  pure function premium_waiver(tariff) result(out)
    type(insurance_tariff),intent(in)::tariff
    type(insurance_tariff)::out
    out=tariff; out%premium_waiver=.true.; out%premium_period=0
  end function

  pure function extend_contract(tariff,additional_period,new_sum_insured,new_interest) result(out)
    type(insurance_tariff),intent(in)::tariff
    integer,intent(in)::additional_period
    real(dp),intent(in),optional::new_sum_insured,new_interest
    type(insurance_tariff)::out
    out=tariff; out%policy_period=max(1,tariff%policy_period+additional_period)
    if(present(new_sum_insured)) out%sum_insured=new_sum_insured
    if(present(new_interest)) out%interest=new_interest
  end function

  subroutine contract_grid_premium(base,mortality,ages,periods,premium,status)
    type(insurance_tariff),intent(in)::base
    type(mortality_table),intent(in)::mortality
    integer,intent(in)::ages(:),periods(:)
    real(dp),allocatable,intent(out)::premium(:,:)
    integer,intent(out)::status
    type(insurance_tariff)::t
    type(mortality_table)::mt
    type(contract_result)::r
    integer::i,j,nmax
    allocate(premium(size(ages),size(periods))); premium=0.0_dp
    nmax=size(mortality%qx)
    do i=1,size(ages)
      do j=1,size(periods)
        t=base; t%policy_period=periods(j); t%premium_period=min(t%premium_period,t%policy_period)
        if(ages(i)+t%policy_period>nmax) then; status=1; return; end if
        allocate(mt%qx(t%policy_period)); mt%qx=mortality%qx(ages(i)+1:ages(i)+t%policy_period)
        if(allocated(mortality%ix)) then
          allocate(mt%ix(t%policy_period)); mt%ix=mortality%ix(ages(i)+1:ages(i)+t%policy_period)
        end if
        call calculate_contract(t,mt,r)
        if(r%status/=lir_success) then; status=r%status; return; end if
        premium(i,j)=r%premiums%written
        if(allocated(mt%qx)) deallocate(mt%qx)
        if(allocated(mt%ix)) deallocate(mt%ix)
      end do
    end do
    status=lir_success
  end subroutine
end module lifeinsurer_contract
