! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_actuarial
  use lifeinsurer_kinds, only : dp, lir_success, lir_invalid_argument, lir_zero_denominator
  use lifeinsurer_types, only : insurance_tariff, present_value_set, premium_values, reserve_values
  use lifeinsurer_helpers, only : round_value
  implicit none
  private
  public :: calculate_premiums, calculate_sum_insured, calculate_reserves
  public :: balance_sheet_reserve, premium_decomposition
contains
  subroutine calculate_premiums(t,pv,p,status)
    type(insurance_tariff), intent(in) :: t
    type(present_value_set), intent(in) :: pv
    type(premium_values), intent(out) :: p
    integer, intent(out) :: status
    real(dp) :: prem_pv,benefit_pv,refund_pv,den_net,den_gross,cost_si,cost_ann
    real(dp) :: before_tax
    if(.not.allocated(pv%premiums) .or. size(pv%premiums)==0) then
      status=lir_invalid_argument; return
    end if
    prem_pv=pv%premiums(1); benefit_pv=pv%benefits(1)-pv%additional_capital(1)/max(t%sum_insured,tiny(1.0_dp))
    refund_pv=pv%death_gross_premium(1)*t%premium_refund
    den_net=prem_pv-(1.0_dp+t%costs%security)*refund_pv
    den_gross=den_net-t%costs%beta*prem_pv
    cost_si=t%costs%alpha+t%costs%gamma_contract*max(pv%survival(1),0.0_dp)
    cost_ann=t%costs%unit_cost/max(t%sum_insured,tiny(1.0_dp)) &
      +t%costs%unit_cost_policy*max(pv%premiums(1),0.0_dp) &
      /max(t%sum_insured,tiny(1.0_dp))
    if(abs(den_net)<=tiny(1.0_dp) .or. abs(den_gross)<=tiny(1.0_dp)) then
      status=lir_zero_denominator; return
    end if
    p%unit_net=(1.0_dp+t%costs%security)*benefit_pv/den_net
    p%unit_zillmer=((1.0_dp+t%costs%security)*benefit_pv+t%costs%zillmer)/den_net
    p%unit_gross=((1.0_dp+t%costs%security)*benefit_pv+cost_si+cost_ann)/den_gross
    p%unit_gross=p%unit_gross*(1.0_dp+t%costs%ongoing_alpha_gross)
    p%net=round_value(p%unit_net*t%sum_insured,t%rounding%premium_net)
    p%zillmer=round_value(p%unit_zillmer*t%sum_insured,t%rounding%premium_net)
    p%gross=round_value(p%unit_gross*t%sum_insured,t%rounding%premium_gross)
    p%unit_cost=t%costs%unit_cost+t%costs%unit_cost_policy
    before_tax=(p%gross+p%unit_cost)*(1.0_dp+t%costs%frequency_loading)
    before_tax=before_tax*(1.0_dp-t%costs%premium_rebate-t%costs%partner_rebate)
    p%written_yearly=before_tax*real(max(1,t%premium_frequency),dp)
    p%written_before_tax=round_value(before_tax/real(max(1,t%premium_frequency),dp),t%rounding%premium_written)
    p%tax=round_value(p%written_before_tax*t%costs%tax,t%rounding%premium_written)
    p%written=round_value(p%written_before_tax+p%tax,t%rounding%premium_written)
    status=lir_success
  end subroutine

  subroutine calculate_sum_insured(t,pv,premium,premium_kind,sum_insured,status)
    type(insurance_tariff), intent(in) :: t
    type(present_value_set), intent(in) :: pv
    real(dp), intent(in) :: premium
    integer, intent(in) :: premium_kind
    real(dp), intent(out) :: sum_insured
    integer, intent(out) :: status
    type(insurance_tariff) :: u
    type(premium_values) :: p
    real(dp) :: unit
    u=t; u%sum_insured=1.0_dp
    call calculate_premiums(u,pv,p,status)
    if(status/=lir_success) return
    select case(premium_kind)
    case(1); unit=p%unit_net
    case(2); unit=p%unit_zillmer
    case(3); unit=p%unit_gross
    case default; unit=p%written
    end select
    if(abs(unit)<=tiny(1.0_dp)) then; status=lir_zero_denominator; return; end if
    sum_insured=round_value(premium/unit,t%rounding%sum_insured)
  end subroutine

  subroutine calculate_reserves(t,pv,p,res,status)
    type(insurance_tariff), intent(in) :: t
    type(present_value_set), intent(in) :: pv
    type(premium_values), intent(in) :: p
    type(reserve_values), intent(out) :: res
    integer, intent(out) :: status
    integer :: n
    real(dp) :: sf
    if(.not.allocated(pv%benefits)) then; status=lir_invalid_argument; return; end if
    n=size(pv%benefits); sf=1.0_dp+t%costs%security
    allocate(res%net(n),res%zillmer(n),res%adequate(n),res%gamma(n),res%contractual(n))
    allocate(res%conversion(n),res%reduction(n),res%surrender(n),res%premium_free_sum_insured(n))
    res%net=t%sum_insured*sf*pv%benefits-p%net*pv%premiums
    res%zillmer=t%sum_insured*sf*pv%benefits-p%zillmer*pv%premiums
    res%gamma=t%sum_insured*t%costs%gamma_contract*pv%survival
    res%adequate=t%sum_insured*sf*pv%benefits-p%gross*pv%premiums+res%gamma
    res%contractual=max(0.0_dp,res%zillmer+res%gamma)
    res%conversion=res%contractual
    res%reduction=max(0.0_dp,res%zillmer+t%sum_insured*t%costs%alpha &
      *max(0.0_dp,1.0_dp-pv%premiums/max(pv%premiums(1),tiny(1.0_dp))))
    res%surrender=round_value(max(0.0_dp,res%reduction*t%surrender_factor),t%rounding%surrender)
    where(pv%benefits>tiny(1.0_dp))
      res%premium_free_sum_insured=res%surrender/(sf*pv%benefits)*t%sum_insured/max(t%sum_insured,tiny(1.0_dp))
    elsewhere
      res%premium_free_sum_insured=0.0_dp
    end where
    res%net=round_value(res%net,t%rounding%reserve)
    res%zillmer=round_value(res%zillmer,t%rounding%reserve)
    res%adequate=round_value(res%adequate,t%rounding%reserve)
    res%gamma=round_value(res%gamma,t%rounding%reserve)
    res%contractual=round_value(res%contractual,t%rounding%reserve)
    res%reduction=round_value(res%reduction,t%rounding%reserve)
    status=lir_success
  end subroutine

  pure real(dp) function balance_sheet_reserve(v0,v1,factor,premium,unearned)
    real(dp), intent(in) :: v0,v1,factor,premium
    logical, intent(in), optional :: unearned
    logical :: u
    u=.false.; if(present(unearned)) u=unearned
    balance_sheet_reserve=(1.0_dp-factor)*(v0+merge(premium,0.0_dp,.not.u))+factor*v1
  end function

  subroutine premium_decomposition(p, savings, risk, costs)
    type(premium_values), intent(in) :: p
    real(dp), intent(in) :: savings,risk
    real(dp), intent(out) :: costs
    costs=p%gross-savings-risk
  end subroutine
end module lifeinsurer_actuarial
