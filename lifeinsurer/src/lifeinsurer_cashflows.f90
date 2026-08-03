! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_cashflows
  use lifeinsurer_kinds, only : dp, lir_success, lir_invalid_argument, &
    tariff_annuity, tariff_term_fix, tariff_dread_disease, tariff_endowment, &
    tariff_pure_endowment, tariff_whole_life, tariff_endowment_dread, payment_advance
  use lifeinsurer_types, only : insurance_tariff, mortality_table, transition_probabilities, cash_flow_set
  use lifeinsurer_helpers, only : premium_refund_period_default
  implicit none
  private
  public :: build_transition_probabilities, build_cash_flows
contains
  subroutine alloc_cf(cf,n)
    type(cash_flow_set), intent(out) :: cf
    integer, intent(in) :: n
    allocate(cf%premiums_advance(n),cf%premiums_arrears(n),cf%additional_capital(n))
    allocate(cf%guaranteed_advance(n),cf%guaranteed_arrears(n))
    allocate(cf%survival_advance(n),cf%survival_arrears(n))
    allocate(cf%death_sum_insured(n),cf%disease_sum_insured(n))
    allocate(cf%death_gross_premium(n),cf%death_refund_past(n),cf%death_premium_free(n))
    cf%premiums_advance=0.0_dp; cf%premiums_arrears=0.0_dp; cf%additional_capital=0.0_dp
    cf%guaranteed_advance=0.0_dp; cf%guaranteed_arrears=0.0_dp
    cf%survival_advance=0.0_dp; cf%survival_arrears=0.0_dp
    cf%death_sum_insured=0.0_dp; cf%disease_sum_insured=0.0_dp
    cf%death_gross_premium=0.0_dp; cf%death_refund_past=0.0_dp; cf%death_premium_free=0.0_dp
  end subroutine

  subroutine build_transition_probabilities(tab,policy_period,invalidity_ends,tr,status)
    type(mortality_table), intent(in) :: tab
    integer, intent(in) :: policy_period
    logical, intent(in) :: invalidity_ends
    type(transition_probabilities), intent(out) :: tr
    integer, intent(out) :: status
    integer :: i,n
    real(dp) :: q,z
    if(policy_period<1 .or. .not.allocated(tab%qx) .or. size(tab%qx)<policy_period) then
      status=lir_invalid_argument; return
    end if
    n=policy_period
    allocate(tr%qx(n),tr%ix(n),tr%px(n)); tr%ix=0.0_dp
    do i=1,n
      q=max(0.0_dp,min(1.0_dp,tab%qx(i))); z=0.0_dp
      if(allocated(tab%ix)) then
        if(i<=size(tab%ix)) z=max(0.0_dp,min(1.0_dp,tab%ix(i)))
      end if
      if(invalidity_ends) then
        if(q+z>1.0_dp) z=1.0_dp-q
        tr%px(i)=1.0_dp-q-z
      else
        tr%px(i)=1.0_dp-q
      end if
      tr%qx(i)=q; tr%ix(i)=z
    end do
    status=lir_success
  end subroutine

  subroutine build_cash_flows(tariff,cf,status)
    type(insurance_tariff), intent(in) :: tariff
    type(cash_flow_set), intent(out) :: cf
    integer, intent(out) :: status
    integer :: n,i,prem,def,prot,guar,refund
    real(dp) :: cum
    if(tariff%policy_period<1 .or. tariff%premium_period<0 .or. &
       tariff%deferral_period<0 .or. tariff%deferral_period>tariff%policy_period) then
      status=lir_invalid_argument; return
    end if
    n=tariff%policy_period+1; call alloc_cf(cf,n)
    cf%additional_capital(1)=tariff%initial_capital
    prem=min(tariff%premium_period,tariff%policy_period)
    if(.not.tariff%premium_waiver) then
      do i=1,prem
        if(tariff%premium_payment_time==payment_advance) then
          cf%premiums_advance(i)=tariff%premium_increase**real(i-1,dp)
        else
          cf%premiums_arrears(i)=tariff%premium_increase**real(i-1,dp)
        end if
      end do
    end if
    def=tariff%deferral_period; prot=tariff%policy_period-def; guar=min(tariff%guaranteed_period,prot)
    select case(tariff%tariff_type)
    case(tariff_annuity)
      do i=def+1,tariff%policy_period
        if(i<=def+guar) then
          if(tariff%benefit_payment_time==payment_advance) then
            cf%guaranteed_advance(i)=tariff%annuity_increase**real(i-def-1,dp)
          else
            cf%guaranteed_arrears(i)=tariff%annuity_increase**real(i-def-1,dp)
          end if
        else
          if(tariff%benefit_payment_time==payment_advance) then
            cf%survival_advance(i)=tariff%annuity_increase**real(i-def-1,dp)
          else
            cf%survival_arrears(i)=tariff%annuity_increase**real(i-def-1,dp)
          end if
        end if
      end do
    case(tariff_term_fix)
      if(tariff%benefit_payment_time==payment_advance) then
        cf%guaranteed_advance(n)=1.0_dp
      else
        cf%guaranteed_arrears(n)=1.0_dp
      end if
    case(tariff_dread_disease)
      if(def+1<=n-1) cf%disease_sum_insured(def+1:n-1)=1.0_dp
    case(tariff_endowment,tariff_endowment_dread)
      if(tariff%benefit_payment_time==payment_advance) then
        cf%survival_advance(n)=1.0_dp
      else
        cf%survival_arrears(n)=1.0_dp
      end if
      if(def+1<=n-1) cf%death_sum_insured(def+1:n-1)=1.0_dp
      if(tariff%tariff_type==tariff_endowment_dread .and. def+1<=n-1) &
        cf%disease_sum_insured(def+1:n-1)=1.0_dp
    case(tariff_pure_endowment)
      if(tariff%benefit_payment_time==payment_advance) then
        cf%survival_advance(n)=1.0_dp
      else
        cf%survival_arrears(n)=1.0_dp
      end if
    case(tariff_whole_life)
      if(def+1<=n-1) cf%death_sum_insured(def+1:n-1)=1.0_dp
    case default
      status=lir_invalid_argument; return
    end select
    cf%death_premium_free=cf%death_sum_insured
    if(abs(tariff%premium_refund)>epsilon(1.0_dp)) then
      refund=tariff%premium_refund_period
      if(refund<0) refund=premium_refund_period_default(tariff%policy_period,tariff%deferral_period)
      refund=min(refund,tariff%policy_period); cum=0.0_dp
      do i=1,refund
        cum=cum+cf%premiums_advance(i)
        if(i>1) cum=cum+cf%premiums_arrears(i-1)
        cf%death_gross_premium(i)=cum
        if(cum>0.0_dp) cf%death_refund_past(i)=1.0_dp
      end do
    end if
    status=lir_success
  end subroutine
end module lifeinsurer_cashflows
