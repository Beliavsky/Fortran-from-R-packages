! SPDX-License-Identifier: GPL-2.0-or-later
module lifeinsurer_types
  use lifeinsurer_kinds, only : dp, tariff_endowment, payment_advance
  implicit none
  private

  type, public :: civil_date
    integer :: year = 1970
    integer :: month = 1
    integer :: day = 1
  end type civil_date

  type, public :: mortality_table
    real(dp), allocatable :: qx(:)
    real(dp), allocatable :: ix(:)
  end type mortality_table

  type, public :: frequency_correction
    real(dp) :: alpha = 1.0_dp
    real(dp) :: beta = 0.0_dp
  end type frequency_correction

  type, public :: expense_loadings
    real(dp) :: alpha = 0.0_dp
    real(dp) :: zillmer = 0.0_dp
    real(dp) :: beta = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: gamma_premium_free = 0.0_dp
    real(dp) :: gamma_contract = 0.0_dp
    real(dp) :: gamma_after_death = 0.0_dp
    real(dp) :: unit_cost = 0.0_dp
    real(dp) :: unit_cost_policy = 0.0_dp
    real(dp) :: security = 0.0_dp
    real(dp) :: tax = 0.0_dp
    real(dp) :: ongoing_alpha_gross = 0.0_dp
    real(dp) :: premium_rebate = 0.0_dp
    real(dp) :: partner_rebate = 0.0_dp
    real(dp) :: frequency_loading = 0.0_dp
    integer :: commission_period = 1
  end type expense_loadings

  type, public :: rounding_rules
    integer :: premium_gross = 8
    integer :: premium_net = 8
    integer :: premium_written = 8
    integer :: sum_insured = 8
    integer :: reserve = 8
    integer :: surrender = 8
  end type rounding_rules

  type, public :: insurance_tariff
    integer :: tariff_type = tariff_endowment
    integer :: policy_period = 20
    integer :: deferral_period = 0
    integer :: premium_period = 20
    integer :: guaranteed_period = 0
    integer :: premium_refund_period = -1
    integer :: premium_frequency = 1
    integer :: benefit_frequency = 1
    integer :: premium_payment_time = payment_advance
    integer :: benefit_payment_time = payment_advance
    real(dp) :: interest = 0.0_dp
    real(dp) :: sum_insured = 1.0_dp
    real(dp) :: initial_capital = 0.0_dp
    real(dp) :: premium_refund = 0.0_dp
    real(dp) :: premium_increase = 1.0_dp
    real(dp) :: annuity_increase = 1.0_dp
    real(dp) :: surrender_factor = 1.0_dp
    logical :: invalidity_ends_contract = .true.
    logical :: premium_waiver = .false.
    type(expense_loadings) :: costs
    type(rounding_rules) :: rounding
  end type insurance_tariff

  type, public :: transition_probabilities
    real(dp), allocatable :: qx(:)
    real(dp), allocatable :: ix(:)
    real(dp), allocatable :: px(:)
  end type transition_probabilities

  type, public :: cash_flow_set
    real(dp), allocatable :: premiums_advance(:)
    real(dp), allocatable :: premiums_arrears(:)
    real(dp), allocatable :: additional_capital(:)
    real(dp), allocatable :: guaranteed_advance(:)
    real(dp), allocatable :: guaranteed_arrears(:)
    real(dp), allocatable :: survival_advance(:)
    real(dp), allocatable :: survival_arrears(:)
    real(dp), allocatable :: death_sum_insured(:)
    real(dp), allocatable :: disease_sum_insured(:)
    real(dp), allocatable :: death_gross_premium(:)
    real(dp), allocatable :: death_refund_past(:)
    real(dp), allocatable :: death_premium_free(:)
  end type cash_flow_set

  type, public :: present_value_set
    real(dp), allocatable :: premiums(:)
    real(dp), allocatable :: additional_capital(:)
    real(dp), allocatable :: guaranteed(:)
    real(dp), allocatable :: survival(:)
    real(dp), allocatable :: death_sum_insured(:)
    real(dp), allocatable :: disease_sum_insured(:)
    real(dp), allocatable :: death_gross_premium(:)
    real(dp), allocatable :: death_refund_past(:)
    real(dp), allocatable :: benefits(:)
  end type present_value_set

  type, public :: premium_values
    real(dp) :: unit_net = 0.0_dp
    real(dp) :: unit_zillmer = 0.0_dp
    real(dp) :: unit_gross = 0.0_dp
    real(dp) :: net = 0.0_dp
    real(dp) :: zillmer = 0.0_dp
    real(dp) :: gross = 0.0_dp
    real(dp) :: unit_cost = 0.0_dp
    real(dp) :: written_yearly = 0.0_dp
    real(dp) :: written_before_tax = 0.0_dp
    real(dp) :: tax = 0.0_dp
    real(dp) :: written = 0.0_dp
  end type premium_values

  type, public :: reserve_values
    real(dp), allocatable :: net(:)
    real(dp), allocatable :: zillmer(:)
    real(dp), allocatable :: adequate(:)
    real(dp), allocatable :: gamma(:)
    real(dp), allocatable :: contractual(:)
    real(dp), allocatable :: conversion(:)
    real(dp), allocatable :: reduction(:)
    real(dp), allocatable :: surrender(:)
    real(dp), allocatable :: premium_free_sum_insured(:)
  end type reserve_values

  type, public :: profit_rate_table
    real(dp), allocatable :: guaranteed_interest(:)
    real(dp), allocatable :: interest_profit(:)
    real(dp), allocatable :: interest_profit2(:)
    real(dp), allocatable :: mortality_profit(:)
    real(dp), allocatable :: expense_profit(:)
    real(dp), allocatable :: sum_profit(:)
    real(dp), allocatable :: terminal_bonus(:)
    real(dp), allocatable :: terminal_bonus_fund(:)
  end type profit_rate_table

  type, public :: profit_values
    real(dp), allocatable :: interest(:)
    real(dp), allocatable :: risk(:)
    real(dp), allocatable :: expense(:)
    real(dp), allocatable :: sum_component(:)
    real(dp), allocatable :: interest_on_profit(:)
    real(dp), allocatable :: terminal_bonus(:)
    real(dp), allocatable :: terminal_bonus_fund(:)
    real(dp), allocatable :: total_assignment(:)
    real(dp), allocatable :: accumulated(:)
    real(dp), allocatable :: benefit(:)
  end type profit_values

  type, public :: contract_result
    type(insurance_tariff) :: tariff
    type(transition_probabilities) :: transition
    type(cash_flow_set) :: cash_flows
    type(present_value_set) :: present_values
    type(premium_values) :: premiums
    type(reserve_values) :: reserves
    type(profit_values) :: profits
    integer :: status = 0
    character(len=160) :: message = ''
  end type contract_result
end module lifeinsurer_types
