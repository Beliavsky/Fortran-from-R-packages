# API summary

The aggregate module is `lifeinsurer`.

## Main types

- `insurance_tariff`: contract type, periods, frequencies, interest, sum insured, state, costs, and rounding.
- `mortality_table`: annual death probabilities `qx` and optional disease/invalidity probabilities `ix`.
- `contract_result`: transitions, cash flows, present values, premiums, reserves, profits, and status.
- `profit_rate_table`: yearly guaranteed, interest, mortality, expense, sum, and terminal-bonus rates.

## Main procedures

- `calculate_contract`
- `build_transition_probabilities`
- `build_cash_flows`
- `calculate_present_values`
- `calculate_premiums`
- `calculate_sum_insured`
- `calculate_reserves`
- `calculate_profit_participation`
- `contract_grid_premium`
- `premium_waiver`
- `extend_contract`

## PVfactory equivalents

- `pv_guaranteed`
- `pv_survival`
- `pv_death`
- `pv_disease`
- `pv_after_death`

## Helper equivalents

- `correction_payment_frequency`
- `death_benefit_linear_decreasing`
- `death_benefit_annuity_decreasing`
- `premium_refund_period_default`
- `frequency_charge`
- `round_value`
- `get_savings_premium`
- `initialize_costs`

R names containing dots are written with underscores.
