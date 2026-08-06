module saccr_types
  use trading, only : dp, str_len
  implicit none
  private

  type, public :: supervisory_record_t
    character(len=str_len) :: asset_class = ""
    character(len=str_len) :: subclass = ""
    real(dp) :: supervisory_factor = 0.0_dp
    real(dp) :: correlation = 0.0_dp
    real(dp) :: option_volatility = 0.0_dp
  end type supervisory_record_t

  type, public :: single_trade_addon_t
    character(len=str_len) :: external_id = ""
    character(len=str_len) :: asset_class = ""
    character(len=str_len) :: hedging_set = ""
    real(dp) :: adjusted_notional = 0.0_dp
    real(dp) :: maturity_factor = 0.0_dp
    real(dp) :: supervisory_delta = 0.0_dp
    real(dp) :: effective_notional = 0.0_dp
    real(dp) :: volatility = 0.0_dp
    logical :: has_volatility = .false.
  end type single_trade_addon_t

  type, public :: hedging_set_result_t
    character(len=str_len) :: asset_class = ""
    character(len=str_len) :: name = ""
    character(len=str_len) :: subgroup = ""
    real(dp) :: effective_notional = 0.0_dp
    real(dp) :: supervisory_factor = 0.0_dp
    real(dp) :: correlation = 0.0_dp
    real(dp) :: addon = 0.0_dp
  end type hedging_set_result_t

  type, public :: asset_class_result_t
    character(len=str_len) :: name = ""
    real(dp) :: systematic_component = 0.0_dp
    real(dp) :: idiosyncratic_component = 0.0_dp
    real(dp) :: addon = 0.0_dp
  end type asset_class_result_t

  type, public :: addon_result_t
    real(dp) :: addon = 0.0_dp
    real(dp) :: maturity_factor = 0.0_dp
    logical :: has_maturity_factor = .false.
    type(single_trade_addon_t), allocatable :: trades(:)
    type(hedging_set_result_t), allocatable :: hedging_sets(:)
    type(asset_class_result_t), allocatable :: asset_classes(:)
  end type addon_result_t

  type, public :: replacement_cost_t
    real(dp) :: v_c = 0.0_dp
    real(dp) :: rc = 0.0_dp
    real(dp) :: v = 0.0_dp
  end type replacement_cost_t

  type, public :: exposure_result_t
    character(len=str_len) :: counterparty = ""
    character(len=str_len) :: csa_id = ""
    logical :: margined = .false.
    real(dp) :: maturity_factor = 0.0_dp
    type(addon_result_t) :: addon
    type(replacement_cost_t) :: replacement_cost
    real(dp) :: pfe = 0.0_dp
    real(dp) :: ead = 0.0_dp
    real(dp) :: unmargined_ead = 0.0_dp
  end type exposure_result_t

  type, public :: portfolio_result_t
    type(exposure_result_t), allocatable :: exposures(:)
    real(dp) :: total_ead = 0.0_dp
  end type portfolio_result_t

  type, public :: fx_hedge_result_t
    real(dp) :: fx_rwa = 0.0_dp
    real(dp) :: cds_rwa = 0.0_dp
    real(dp) :: protected_amount = 0.0_dp
    real(dp) :: unprotected_amount = 0.0_dp
    real(dp) :: multiplier = 1.0_dp
  end type fx_hedge_result_t

end module saccr_types
