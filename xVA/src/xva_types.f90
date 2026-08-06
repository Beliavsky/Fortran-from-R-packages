module xva_types
  use trading, only : dp, str_len
  use saccr, only : exposure_result_t
  implicit none
  private

  type, public :: simulation_data_t
    real(dp) :: pfe_percentile = 0.90_dp
    integer :: num_simulations = 1000
    real(dp) :: mean_reversion_a = 0.001_dp
    real(dp) :: volatility = 0.01_dp
  end type simulation_data_t

  type, public :: regulatory_data_t
    character(len=str_len) :: ccr_framework = "SA-CCR"
    character(len=str_len) :: sa_ccr_simplified = ""
    character(len=str_len) :: cva_framework = "BA-CVA"
    real(dp) :: sa_cva_multiplier = 1.0_dp
    character(len=str_len) :: counterparty_sector = &
      "Sovereigns including central banks and multilateral development banks"
    logical :: ignore_default_charge = .true.
    real(dp) :: pd_counterparty = 0.002_dp
    real(dp) :: pd_processing_organization = 0.005_dp
    real(dp) :: pd_funding = 0.001_dp
    real(dp) :: lgd = 0.45_dp
    real(dp) :: return_on_capital = 0.03_dp
    character(len=str_len) :: counterparty_rating = "A"
    real(dp) :: mva_days = 10.0_dp
    real(dp) :: mva_percentile = 0.99_dp
  end type regulatory_data_t

  type, public :: exposure_profile_t
    real(dp), allocatable :: ee_uncollateralized(:)
    real(dp), allocatable :: nee_uncollateralized(:)
    real(dp), allocatable :: pfe_uncollateralized(:)
    real(dp), allocatable :: ee(:)
    real(dp), allocatable :: nee(:)
    real(dp), allocatable :: pfe(:)
    real(dp), allocatable :: eee(:)
  end type exposure_profile_t

  type, public :: ead_result_t
    real(dp) :: ead_value = 0.0_dp
    type(exposure_result_t) :: exposure
    logical :: has_exposure = .false.
  end type ead_result_t

  type, public :: cva_sensitivity_t
    real(dp), allocatable :: credit_spread_delta(:)
    real(dp), allocatable :: credit_spread_tenors(:)
    real(dp), allocatable :: interest_rate_delta(:)
    real(dp), allocatable :: interest_rate_tenors(:)
  end type cva_sensitivity_t

  type, public :: cva_capital_result_t
    real(dp) :: interest_rate_charge = 0.0_dp
    real(dp) :: credit_spread_charge = 0.0_dp
    real(dp) :: total_charge = 0.0_dp
  end type cva_capital_result_t


  type, public :: labeled_matrix_t
    character(len=str_len), allocatable :: row_labels(:)
    character(len=str_len), allocatable :: column_labels(:)
    real(dp), allocatable :: values(:,:)
  end type labeled_matrix_t

  type, public :: relationship_correlation_t
    character(len=str_len) :: relationship = ""
    real(dp) :: correlation = 0.0_dp
  end type relationship_correlation_t

  type, public :: commodity_risk_weight_t
    integer :: bucket_number = 0
    character(len=str_len) :: group = ""
    real(dp) :: risk_weight = 0.0_dp
  end type commodity_risk_weight_t

  type, public :: equity_risk_weight_t
    integer :: bucket_number = 0
    character(len=str_len) :: size_class = ""
    character(len=str_len) :: region = ""
    character(len=str_len) :: sector = ""
    real(dp) :: risk_weight = 0.0_dp
  end type equity_risk_weight_t

  type, public :: supervisory_cva_data_t
    real(dp), allocatable :: ir_eligible_tenors(:)
    real(dp), allocatable :: ir_eligible_correlation(:,:)
    real(dp), allocatable :: ir_other_tenors(:)
    real(dp), allocatable :: ir_other_correlation(:,:)
    real(dp), allocatable :: ir_risk_weight_eligible(:)
    real(dp), allocatable :: ir_risk_weight_other(:)
    real(dp), allocatable :: cs_tenors(:)
    real(dp), allocatable :: cs_tenor_correlation(:,:)
    character(len=str_len), allocatable :: sectors(:)
    real(dp), allocatable :: sector_risk_weight_ig(:)
    real(dp), allocatable :: sector_risk_weight_hy_nr(:)
    character(len=str_len), allocatable :: ratings(:)
    real(dp), allocatable :: rating_weights(:)
    type(labeled_matrix_t) :: ir_eligible_full
    type(labeled_matrix_t) :: ir_other_full
    type(labeled_matrix_t) :: cs_correlation_by_sector
    type(labeled_matrix_t) :: cs_correlation_by_tenor_full
    type(labeled_matrix_t) :: cs_reference_sector_correlation
    type(labeled_matrix_t) :: cs_sector_counterparty_correlation
    type(relationship_correlation_t), allocatable :: hedge_counterparty_correlations(:)
    type(commodity_risk_weight_t), allocatable :: commodity_risk_weights(:)
    type(equity_risk_weight_t), allocatable :: equity_risk_weights(:)
  end type supervisory_cva_data_t

  type, public :: xva_result_t
    real(dp) :: cva_simulated = 0.0_dp
    real(dp) :: dva_simulated = 0.0_dp
    real(dp) :: fca_simulated = 0.0_dp
    real(dp) :: fba_simulated = 0.0_dp
    real(dp) :: mva_simulated = 0.0_dp
    logical :: has_simulated_values = .false.
    real(dp) :: cva_saccr = 0.0_dp
    real(dp) :: dva_saccr = 0.0_dp
    real(dp) :: fca_saccr = 0.0_dp
    real(dp) :: fba_saccr = 0.0_dp
    real(dp) :: mva_saccr = 0.0_dp
    logical :: has_saccr_values = .false.
    type(cva_capital_result_t) :: cva_capital
    real(dp) :: kva = 0.0_dp
    real(dp) :: ead = 0.0_dp
    real(dp) :: effective_maturity = 0.0_dp
    type(exposure_profile_t) :: exposure_profile
  end type xva_result_t

end module xva_types
