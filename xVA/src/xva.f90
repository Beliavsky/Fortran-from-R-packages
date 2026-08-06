module xva
  use trading, only : dp, str_len, trade_t, csa_t, collateral_t, curve_t, &
    load_curve_csv
  use saccr, only : exposure_result_t
  use xva_types
  use xva_math, only : normal_pdf, inverse_normal_cdf
  use xva_core
  use xva_supervisory
  use xva_exposure
  use xva_regulatory
  use xva_calculator_mod
  implicit none
  public
end module xva
